#!/usr/bin/env python3
"""Figures of the experiments beyond the main grid, one subfigure per PDF, as
grouped bars (plots/bars.py): a panel per platform, a group per phaser, a bar
per thing compared, its value written over it. Each figure draws what
is computed and leaves the rest out, so it can be run again as runs finish.

A  extra variants in the phaser input (bin/make_truth_variants.sh), 30x: the
   same phaser on the same cell with the GIAB truth het SNPs (pale bar) and
   with the het indels added (solid bar) and with every het variant (truthall,
   hatched); --inputs truthindel leaves truthall out.
   Phasing is scored on the SNPs only, so a change is what the added variants
   do to the phase of the SNPs.
     figA1_switch_<ref>          switch error rate, chr1-22
     figA2_hamming_<ref>         blockwise Hamming error rate
     figA3_ng50_<ref>            block NG50 over chr1-22
     figA4_reads_missed_<ref>    truth-tagged reads not tagged to their haplotype
     figA5_indels_phased_<ref>   input indels the phaser phased (bin/variants_eval.sh)
     figA6_indels_wrong_<ref>    oriented indels on the wrong haplotype
     figA7_reads_<ref>           reads tagged against the SNP input: gained, lost
B  methylation (meth_run.sh), truth input, 30x: each phaser alone (pale), with
   MethPhaser on top of its VCF (solid), and LongPhase with its own 5mC mode
   (hatched)
     figB1_switch_<ref>, figB2_hamming_<ref>, figB3_ng50_<ref>,
     figB4_read_recall_matpat_<ref>   read haplotagging recall: of the chr1-22
                                      reads the matpat truth labels (alignment to
                                      the Q100 parental copies), the share tagged
                                      to their haplotype; the reads of the run's
                                      own reference, no truth SNP minimum
                                      (READ_RECALL_COMMON=0 READ_RECALL_MIN_SNPS=0
                                      READ_TRUTH=matpat bin/read_recall_strata.sh)
C  assemblies (dipcall_run.sh) next to the read-based phasers on the truth
   input, one coverage (--asm-depth, 30x). An assembly brings its own SNPs
   (dipcall, inside dip.bed). hifiasm without parents is shown for the switch
   rate only: its blocks (a pair of contigs) are not phase-consistent, so its
   Hamming, NG50 and read accuracy say nothing about the assembly (D4).
     figC1_switch_<depth>x_<ref>, figC2_hamming_..., figC3_ng50_...,
     figC4_reads_wrong_...       tagged reads on the wrong haplotype
     figC5_reads_missed_...      truth-tagged reads not tagged to their haplotype

HG002 files carry the reference only; another sample's carry its name too
(figA1_switch_HG001_GRCh38).

  python plots/figures_experiments.py                        all there is
  python plots/figures_experiments.py --only A,B --png
  python plots/figures_experiments.py --inputs truthindel
  python plots/figures_experiments.py --only C --asm-depth 50
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import style as st
from bars import bar, group_ticks, label_bar, offsets, scale_axes, value_text
from data import Metrics, _read_table, num, ratio, sample_root
from figures import (NoData, finite, matpat_recall, n50_mb, panel_row, positive, read_precision_pct,
                     read_recall_pct, save, switch_pct, top_legend)

DEPTH = 30
ASM_DEPTHS = [20, 30, 50]

INPUTS_A = ["truth", "truthindel", "truthall"]
INPUT_LABEL = {"truth": "HET SNPs", "truthindel": "HET SNPs + HET Indels",
               "truthall": "HET SNPs + HET Indels + multi-allelic variants"}

ASSEMBLERS = ["hifiasm", "hifiasmtrio", "verkkotrio"]
ASM_LABEL = {"hifiasm": "hifiasm", "hifiasmtrio": "hifiasm trio", "verkkotrio": "Verkko trio"}
# one sage family, light to dark, apart from the phasers' and platforms' hues
ASM_COLOR = {"hifiasm": "#A9BC8F", "hifiasmtrio": "#5F8A5B", "verkkotrio": "#34573F"}

PALE = 0.35          # the reference bar of a group: the same colour, lighter


def hamming_pct(db, key):
    return 100 * db.get(key, "hamming_rate")


def reads_missed_pct(db, key):
    """truth-tagged reads the run does not tag to their haplotype: 1 - recall"""
    return 100 - read_recall_pct(db, key)


def reads_wrong_pct(db, key):
    """tagged reads on the wrong haplotype: 1 - precision"""
    return 100 - read_precision_pct(db, key)


# (file tag, axis label, value, log axis)
METRICS = [
    ("switch", "Switch error rate (%)", switch_pct, True),
    ("hamming", "Hamming error rate (%)", hamming_pct, True),
    ("ng50", "Block NG50 (Mb)", n50_mb, True),
    ("reads_missed", "Truth-tagged reads not tagged\nto their haplotype (%)", reads_missed_pct, False),
]
READS_WRONG = ("reads_wrong", "Tagged reads on the\nwrong haplotype (%)", reads_wrong_pct, False)


B_RECALL_SFX = "_matpat_own_min0"
_b_recall = {}


def read_recall_matpat_pct(db, key):
    """reads tagged to their haplotype / reads the matpat truth labels"""
    if db.sample_root not in _b_recall:
        _b_recall[db.sample_root] = matpat_recall(db, B_RECALL_SFX)
    return _b_recall[db.sample_root].get(key, float("nan"))


# B scores reads against the Q100 parental copies, not the truth VCF's tags
READ_RECALL_MATPAT = ("read_recall_matpat", "Read haplotagging recall (%)",
                      read_recall_matpat_pct, False)
METRICS_B = METRICS[:3] + [READ_RECALL_MATPAT]


def fig_name(prefix, sample, ref):
    return "%s_%s" % (prefix, ref) if sample == "HG002" else "%s_%s_%s" % (prefix, sample, ref)


def phasers_on(platform):
    return [p for p in st.PHASERS if st.applies(platform, p)]


def legend_patches(items):
    """[(label, alpha, hatch)] in a neutral grey"""
    return [Patch(facecolor=st.MUTED, alpha=alpha, hatch=hatch, edgecolor="white" if hatch else "none",
                  label=label) for label, alpha, hatch in items]


def group_panels(platforms, height):
    """A panel per platform, each as wide as its number of phaser groups, so
    a group is as wide on PacBio (five phasers) as on ONT (four); one shared
    y axis, the legend strip under the panels."""
    widths = [len(phasers_on(p)) for p in platforms]
    fig = plt.figure(figsize=(st.FULL_WIDTH * st.MM, height * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(2, len(platforms), width_ratios=widths, height_ratios=[1, 0.1])
    axes = []
    for i in range(len(platforms)):
        ax = fig.add_subplot(grid[0, i], sharey=axes[0] if axes else None)
        if i:
            ax.tick_params(labelleft=False)
        axes.append(ax)
    legend_ax = fig.add_subplot(grid[1, :])
    legend_ax.axis("off")
    return fig, legend_ax, axes


def compare_groups(db, ctx, name, ref, metric, platforms, bars_of, legend):
    """A panel per platform, a group per phaser; bars_of(platform, phaser) ->
    [(key, alpha, hatch)], the first the reference the others are compared with."""
    tag, ylabel, value, log = metric
    fig, legend_ax, axes = group_panels(platforms, st.PANEL_HEIGHT + 14)
    drawn = []
    for ax, platform in zip(axes, platforms):
        phasers = phasers_on(platform)
        for x, phaser in enumerate(phasers):
            group = bars_of(platform, phaser)
            xs, w = offsets(len(group))
            values = [value(db, k) for k, _, _ in group]
            values = [positive(v) for v in values] if log else values
            for dx, v, (_, alpha, hatch) in zip(xs, values, group):
                if bar(ax, x + dx, v, w, st.PHASER_COLOR[phaser], alpha, hatch) is None:
                    continue
                drawn.append(v)
                # the value alone, no change against the first bar of the group
                label_bar(ax, x + dx, v, w, value_text(v))
        group_ticks(ax, [st.PHASER_LABEL[p] for p in phasers])
        ax.set_title(st.PLATFORM_LABEL[platform])
    if not drawn:
        plt.close(fig)
        raise NoData("no %s values" % tag)
    scale_axes(axes, drawn, log, ylabel)
    if tag.startswith("read_recall"):
        # a recall stops at 100%; the room over it is for the bar labels
        axes[0].set_ylim(0, 112)
        axes[0].set_yticks(range(0, 101, 20))
    top_legend(legend_ax, legend_patches(legend))
    return save(fig, ctx, name)


# ------------------------------------------------------- A: extra variants

def figA_metric(db, ctx, sample, ref, number, metric, extra):
    keys = lambda p, ph: [(ref, c, p, DEPTH, ph) for c in ["truth"] + extra]
    platforms = [p for p in st.PLATFORMS
                 if any(db.has(k) for ph in phasers_on(p) for k in keys(p, ph)[1:])]
    if not platforms:
        raise NoData("no %s runs of %s on %s" % (" / ".join(extra), sample, ref))
    hatch = {"truthindel": None, "truthall": "////"}
    bars_of = lambda p, ph: [(k, PALE if k[1] == "truth" else 1.0, hatch.get(k[1]))
                             for k in keys(p, ph)]
    legend = [(INPUT_LABEL["truth"], PALE, None)] + [(INPUT_LABEL[c], 1.0, hatch[c]) for c in extra]
    return compare_groups(db, ctx, fig_name("figA%d_%s" % (number, metric[0]), sample, ref), ref,
                          metric, platforms, bars_of, legend)


def variants_table(root, sample):
    """bin/variants_eval.sh: (ref, caller, platform, depth, phaser) -> numbers"""
    out = {}
    path = os.path.join(sample_root(root, sample), "metrics", "summary_variants.tsv")
    for line in _read_table(path):
        try:
            key = (line["reference"], line["caller"], line["platform"],
                   int(float(line["depth"])), line["phaser"])
        except (KeyError, ValueError):
            continue
        out[key] = {k: num(v) for k, v in line.items()}
    return out


def variant_bars(ctx, sample, ref, name, table, extra, value, ylabel, note=None, signed=False):
    """A panel per platform, a group per phaser, a bar per extra input
    (truthindel solid, truthall hatched). value(row) -> [numbers]: one bar per
    number from the same x (signed: one up, one down)."""
    platforms = [p for p in st.PLATFORMS
                 if any((ref, c, p, DEPTH, ph) in table for c in extra for ph in phasers_on(p))]
    if not platforms:
        raise NoData("no variants_eval rows for %s %s" % (sample, ref))
    hatch = {"truthindel": None, "truthall": "////"}
    fig, legend_ax, axes = group_panels(platforms, st.PANEL_HEIGHT + 14)
    xs, w = offsets(len(extra))
    drawn = []
    for ax, platform in zip(axes, platforms):
        phasers = phasers_on(platform)
        for x, phaser in enumerate(phasers):
            for dx, c in zip(xs, extra):
                row = table.get((ref, c, platform, DEPTH, phaser))
                for v in (value(row) if row else []):
                    if bar(ax, x + dx, v, w, st.PHASER_COLOR[phaser], 1.0, hatch[c]) is not None:
                        drawn.append(v)
                        label_bar(ax, x + dx, v, w)
        if signed:
            ax.axhline(0, color=st.INK, linewidth=0.6, zorder=4)
        group_ticks(ax, [st.PHASER_LABEL[p] for p in phasers])
        ax.set_title(st.PLATFORM_LABEL[platform])
    if not drawn:
        plt.close(fig)
        raise NoData("no values for %s %s %s" % (name, sample, ref))
    scale_axes(axes, drawn, False, ylabel)
    # a note that says how to read the bars (A6, A7) heads the legend under the panels
    top_legend(legend_ax, legend_patches([(INPUT_LABEL[c], 1.0, hatch[c]) for c in extra]), title=note)
    return save(fig, ctx, fig_name(name, sample, ref))


def figA5(ctx, sample, ref, table, extra):
    return variant_bars(ctx, sample, ref, "figA5_indels_phased", table, extra,
                        lambda r: [100 * ratio(r.get("indel_phased"), r.get("indel_in_input"))],
                        "Input indels phased (%)")


def figA6(ctx, sample, ref, table, extra):
    return variant_bars(ctx, sample, ref, "figA6_indels_wrong", table, extra,
                        lambda r: [100 * ratio(r.get("indel_wrong"), r.get("indel_oriented"))],
                        "Indels on the wrong haplotype (%)",
                        "of the indels in blocks that the truth SNPs orient")


def figA7(ctx, sample, ref, table, extra):
    """Reads tagged with the extra variants against the SNP input of the same
    phaser, as a share of the reads the SNP input tagged: gained above zero,
    lost below."""
    def gained_lost(r):
        base = r.get("reads_tagged_base")
        lost = 100 * ratio(r.get("reads_lost"), base)
        return [100 * ratio(r.get("reads_rescued"), base), -lost if finite(lost) else lost]
    return variant_bars(ctx, sample, ref, "figA7_reads", table, extra, gained_lost,
                        "Reads tagged against the SNP input\n(% of the reads it tagged)",
                        "up: reads gained, down: reads lost", signed=True)


# -------------------------------------------------------- B: methylation

def figB_metric(db, ctx, sample, ref, number, metric):
    keys = lambda p: [(ref, "truth", p, DEPTH, m) for m in st.METH_PHASERS if st.applies(p, m)]
    platforms = [p for p in st.PLATFORMS if any(db.has(k) for k in keys(p))]
    if not platforms:
        raise NoData("no methylation runs of %s on %s" % (sample, ref))

    def bars_of(platform, phaser):
        out = [((ref, "truth", platform, DEPTH, phaser), PALE, None),
               ((ref, "truth", platform, DEPTH, "mp" + phaser), 1.0, None)]
        if phaser == "longphase":
            out.append(((ref, "truth", platform, DEPTH, "lpmod"), 1.0, "////"))
        return out

    legend = [("phaser alone", PALE, None), ("+ MethPhaser", 1.0, None),
              ("LongPhase + 5mC (modcall)", 1.0, "////")]
    return compare_groups(db, ctx, fig_name("figB%d_%s" % (number, metric[0]), sample, ref), ref,
                          metric, platforms, bars_of, legend)


# ---------------------------------------------------------- C: assemblies

def figC_metric(db, ctx, sample, ref, number, metric, depth, with_hifiasm):
    """A panel per platform with an assembly: the read-based phasers on the
    truth input (their colours), then the assemblies (a sage family)."""
    tag, ylabel, value, log = metric
    assemblers = ASSEMBLERS if with_hifiasm else ASSEMBLERS[1:]

    def val(k):
        v = value(db, k)
        return positive(v) if log else v

    platforms = [p for p in st.PLATFORMS
                 if any(finite(val((ref, a, p, depth, "dipcall"))) for a in assemblers)]
    if not platforms:
        raise NoData("no assembly runs of %s at %dx on %s" % (sample, depth, ref))
    fig, legend_ax, axes = panel_row(len(platforms), height=st.PANEL_HEIGHT + 24)
    drawn = []
    for ax, platform in zip(axes, platforms):
        groups = [(st.PHASER_LABEL[p], (ref, "truth", platform, depth, p), st.PHASER_COLOR[p])
                  for p in phasers_on(platform)]
        n_phasers = len(groups)
        groups += [(ASM_LABEL[a], (ref, a, platform, depth, "dipcall"), ASM_COLOR[a])
                   for a in assemblers if db.has((ref, a, platform, depth, "dipcall"))]
        # the assemblies a little apart from the phasers, a dotted line between
        xs = [i + (0.6 if i >= n_phasers else 0) for i in range(len(groups))]
        # many bars: the values turned upright so neighbours do not overlap
        label_width = 1.0 if len(groups) <= 5 else 0.4
        for x, (_, k, color) in zip(xs, groups):
            v = val(k)
            if bar(ax, x, v, 1.0, color) is not None:
                drawn.append(v)
                label_bar(ax, x, v, label_width)
        if n_phasers < len(groups):
            ax.axvline(n_phasers - 0.2, color=st.MUTED, linestyle=":", linewidth=0.6, zorder=1)
        ax.set_xticks(xs)
        ax.set_xticklabels([label for label, _, _ in groups], rotation=60, ha="right",
                           rotation_mode="anchor")
        ax.set_xlim(xs[0] - 0.6, xs[-1] + 0.6)
        ax.set_title(st.PLATFORM_LABEL[platform])
    if not drawn:
        plt.close(fig)
        raise NoData("no %s values" % tag)
    scale_axes(axes, drawn, log, "%s, %d×" % (ylabel, depth))
    note = "read-based phasers: GIAB truth SNPs as input; assemblies: their own SNPs (dipcall)"
    if not with_hifiasm:
        note += "; hifiasm without parents: switch rate only"
    legend_ax.text(0.5, 0.5, note, ha="center", va="center", transform=legend_ax.transAxes)
    return save(fig, ctx, fig_name("figC%d_%s_%dx" % (number, tag, depth), sample, ref))


# --------------------------------------------------------------------- main

def samples_of(ctx):
    if ctx.samples:
        return ctx.samples.split()
    found = ["HG002"]
    sdir = os.path.join(ctx.root, "samples")
    if os.path.isdir(sdir):
        found += sorted(d for d in os.listdir(sdir) if os.path.isdir(os.path.join(sdir, d)))
    return found


def run(label, draw, *args):
    try:
        print("  ok      %-34s %s" % (label, draw(*args)))
        return "ok"
    except NoData as why:
        print("  skip    %-34s %s" % (label, why))
        return "skip"
    except Exception:
        import traceback
        print("  FAILED  %s" % label)
        traceback.print_exc()
        return "failed"
    finally:
        plt.close("all")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--out", default=None, help="default <root>/plots/pdf")
    parser.add_argument("--only", default="A,B,C", help="experiments, comma-separated: A, B, C")
    parser.add_argument("--samples", default="", help="default: HG002 and samples/*")
    parser.add_argument("--refs", default=" ".join(st.REFERENCES),
                        help="references of HG002; the other samples are GRCh38 only")
    parser.add_argument("--inputs", default="truthindel truthall",
                        help="experiment A inputs set against the truth SNPs: truthindel, truthall")
    parser.add_argument("--asm-depth", type=int, default=30, choices=ASM_DEPTHS)
    parser.add_argument("--rules", default="auto_min2")
    parser.add_argument("--png", action="store_true")
    ctx = parser.parse_args()
    ctx.out = ctx.out or os.path.join(ctx.root, "plots", "pdf")
    only = {s.strip().upper() for s in ctx.only.split(",") if s.strip()}
    extra = [c for c in INPUTS_A[1:] if c in ctx.inputs.split()]
    if not extra:
        parser.error("--inputs names none of truthindel, truthall")

    font = st.apply_style(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    if font != st.FONT_STACK[0]:
        print("note: %s not found, using %s" % (st.FONT_STACK[0], font))

    results = []
    for sample in samples_of(ctx):
        if not os.path.isdir(sample_root(ctx.root, sample)):
            continue
        refs = ctx.refs.split() if sample == "HG002" else [r for r in ctx.refs.split() if r == "GRCh38"]
        print("%s %s" % (sample, " ".join(refs)))
        if "A" in only or "B" in only:
            db = Metrics(ctx.root, sample, st.REGION, readacc_rules=ctx.rules).load(
                refs, ["truth"] + extra, st.PLATFORMS, [DEPTH], st.PHASERS + st.METH_PHASERS)
            table = variants_table(ctx.root, sample) if "A" in only else {}
            for ref in refs:
                if "A" in only:
                    for number, metric in enumerate(METRICS, 1):
                        results.append(run("A %s %s %s" % (sample, ref, metric[0]),
                                           figA_metric, db, ctx, sample, ref, number, metric, extra))
                    for label, draw in (("indels phased", figA5), ("indels wrong", figA6),
                                        ("reads", figA7)):
                        results.append(run("A %s %s %s" % (sample, ref, label),
                                           draw, ctx, sample, ref, table, extra))
                if "B" in only:
                    for number, metric in enumerate(METRICS_B, 1):
                        results.append(run("B %s %s %s" % (sample, ref, metric[0]),
                                           figB_metric, db, ctx, sample, ref, number, metric))
        if "C" in only:
            db = Metrics(ctx.root, sample, st.REGION, readacc_rules=ctx.rules).load(
                refs, ["truth"] + ASSEMBLERS, st.PLATFORMS, [ctx.asm_depth], st.PHASERS + ["dipcall"])
            c_metrics = [(METRICS[0], True), (METRICS[1], False), (METRICS[2], False),
                         (READS_WRONG, False), (METRICS[3], False)]
            for ref in refs:
                for number, (metric, with_hifiasm) in enumerate(c_metrics, 1):
                    results.append(run("C %s %s %s" % (sample, ref, metric[0]), figC_metric, db, ctx,
                                       sample, ref, number, metric, ctx.asm_depth, with_hifiasm))
    print("%d written, %d skipped for lack of data, %d failed -> %s" % (
        results.count("ok"), results.count("skip"), results.count("failed"), ctx.out))
    return 1 if "failed" in results else 0


if __name__ == "__main__":
    sys.exit(main())
