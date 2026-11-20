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
     figA4_read_recall_matpat_<ref>   read haplotagging recall, as figB4
                                      (the same command with
                                      PHASING_CALLERS="truthindel truthall")
     figA5_indels_phased_<ref>   input indels the phaser phased (bin/variants_eval.sh)
     figA6_indels_wrong_<ref>    oriented indels on the wrong haplotype
     figA7_reads_<ref>           reads tagged against the SNP input: gained, lost
B  methylation (meth_run.sh), truth input, 30x: each phaser alone (pale) and
   with MethPhaser on top of its VCF (solid)
     figB1_switch_<ref>, figB2_hamming_<ref>, figB3_ng50_<ref>,
     figB4_read_recall_matpat_<ref>   read haplotagging recall: of the chr1-22
                                      reads the matpat truth labels (alignment to
                                      the Q100 parental copies), the share tagged
                                      to their haplotype; the reads of the run's
                                      own reference, no truth SNP minimum
                                      (READ_RECALL_COMMON=0 READ_RECALL_MIN_SNPS=0
                                      READ_TRUTH=matpat bin/read_recall_strata.sh)
     figB5_longphase_methylation_<ref>   LongPhase alone, with MethPhaser and
                                      with its own 5mC mode (longphase modcall),
                                      a panel per metric of B1-B4, a group per
                                      platform
C  assemblies (dipcall_run.sh) next to the read-based phasers on the truth
   input, one coverage (--asm-depth, 30x). An assembly brings its own SNPs
   (dipcall, inside dip.bed). hifiasm without parents is shown for the switch
   rate only: its blocks (a pair of contigs) are not phase-consistent, so its
   Hamming, NG50 and read accuracy say nothing about the assembly (D4).
     figC1_switch_<depth>x_<ref>, figC2_hamming_..., figC3_ng50_...,
     figC4_reads_wrong_matpat_...    tagged reads on the wrong parental haplotype
     figC5_reads_missed_matpat_...   labelled reads not tagged to their parental
                                     haplotype (Q100 labels, Metrics.matpat)
     figC1_side_switch_ng50_...      the phasers on the truth SNPs and on the
                                     hifiasm SNPs, horizontal bars
     figC6_switch_ng50_...           the same as switch rate against NG50, as B6,
                                     on the SNPs phased in both runs
                                     (asm_common_hifiasm/summary_<depth>x.txt)

HG002 files carry the reference only; another sample's carry its name too
(figA1_switch_HG001_GRCh38).

  python plots/figures_experiments.py                        all there is
  python plots/figures_experiments.py --only A,B --png
  python plots/figures_experiments.py --inputs truthindel
  python plots/figures_experiments.py --only C --asm-depth 50
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator, NullLocator

import style as st
from bars import bar, group_ticks, label_bar, offsets, scale_axes, value_text
from data import READ_TRUTH, TRUTH_SFX, Metrics, _read_table, num, ratio, sample_root
from figures import (NoData, empty_panel, finite, log_ticks, matpat_precision_pct, matpat_recall, matpat_recall_pct,
                     n50_mb, panel_row, paper_title, positive, save, switch_pct, top_legend)

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
    """reads with a Q100 parental label the run does not tag to their
    haplotype: 1 - recall (Metrics.matpat)"""
    return 100 - matpat_recall_pct(db, key)


def reads_wrong_pct(db, key):
    """tagged reads on the wrong parental haplotype: 1 - precision (Metrics.matpat)"""
    return 100 - matpat_precision_pct(db, key)


# (file tag, axis label, value, log axis)
METRICS = [
    ("switch", "Switch error rate (%)", switch_pct, True),
    ("hamming", "Hamming error rate (%)", hamming_pct, True),
    ("ng50", "Phase-block NG50 (Mb)", n50_mb, True),
    ("reads_missed" + TRUTH_SFX, "Labelled reads not tagged to\ntheir parental haplotype (%)", reads_missed_pct, False),
]
READS_WRONG = ("reads_wrong" + TRUTH_SFX, "Tagged reads on the wrong\nparental haplotype (%)", reads_wrong_pct, False)


B_RECALL_SFX = "_%s_own_min0" % READ_TRUTH
_b_recall = {}


def read_recall_matpat_pct(db, key):
    """reads tagged to their haplotype / reads the matpat truth labels"""
    if db.sample_root not in _b_recall:
        _b_recall[db.sample_root] = matpat_recall(db, B_RECALL_SFX)
    return _b_recall[db.sample_root].get(key, float("nan"))


# B scores reads against the Q100 parental copies, not the truth VCF's tags
READ_RECALL_MATPAT = ("read_recall" + TRUTH_SFX, "Read haplotagging recall (%)",
                      read_recall_matpat_pct, False)
METRICS_A = METRICS[:3] + [READ_RECALL_MATPAT]
METRICS_B = METRICS[:3] + [READ_RECALL_MATPAT]


def fig_name(prefix, sample, ref):
    return "%s_%s" % (prefix, ref) if sample == "HG002" else "%s_%s_%s" % (prefix, sample, ref)


def phasers_on(platform):
    return [p for p in st.PHASERS if st.applies(platform, p)]


def legend_patches(items):
    """[(label, alpha, hatch)] in a neutral grey"""
    return [Patch(facecolor=st.MUTED, alpha=alpha, hatch=hatch, edgecolor="white" if hatch else "none",
                  label=label) for label, alpha, hatch in items]


def group_panels(platforms, width, height, legend=True):
    """A panel per platform, each as wide as its number of phaser groups, so
    a group is as wide on PacBio (five phasers) as on ONT (four); one shared
    y axis, the legend strip under the panels (legend_ax None without)."""
    widths = [len(phasers_on(p)) for p in platforms]
    fig = plt.figure(figsize=(width * st.MM, height * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(2 if legend else 1, len(platforms), width_ratios=widths,
                            height_ratios=[1, 0.1] if legend else None)
    axes = []
    for i in range(len(platforms)):
        ax = fig.add_subplot(grid[0, i], sharey=axes[0] if axes else None)
        if i:
            ax.tick_params(labelleft=False)
        axes.append(ax)
    legend_ax = None
    if legend:
        legend_ax = fig.add_subplot(grid[1, :])
        legend_ax.axis("off")
    return fig, legend_ax, axes


def group_panels_2x2(platforms, width, height, legend=True):
    """As group_panels, the platforms two to a row (ONT over PacBio), each row
    with its own widths, one y axis shared by all; for a page too narrow for
    four panels side by side."""
    rows = [platforms[i:i + 2] for i in range(0, len(platforms), 2)]
    fig = plt.figure(figsize=(width * st.MM, height * st.MM), constrained_layout=True)
    outer = fig.add_gridspec(len(rows) + (1 if legend else 0), 1,
                             height_ratios=[1] * len(rows) + ([0.12] if legend else []))
    axes = []
    for r, row in enumerate(rows):
        inner = outer[r, 0].subgridspec(1, len(row), width_ratios=[len(phasers_on(p)) for p in row])
        for i in range(len(row)):
            ax = fig.add_subplot(inner[0, i], sharey=axes[0] if axes else None)
            if i:
                ax.tick_params(labelleft=False)
            axes.append(ax)
    legend_ax = None
    if legend:
        legend_ax = fig.add_subplot(outer[len(rows), 0])
        legend_ax.axis("off")
    return fig, legend_ax, axes


# the panels of the paper: (st.PAPER_PANELS key, legend?, title) by file
# prefix; any other file is as wide as the page with its legend, no title.
# The legend goes under the last panel of each figure only (A4, B4).
PAPER = {
    "figA1": ("experiment_narrow", False, "Impact of additional variant types on switch error rate"),
    "figA3": ("experiment_narrow_legend", True, "Impact of additional variant types on phase-block NG50"),
    "figA4": ("experiment_legend", True, "Impact of additional variant types on read haplotagging recall"),
    "figB1": ("experiment_meth", False, "Effect of DNA methylation integration on switch error rate"),
    "figB3": ("experiment_meth", False, "Effect of DNA methylation integration on phase-block NG50"),
    "figB4": ("experiment_meth_legend", True, "Effect of DNA methylation integration on read haplotagging recall"),
    "figB5": ("meth_longphase", True, "Effects of methylation-assisted phasing on accuracy and continuity"),
}


def zoom_recall(ax, drawn):
    """A recall axis zoomed on the bars so the differences show: from under the
    lowest to over the highest, with room for the turned bar labels; no tick
    over 100%."""
    values = [v for v in drawn if finite(v)]
    lo, hi = min(values), max(values)
    span = max(hi - lo, 1.0)
    bottom = max(0.0, lo - 0.15 * span)
    top = hi + 0.35 * span
    ax.set_ylim(bottom, top)
    ticks = MaxNLocator(5).tick_values(bottom, min(top, 100))
    ax.set_yticks([t for t in ticks if bottom <= t <= 100])


def compare_groups(db, ctx, name, ref, metric, platforms, bars_of, legend, paper=None):
    """A panel per platform, a group per phaser; bars_of(platform, phaser) ->
    [(key, alpha, hatch)], the first the reference the others are compared with.
    paper: an entry of PAPER, the size, legend and title of a panel of the paper."""
    tag, ylabel, value, log = metric
    size, with_legend, title = paper or (None, True, None)
    width, height = st.PAPER_PANELS[size] if size else (st.FULL_WIDTH, st.PANEL_HEIGHT + 14)
    grid2 = bool(size) and size.startswith("experiment_grid")
    fig, legend_ax, axes = (group_panels_2x2 if grid2 else group_panels)(platforms, width, height, with_legend)
    # the labels are turned: neighbouring bars closer than ~2.4 mm take a smaller font
    n_bars = max(len(bars_of(platforms[0], phasers_on(platforms[0])[0])), 1)
    per_row = [platforms[i:i + 2] for i in range(0, len(platforms), 2)] if grid2 else [platforms]
    groups = max(sum(len(phasers_on(p)) for p in row) for row in per_row)
    # a narrow figure: the bars of a group wider apart, so their turned labels do not touch
    group_width = 0.72 if grid2 else 0.88 if size and size.startswith("experiment_narrow") else 0.8
    bar_mm = 0.85 * width / groups * group_width / n_bars
    label_size = 4 if bar_mm < 2.4 else None
    drawn = []
    for ax, platform in zip(axes, platforms):
        phasers = phasers_on(platform)
        for x, phaser in enumerate(phasers):
            group = bars_of(platform, phaser)
            xs, w = offsets(len(group), group_width)
            values = [value(db, k) for k, _, _ in group]
            values = [positive(v) for v in values] if log else values
            for dx, v, (_, alpha, hatch) in zip(xs, values, group):
                if bar(ax, x + dx, v, w, st.PHASER_COLOR[phaser], alpha, hatch) is None:
                    continue
                drawn.append(v)
                # the value alone, no change against the first bar of the group
                label_bar(ax, x + dx, v, w, value_text(v), label_size)
        group_ticks(ax, [st.PHASER_LABEL[p] for p in phasers])
        ax.set_title(st.PLATFORM_LABEL[platform])
    if not drawn:
        plt.close(fig)
        raise NoData("no %s values" % tag)
    scale_axes(axes, drawn, log, ylabel)
    if grid2:
        # the y label once per row, on its first panel
        for i, ax in enumerate(axes):
            ax.set_ylabel(ylabel if i % 2 == 0 else "")
            if i % 2 == 0:
                ax.tick_params(labelleft=True)
    if tag.startswith("read_recall"):
        zoom_recall(axes[0], drawn)
    top_legend(legend_ax, legend_patches(legend))
    paper_title(fig, title)
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
                          metric, platforms, bars_of, legend, PAPER.get("figA%d" % number))


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
    fig, legend_ax, axes = group_panels(platforms, st.FULL_WIDTH, st.PANEL_HEIGHT + 14)
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
        return [((ref, "truth", platform, DEPTH, phaser), PALE, None),
                ((ref, "truth", platform, DEPTH, "mp" + phaser), 1.0, None)]

    legend = [("Phasing method alone", PALE, None), ("+ MethPhaser", 1.0, None)]
    return compare_groups(db, ctx, fig_name("figB%d_%s" % (number, metric[0]), sample, ref), ref,
                          metric, platforms, bars_of, legend, PAPER.get("figB%d" % number))


# B5 is all LongPhase: a lighter blue of its hue than its own, on a small panel
LONGPHASE_LIGHT = "#4A7FB5"
LONGPHASE_METH = [("longphase", PALE, None, "LongPhase alone"),
                  ("mplongphase", 1.0, None, "LongPhase + MethPhaser"),
                  ("lpmod", 1.0, "////", "LongPhase + 5mC (modcall)")]


def figB5(db, ctx, sample, ref):
    """The two ways of adding 5mC to LongPhase side by side: a panel per metric
    of B1-B4, a group per platform, LongPhase alone (pale), with MethPhaser on
    top of its VCF (solid) and in its own 5mC mode (hatched)."""
    key = lambda p, ph: (ref, "truth", p, DEPTH, ph)
    platforms = [p for p in st.PLATFORMS
                 if any(db.has(key(p, ph)) for ph in ("mplongphase", "lpmod"))]
    if not platforms:
        raise NoData("no LongPhase methylation runs of %s on %s" % (sample, ref))
    size, _, title = PAPER["figB5"]
    w, h = st.PAPER_PANELS[size]
    fig = plt.figure(figsize=(w * st.MM, h * st.MM), constrained_layout=True)
    metrics = METRICS_B[:3]          # switch, Hamming, NG50: no read recall
    grid = fig.add_gridspec(2, len(metrics), height_ratios=[1, 0.1])
    xs, bw = offsets(len(LONGPHASE_METH))
    color = LONGPHASE_LIGHT
    drawn_any = False
    for i, (tag, ylabel, value, log) in enumerate(metrics):
        ax = fig.add_subplot(grid[0, i])
        drawn = []
        for x, platform in enumerate(platforms):
            for dx, (phaser, alpha, hatch, _) in zip(xs, LONGPHASE_METH):
                v = value(db, key(platform, phaser))
                v = positive(v) if log else v
                if bar(ax, x + dx, v, bw, color, alpha, hatch) is None:
                    continue
                drawn.append(v)
                label_bar(ax, x + dx, v, bw, value_text(v))
        group_ticks(ax, [st.PLATFORM_LABEL[p] for p in platforms])
        if not drawn:
            empty_panel(ax)
            continue
        drawn_any = True
        scale_axes([ax], drawn, log, ylabel)
        if tag.startswith("read_recall"):
            zoom_recall(ax, drawn)
    if not drawn_any:
        plt.close(fig)
        raise NoData("no LongPhase methylation values of %s on %s" % (sample, ref))
    legend_ax = fig.add_subplot(grid[1, :])
    legend_ax.axis("off")
    top_legend(legend_ax, [Patch(facecolor=color, alpha=alpha, hatch=hatch, edgecolor="white" if hatch else "none",
                                 label=label) for _, alpha, hatch, label in LONGPHASE_METH])
    paper_title(fig, title)
    return save(fig, ctx, fig_name("figB5_longphase_methylation", sample, ref))


def log_ticks_13(axis, lo, hi):
    """1 and 3 of each decade between lo and hi, plain numbers: few enough for a narrow panel"""
    ticks = [m * 10.0 ** k for k in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi)) + 1)
             for m in (1, 3) if lo <= m * 10.0 ** k <= hi]
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    axis.set_minor_locator(NullLocator())


def figB6_scatter(db, ctx, sample, ref):
    """Switch error rate against phase-block NG50 with and without DNA
    methylation, ONT 5 kHz over PacBio Revio, each with its own axes, a tall
    narrow column beside the region view and B5. A phaser has its own marker and colour: open alone, an arrow to the
    filled + MethPhaser point; LongPhase in its own 5mC mode is a star at the
    end of a dotted arrow. The switch axis runs from more to fewer errors, so
    the top right corner is best."""
    platforms = [p for p in ("ont5k", "revio")
                 if any(db.has((ref, "truth", p, DEPTH, "mp" + ph)) for ph in phasers_on(p))]
    if not platforms:
        raise NoData("no MethPhaser runs of %s on %s" % (sample, ref))
    w, h = st.PAPER_PANELS["meth_scatter"]
    fig = plt.figure(figsize=(w * st.MM, h * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(len(platforms) + 1, 1, height_ratios=[1] * len(platforms) + [0.3])
    shown = []
    for i, platform in enumerate(platforms):
        ax = fig.add_subplot(grid[i, 0])
        xs, ys = [], []
        for ph in phasers_on(platform):
            k0 = (ref, "truth", platform, DEPTH, ph)
            x0, y0 = switch_pct(db, k0), n50_mb(db, k0)
            if not (finite(x0) and finite(y0) and x0 > 0 and y0 > 0):
                continue
            c, m = st.PHASER_COLOR[ph], st.PHASER_MARKER[ph]
            if ph not in shown:
                shown.append(ph)
            ax.plot(x0, y0, marker=m, markersize=3.6, markerfacecolor="white", markeredgecolor=c,
                    markeredgewidth=0.9, linestyle="none", zorder=4)
            xs.append(x0)
            ys.append(y0)
            for other, marker, style in (("mp" + ph, m, "-"), ("lpmod" if ph == "longphase" else None, "*", ":")):
                if not other:
                    continue
                k1 = (ref, "truth", platform, DEPTH, other)
                x1, y1 = switch_pct(db, k1), n50_mb(db, k1)
                if not (finite(x1) and finite(y1) and x1 > 0 and y1 > 0):
                    continue
                ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                            arrowprops=dict(arrowstyle="-|>", color=c, linewidth=0.7, linestyle=style,
                                            shrinkA=2.5, shrinkB=2.5, mutation_scale=5), zorder=3)
                ax.plot(x1, y1, marker=marker, markersize=3.6 if marker != "*" else 5.5, color=c,
                        markeredgecolor="white", markeredgewidth=0.3, linestyle="none", zorder=5)
                xs.append(x1)
                ys.append(y1)
        if not xs:
            empty_panel(ax)
            continue
        ax.set_title(st.PLATFORM_LABEL[platform])
        ax.set_xscale("log")
        ax.set_yscale("log")
        xlo, xhi, ylo, yhi = min(xs) / 1.5, max(xs) * 1.5, min(ys) / 1.2, max(ys) * 1.2
        ax.set_xlim(xhi, xlo)          # reversed: fewer switches to the right, the top right corner is best
        ax.set_ylim(ylo, yhi)
        log_ticks_13(ax.xaxis, xlo, xhi)
        log_ticks_13(ax.yaxis, ylo, yhi)
        if len(ax.get_yticks()) < 2:
            ax.yaxis.set_major_locator(FixedLocator([float("%.2g" % min(ys)), float("%.2g" % max(ys))]))
        ax.grid(color=st.GRID, linewidth=0.5)
        ax.set_axisbelow(True)
        ax.set_xlabel("Switch error rate (%)")
        ax.set_ylabel("Phase-block NG50 (Mb)")
    if not shown:
        plt.close(fig)
        raise NoData("no switch rates of the MethPhaser runs of %s on %s" % (sample, ref))
    legend_ax = fig.add_subplot(grid[len(platforms), 0])
    legend_ax.axis("off")
    handles = [Line2D([], [], marker=st.PHASER_MARKER[ph], linestyle="none", markersize=3,
                      color=st.PHASER_COLOR[ph], markeredgecolor="white", markeredgewidth=0.3,
                      label=st.PHASER_LABEL[ph]) for ph in shown]
    handles += [Line2D([], [], marker="o", linestyle="none", markersize=3, markerfacecolor="white",
                       markeredgecolor=st.MUTED, label="Alone"),
                Line2D([], [], marker="o", linestyle="none", markersize=3, color=st.MUTED, label="+ MethPhaser"),
                Line2D([], [], marker="*", linestyle="none", markersize=5, color=st.MUTED, label="LongPhase + 5mC")]
    legend_ax.legend(handles=handles, loc="center", ncol=2, borderaxespad=0, columnspacing=0.8,
                     handletextpad=0.2, handlelength=1.0, labelspacing=0.3, frameon=False)
    paper_title(fig, "Switch error rate and NG50 with and without methylation")
    return save(fig, ctx, fig_name("figB6_meth_scatter", sample, ref))


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


# ------------------------------------------- horizontal bars, beside a figure

def side_bars(ctx, name, size, title, platforms, columns, legend=None, xmax=None):
    """Horizontal bars on log x axes for a narrow column beside another figure:
    a row of panels per platform (each as tall as its rows), a column per
    metric with its own x axis, the row labels on the first column only.
    columns: [(x label, rows_of)], rows_of(platform) -> [(row label,
    [(value, colour, alpha, hatch)])], the same rows in every column.
    legend: [(label, alpha, hatch)] in a strip under the panels; xmax: the
    right end of the x axes (default 8 times the largest value; at least 1.6
    times it)."""
    w, h = st.PAPER_PANELS[size]
    rows = [{p: rows_of(p) for p in platforms} for _, rows_of in columns]
    platforms = [p for p in platforms
                 if any(finite(v) and v > 0 for r in rows for _, bars in r[p] for v, _, _, _ in bars)]
    if not platforms:
        raise NoData("no values for %s" % name)
    fig = plt.figure(figsize=(w * st.MM, h * st.MM), constrained_layout=True)
    ratios = [len(rows[0][p]) for p in platforms]
    # the legend in a row when each entry gets 22 mm, else one under the other
    ncol = (len(legend) if w / len(legend) >= 22 else 1) if legend else 1
    grid = fig.add_gridspec(len(platforms) + bool(legend), len(columns),
                            height_ratios=ratios + [0.07 * sum(ratios) * math.ceil(len(legend) / ncol)]
                            if legend else ratios)
    # the value labels: 5 pt where a bar is 2.2 mm high or more and the column 20 mm wide, else 4 pt
    n_rows = sum(len(rows[0][p]) for p in platforms)
    n_bars = max(len(bars) for p in platforms for _, bars in rows[0][p])
    bar_mm = 0.7 * h / n_rows * 0.8 / n_bars
    label_size = st.SMALL_FONT if bar_mm >= 2.2 and w / len(columns) > 20 else 4
    for c, (xlabel, _) in enumerate(columns):
        axes, drawn = [], []
        for i, platform in enumerate(platforms):
            ax = fig.add_subplot(grid[i, c], sharex=axes[0] if axes else None)
            axes.append(ax)
            for y, (_, bars) in enumerate(rows[c][platform]):
                n = len(bars)
                for j, (v, color, alpha, hatch) in enumerate(bars):
                    if not (finite(v) and v > 0):
                        continue
                    yy = y - 0.4 + 0.8 * (j + 0.5) / n
                    ax.barh(yy, v, height=0.8 / n * 0.9, color=color, alpha=alpha, hatch=hatch,
                            edgecolor="white" if hatch else "none", linewidth=0, zorder=3)
                    ax.annotate(value_text(v), (v, yy), xytext=(1.5, 0), textcoords="offset points",
                                ha="left", va="center", fontsize=label_size,
                                color=st.INK, zorder=4, annotation_clip=False)
                    drawn.append(v)
            ax.set_yticks(range(len(rows[c][platform])))
            ax.set_yticklabels([label for label, _ in rows[c][platform]] if c == 0 else [])
            ax.set_ylim(len(rows[c][platform]) - 0.5, -0.5)
            ax.set_title(st.PLATFORM_LABEL[platform] if c == 0 or len(columns) == 1 else "")
            ax.grid(axis="x", color=st.GRID, linewidth=0.5, zorder=0)
            ax.set_axisbelow(True)
            if i < len(platforms) - 1:
                ax.tick_params(labelbottom=False)
        if not drawn:
            continue
        lo, hi = min(drawn), max(drawn)
        # xmax: room for the value labels kept, should a value come near it
        right = max(xmax, hi * 1.6) if xmax else hi * 8
        axes[0].set_xscale("log")
        axes[0].set_xlim(lo / 2, right)
        log_ticks(axes[0].xaxis, lo / 2, right)
        if xmax:
            # a short axis: the ticks that fit, so its right end is labelled
            ticks = [t for t in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10) if lo / 2 <= t <= right]
            axes[0].xaxis.set_major_locator(FixedLocator(ticks))
            axes[0].xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
        elif w / len(columns) <= 30:
            # a narrow column: decades only, every other one when there are many
            first, last = math.ceil(math.log10(lo / 2)), math.floor(math.log10(hi * 8))
            ticks = [10.0 ** k for k in range(first, last + 1)]
            ticks = ticks[::2] if len(ticks) > 3 else ticks
            axes[0].xaxis.set_major_locator(FixedLocator(ticks))
            axes[0].xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
        axes[-1].set_xlabel(xlabel)
    if legend:
        legend_ax = fig.add_subplot(grid[len(platforms), :])
        legend_ax.axis("off")
        legend_ax.legend(handles=legend_patches(legend), loc="center", ncol=ncol, borderaxespad=0,
                         handlelength=1.2, handletextpad=0.4, labelspacing=0.3)
    paper_title(fig, title)
    return save(fig, ctx, name)


def figB3_side(db, ctx, sample, ref):
    """B3 (phase-block NG50, alone and with MethPhaser) as horizontal bars, a
    narrow column beside A1 and A3."""
    platforms = [p for p in st.PLATFORMS
                 if any(db.has((ref, "truth", p, DEPTH, m)) for m in st.METH_PHASERS if st.applies(p, m))]

    def rows_of(platform):
        return [(st.PHASER_LABEL[ph],
                 [(n50_mb(db, (ref, "truth", platform, DEPTH, ph)), st.PHASER_COLOR[ph], PALE, None),
                  (n50_mb(db, (ref, "truth", platform, DEPTH, "mp" + ph)), st.PHASER_COLOR[ph], 1.0, None)])
                for ph in phasers_on(platform) if db.has((ref, "truth", platform, DEPTH, "mp" + ph))]
    return side_bars(ctx, fig_name("figB3_side_ng50", sample, ref), "meth_side",
                     "Effect of DNA methylation integration on phase-block NG50",
                     platforms, [("Phase-block NG50 (Mb)", rows_of)],
                     legend=[("Phasing method alone", PALE, None), ("+ MethPhaser", 1.0, None)], xmax=2)


def figC1_side(db, ctx, sample, ref, depth):
    """Switch error rate and phase-block NG50 side by side as horizontal bars:
    each read-based phaser on the truth SNPs (light) and on the SNPs of the
    hifiasm assembly without parents (dark), told apart in the legend - a
    narrow column beside the chromosome figure of the switch errors."""
    platforms = [p for p in st.PLATFORMS
                 if any(db.has((ref, "hifiasm", p, depth, ph)) for ph in phasers_on(p))]

    def rows(value):
        def rows_of(platform):
            return [(st.PHASER_LABEL[ph],
                     [(value(db, (ref, "truth", platform, depth, ph)), st.PHASER_COLOR[ph], PALE, None),
                      (value(db, (ref, "hifiasm", platform, depth, ph)), st.PHASER_COLOR[ph], 1.0, None)])
                    for ph in phasers_on(platform) if db.has((ref, "hifiasm", platform, depth, ph))]
        return rows_of
    return side_bars(ctx, fig_name("figC1_side_switch_ng50_%dx" % depth, sample, ref), "asm_side",
                     "Phasers on the truth SNPs and on the hifiasm SNPs, %d×" % depth,
                     platforms, [("Switch error\nrate (%)", rows(switch_pct)),
                                 ("Phase-block\nNG50 (Mb)", rows(n50_mb))],
                     legend=[("Truth SNPs", PALE, None), ("hifiasm SNPs", 1.0, None)])


def hifiasm_common(ctx, sample, ref, depth):
    """(platform, phaser) -> {"truth": (switch %, NG50 Mb), "hifiasm": (...)}:
    each phaser on the truth SNPs and on the hifiasm SNPs, both scored on the
    SNPs phased in both runs (identical position and alleles), from
    <root>/asm_common_hifiasm/summary_<depth>x.txt"""
    path = os.path.join(sample_root(ctx.root, sample), "asm_common_hifiasm", "summary_%dx.txt" % depth)
    out = {}
    if os.path.exists(path):
        with open(path) as fh:
            h = fh.readline().split()
            for line in fh:
                r = dict(zip(h, line.split()))
                if r["reference"] == ref:
                    out[(r["platform"], r["phaser"])] = {
                        "truth": (float(r["switch_truth"]), float(r["ng50_truth"])),
                        "hifiasm": (float(r["switch_hifiasm"]), float(r["ng50_hifiasm"]))}
    return out


def figC6_scatter(db, ctx, sample, ref, depth):
    """Switch error rate against phase-block NG50, the read-based phasers on
    the truth SNPs and on the SNPs of the hifiasm assembly without parents,
    both on the SNPs phased in both runs (hifiasm_common), as B6: a panel per
    platform, each with its own axes, a phaser in its own marker and colour,
    open on the truth SNPs, an arrow to the filled point on the hifiasm SNPs.
    The switch axis runs from more to fewer errors, so the top right corner is
    best. The same slot as figC1_side."""
    common = hifiasm_common(ctx, sample, ref, depth)
    platforms = [p for p in st.PLATFORMS if any((p, ph) in common for ph in phasers_on(p))]
    if not platforms:
        raise NoData("no asm_common_hifiasm/summary_%dx.txt rows of %s on %s - run the common-SNP step"
                     % (depth, sample, ref))
    w, h = st.PAPER_PANELS["asm_side"]
    fig = plt.figure(figsize=(w * st.MM, h * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(len(platforms) + 1, 1, height_ratios=[1] * len(platforms) + [0.3])
    shown = []
    for i, platform in enumerate(platforms):
        ax = fig.add_subplot(grid[i, 0])
        xs, ys = [], []
        for ph in phasers_on(platform):
            pts = []
            for inp in ("truth", "hifiasm"):
                x, y = common.get((platform, ph), {}).get(inp, (math.nan, math.nan))
                pts.append((x, y) if finite(x) and finite(y) and x > 0 and y > 0 else None)
            if not any(pts):
                continue
            c, m = st.PHASER_COLOR[ph], st.PHASER_MARKER[ph]
            if ph not in shown:
                shown.append(ph)
            if all(pts):
                ax.annotate("", xy=pts[1], xytext=pts[0],
                            arrowprops=dict(arrowstyle="-|>", color=c, linewidth=0.7,
                                            shrinkA=2.5, shrinkB=2.5, mutation_scale=5), zorder=3)
            for pt, filled in zip(pts, (False, True)):
                if not pt:
                    continue
                ax.plot(*pt, marker=m, markersize=3.6, linestyle="none", zorder=4 + filled,
                        markerfacecolor=c if filled else "white", markeredgecolor="white" if filled else c,
                        markeredgewidth=0.3 if filled else 0.9)
                xs.append(pt[0])
                ys.append(pt[1])
        if not xs:
            empty_panel(ax)
            continue
        ax.set_title(st.PLATFORM_LABEL[platform])
        ax.set_xscale("log")
        ax.set_yscale("log")
        xlo, xhi, ylo, yhi = min(xs) / 1.5, max(xs) * 1.5, min(ys) / 1.2, max(ys) * 1.2
        ax.set_xlim(xhi, xlo)          # reversed: fewer switches to the right, the top right corner is best
        ax.set_ylim(ylo, yhi)
        log_ticks_13(ax.xaxis, xlo, xhi)
        log_ticks_13(ax.yaxis, ylo, yhi)
        if len(ax.get_xticks()) < 2:
            ax.xaxis.set_major_locator(FixedLocator([float("%.2g" % min(xs)), float("%.2g" % max(xs))]))
        if len(ax.get_yticks()) < 2:
            ax.yaxis.set_major_locator(FixedLocator([float("%.2g" % min(ys)), float("%.2g" % max(ys))]))
        ax.grid(color=st.GRID, linewidth=0.5)
        ax.set_axisbelow(True)
        ax.set_ylabel("NG50 (Mb)")
        if i == len(platforms) - 1:
            ax.set_xlabel("Switch error rate (%)")
    if not shown:
        plt.close(fig)
        raise NoData("no switch rates of the hifiasm SNP runs of %s on %s" % (sample, ref))
    legend_ax = fig.add_subplot(grid[len(platforms), 0])
    legend_ax.axis("off")
    handles = [Line2D([], [], marker=st.PHASER_MARKER[ph], linestyle="none", markersize=3,
                      color=st.PHASER_COLOR[ph], markeredgecolor="white", markeredgewidth=0.3,
                      label=st.PHASER_LABEL[ph]) for ph in shown]
    handles += [Line2D([], [], marker="o", linestyle="none", markersize=3, markerfacecolor="white",
                       markeredgecolor=st.MUTED, label="Truth SNPs"),
                Line2D([], [], marker="o", linestyle="none", markersize=3, color=st.MUTED, label="hifiasm SNPs")]
    legend_ax.legend(handles=handles, loc="center", ncol=2, borderaxespad=0, columnspacing=0.8,
                     handletextpad=0.2, handlelength=1.0, labelspacing=0.3, frameon=False)
    paper_title(fig, "Switch error rate and NG50 with truth and assembly-derived SNPs")
    return save(fig, ctx, fig_name("figC6_switch_ng50_%dx" % depth, sample, ref))


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
                    for number, metric in enumerate(METRICS_A, 1):
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
                    results.append(run("B %s %s LongPhase 5mC" % (sample, ref),
                                       figB5, db, ctx, sample, ref))
                    results.append(run("B %s %s NG50 beside A" % (sample, ref), figB3_side, db, ctx, sample, ref))
                    results.append(run("B %s %s switch vs NG50" % (sample, ref), figB6_scatter, db, ctx, sample, ref))
        if "C" in only:
            db = Metrics(ctx.root, sample, st.REGION, readacc_rules=ctx.rules).load(
                refs, ["truth"] + ASSEMBLERS, st.PLATFORMS, [ctx.asm_depth], st.PHASERS + ["dipcall"])
            c_metrics = [(METRICS[0], True), (METRICS[1], False), (METRICS[2], False),
                         (READS_WRONG, False), (METRICS[3], False)]
            for ref in refs:
                for number, (metric, with_hifiasm) in enumerate(c_metrics, 1):
                    results.append(run("C %s %s %s" % (sample, ref, metric[0]), figC_metric, db, ctx,
                                       sample, ref, number, metric, ctx.asm_depth, with_hifiasm))
                results.append(run("C %s %s switch beside the chromosomes" % (sample, ref), figC1_side, db, ctx,
                                   sample, ref, ctx.asm_depth))
                results.append(run("C %s %s switch vs NG50" % (sample, ref), figC6_scatter, db, ctx,
                                   sample, ref, ctx.asm_depth))
    print("%d written, %d skipped for lack of data, %d failed -> %s" % (
        results.count("ok"), results.count("skip"), results.count("failed"), ctx.out))
    return 1 if "failed" in results else 0


if __name__ == "__main__":
    sys.exit(main())
