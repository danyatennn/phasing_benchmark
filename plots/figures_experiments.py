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
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator

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
    "figB5": ("experiment_legend", True, "MethPhaser and LongPhase 5mC mode (modcall) on LongPhase"),
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
    fig, legend_ax, axes = group_panels(platforms, width, height, with_legend)
    # the labels are turned: neighbouring bars closer than ~2.4 mm take a smaller font
    n_bars = max(len(bars_of(platforms[0], phasers_on(platforms[0])[0])), 1)
    bar_mm = 0.85 * width / sum(len(phasers_on(p)) for p in platforms) * 0.8 / n_bars
    label_size = 4 if bar_mm < 2.4 else None
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
                label_bar(ax, x + dx, v, w, value_text(v), label_size)
        group_ticks(ax, [st.PHASER_LABEL[p] for p in phasers])
        ax.set_title(st.PLATFORM_LABEL[platform])
    if not drawn:
        plt.close(fig)
        raise NoData("no %s values" % tag)
    scale_axes(axes, drawn, log, ylabel)
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
    color = st.PHASER_COLOR["longphase"]
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
    top_legend(legend_ax, legend_patches([(label, alpha, hatch)
                                          for _, alpha, hatch, label in LONGPHASE_METH]))
    paper_title(fig, title)
    return save(fig, ctx, fig_name("figB5_longphase_methylation", sample, ref))


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

def side_bars(ctx, name, size, title, platforms, columns):
    """Horizontal bars on log x axes for a narrow column beside another figure:
    a row of panels per platform (each as tall as its rows), a column per
    metric with its own x axis, the row labels on the first column only.
    columns: [(x label, rows_of)], rows_of(platform) -> [(row label,
    [(value, colour, alpha, hatch)])], the same rows in every column."""
    w, h = st.PAPER_PANELS[size]
    rows = [{p: rows_of(p) for p in platforms} for _, rows_of in columns]
    platforms = [p for p in platforms
                 if any(finite(v) and v > 0 for r in rows for _, bars in r[p] for v, _, _, _ in bars)]
    if not platforms:
        raise NoData("no values for %s" % name)
    fig = plt.figure(figsize=(w * st.MM, h * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(len(platforms), len(columns), height_ratios=[len(rows[0][p]) for p in platforms])
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
        axes[0].set_xscale("log")
        axes[0].set_xlim(lo / 2, hi * 8)
        log_ticks(axes[0].xaxis, lo / 2, hi * 8)
        if w / len(columns) <= 30:
            # a narrow column: decades only, every other one when there are many
            first, last = math.ceil(math.log10(lo / 2)), math.floor(math.log10(hi * 8))
            ticks = [10.0 ** k for k in range(first, last + 1)]
            ticks = ticks[::2] if len(ticks) > 3 else ticks
            axes[0].xaxis.set_major_locator(FixedLocator(ticks))
            axes[0].xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
        axes[-1].set_xlabel(xlabel)
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
                     platforms, [("Phase-block NG50 (Mb)", rows_of)])


def figC1_side(db, ctx, sample, ref, depth):
    """Switch error rate and phase-block NG50 side by side as horizontal bars:
    the read-based phasers on the truth input, then the trio assemblies - a
    narrow column beside the chromosome figure of the switch errors."""
    platforms = [p for p in st.PLATFORMS
                 if any(db.has((ref, a, p, depth, "dipcall")) for a in ASSEMBLERS[1:])]

    def rows(value):
        def rows_of(platform):
            out = [(st.PHASER_LABEL[ph], [(value(db, (ref, "truth", platform, depth, ph)), st.PHASER_COLOR[ph], 1.0, None)])
                   for ph in phasers_on(platform)]
            out += [(ASM_LABEL[a], [(value(db, (ref, a, platform, depth, "dipcall")), ASM_COLOR[a], 1.0, None)])
                    for a in ASSEMBLERS[1:] if db.has((ref, a, platform, depth, "dipcall"))]
            return out
        return rows_of
    return side_bars(ctx, fig_name("figC1_side_switch_ng50_%dx" % depth, sample, ref), "asm_side",
                     "Alignment-based and assembly-based phasing, %d×" % depth,
                     platforms, [("Switch error\nrate (%)", rows(switch_pct)),
                                 ("Phase-block\nNG50 (Mb)", rows(n50_mb))])


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
    print("%d written, %d skipped for lack of data, %d failed -> %s" % (
        results.count("ok"), results.count("skip"), results.count("failed"), ctx.out))
    return 1 if "failed" in results else 0


if __name__ == "__main__":
    sys.exit(main())
