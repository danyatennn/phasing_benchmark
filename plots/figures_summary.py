#!/usr/bin/env python3
"""Summary figures of the paper, one PDF per subfigure, from the merged tables:

  figR5_untagged_reasons_<phaser>      why the remaining reads are not tagged:
      a donut per platform and reference, the untagged reads split by reason
      (metrics_reads/summary_classes.tsv, 30x)
  figR6_<metric>_<ref>                 phasing rate (truth het SNPs phased in a
      block of two or more) and haplotag rate (reads tagged) as bars: a row per
      phaser, coverage on x, a bar per platform; the y axis is shared by both
      references (metrics/summary_auto.tsv, summary_tagbases)
  figR7_difficult_drop_<phaser>        from genome-wide (open) to all difficult
      regions (filled): switch error rate against truth phased, an arrow per
      platform, a panel per reference (summary_auto, summary_phase_strata, 30x)
  figR8_callers_strata_<ref>           F1 of each caller's het SNPs per genomic
      region, a row per platform (metrics_strata_auto/summary_calls_strata.tsv)
  figR9_rescued_<platform>             reads untagged on GRCh38, and those of
      them that CHM13 tags, WhatsHap and LongPhase, as an UpSet plot
      (metrics_reads/rescue/, bin/rescue_reads.sh)

Input: the GIAB truth (--caller truth), so the phaser is the only variable.

  python plots/figures_summary.py --png
  python plots/figures_summary.py --only R5,R9 --png
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Patch
import matplotlib.colors as matplotlib_colors

import style as st
from data import sample_root
from figures import NoData, colorbar, heatmap, log_ticks, save, series_handle, stratum_label, top_legend, y_grid


# ------------------------------------------------------------------ helpers

def table(path):
    """rows of a TSV with a header, [] when the file is not there"""
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return math.nan
    return x if math.isfinite(x) else math.nan


def thousands(n):
    return "{:,}".format(int(round(n)))


def short(n):
    if n >= 1e6:
        return "%.2fM" % (n / 1e6)
    if n >= 1e4:
        return "%.0fk" % (n / 1e3)
    if n >= 1e3:
        return "%.1fk" % (n / 1e3)
    return "%d" % n


def axis_break(ax):
    """Two short slanted strokes at the bottom of the y axis: it does not start at 0"""
    kw = dict(transform=ax.transAxes, color=st.INK, clip_on=False, linewidth=0.6)
    for dy in (0.0, 0.025):
        ax.plot([-0.012, 0.012], [0.02 + dy, 0.045 + dy], **kw)


def path_of(ctx, *parts):
    return os.path.join(sample_root(ctx.root, ctx.sample), *parts)


# ------------------------------------------------- R5: why reads stay untagged

# in drawing order, clockwise from 12 o'clock: the small classes first, so
# their labels share the roomy right side instead of crowding the top
CLASSES = [
    ("lowmapq_paralog", "Low MAPQ, paralogous region", "#9C3D2E", None),
    ("lowmapq_other", "Low MAPQ, elsewhere", "#DA9A82", None),
    ("snps_unphased", "Het SNPs, not phased", "#4F6D80", None),
    ("phased_unassigned", "Phased, no haplotype", "#B5C4CF", None),
    ("no_snp_outside_bed", "No het SNP, outside benchmark", "#EADCA6", "////"),
    ("no_snp", "No het SNP", "#C7A23F", None),
]


def spread(labels, gap, low=-1.35, high=1.35):
    """labels [(y, ...)] of one side, top to bottom, at least gap apart and
    inside [low, high]"""
    labels = sorted(labels, key=lambda t: -t[0])
    ys = []
    for t in labels:
        y = min(t[0], high) if not ys else min(t[0], ys[-1] - gap)
        ys.append(y)
    if ys and ys[-1] < low:
        ys = [y + (low - ys[-1]) for y in ys]
    return [(y,) + t[1:] for y, t in zip(ys, labels)]


def donut(ax, counts, center, min_label, min_slice, classes=None):
    total = sum(counts.values())
    vals = [(k, col, hatch, counts.get(k, 0)) for k, _, col, hatch in (classes or CLASSES)
            if counts.get(k, 0) > 0]
    # a slice of a fraction of a percent is a hairline: drawn at min_slice of
    # the ring at least (the labels keep the true numbers; said in the figure)
    drawn = [max(v[3], min_slice * total) for v in vals]
    wedges, _ = ax.pie(drawn, radius=1.0, startangle=90, counterclock=False,
                       colors=[v[1] for v in vals],
                       wedgeprops=dict(width=0.3, edgecolor="white", linewidth=0.8))
    for w, v in zip(wedges, vals):
        if v[2]:
            w.set_hatch(v[2])
    ax.text(0, 0, center, ha="center", va="center", color=st.INK, fontsize=st.FONT_SIZE - 1,
            linespacing=1.15)
    sides = {1: [], -1: []}
    for w, v in zip(wedges, vals):
        share = v[3] / total
        if share < min_label:
            continue
        ang = math.radians((w.theta1 + w.theta2) / 2)
        x, y = math.cos(ang), math.sin(ang)
        sides[1 if x >= 0 else -1].append((1.15 * y, x, y, "%s (%.1f%%)" % (short(v[3]), 100 * share)))
    for sign, labels in sides.items():
        for ty, x, y, text in spread(labels, 0.34, low=-1.45, high=1.45):
            # the line stops a few points short of the text
            ax.annotate(text, xy=(1.01 * x, 1.01 * y), xytext=(sign * 1.32, ty),
                        ha="left" if sign > 0 else "right", va="center", color=st.INK,
                        fontsize=st.FONT_SIZE - 1,
                        arrowprops=dict(arrowstyle="-", color=st.MUTED, linewidth=0.5, shrinkA=4, shrinkB=0))
    ax.set_xlim(-2.45, 2.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")


def donut_figure(ctx, group_type, classes, center_word, title, name, share=True):
    """Donuts of the reads of genome-wide read classes (bin/read_classes.sh),
    a row per platform, the references side by side. group_type: class (why
    untagged) or confidence (the SNP votes behind the read's haplotype)."""
    rows = [r for r in table(path_of(ctx, "metrics_reads", "summary_classes.tsv"))
            if r.get("caller") == ctx.caller and r.get("depth") == str(ctx.depth)
            and r.get("phaser") == ctx.phaser and r.get("stratum") == "genome"]
    if not rows:
        raise NoData("no read classes for %s %s %dx" % (ctx.caller, ctx.phaser, ctx.depth))
    data = {}
    for r in rows:
        key = (r["reference"], r["platform"])
        if r["group_type"] in (group_type, "all"):
            data.setdefault(key, {})[r["group"]] = num(r["reads"])
    refs = [r for r in st.REFERENCES if any(k[0] == r for k in data)]
    platforms = [p for p in st.PLATFORMS if any(k[1] == p for k in data)]
    # a row per platform, the two references side by side
    top, bottom = 0.9, 0.1
    fig, axes = plt.subplots(len(platforms), len(refs), squeeze=False,
                             figsize=(st.FULL_WIDTH * st.MM, (26 + 40 * len(platforms)) * st.MM))
    fig.subplots_adjust(left=0.1, right=0.99, top=top, bottom=bottom, wspace=0.0, hspace=0.05)
    for i, platform in enumerate(platforms):
        for j, ref in enumerate(refs):
            ax = axes[i][j]
            counts = data.get((ref, platform))
            if not counts:
                ax.axis("off")
                continue
            part = {k: counts.get(k, 0) for k, *_ in classes}
            total_reads = counts.get("all", math.nan)
            n = sum(part.values())
            if share and total_reads == total_reads:
                center = "%s\n%s\n(%.1f%%)" % (short(n), center_word, 100 * n / total_reads)
            else:
                center = "%s\n%s" % (short(n), center_word)
            donut(ax, part, center, ctx.min_label, ctx.min_slice, classes)
            if i == 0:
                ax.set_title(st.REFERENCE_LABEL[ref], pad=4)
        fig.text(0.02, top - (top - bottom) * (i + 0.5) / len(platforms), st.PLATFORM_LABEL[platform],
                 rotation=90, ha="center", va="center")
    handles = [Patch(facecolor=c, hatch=h, edgecolor="white", label=l) for _, l, c, h in classes]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.55, 0.02))
    if ctx.min_slice > 0:
        fig.text(0.55, 0.005, "Slices under %g%% are drawn at %g%% so they can be seen; the labels give the "
                 "true numbers." % (100 * ctx.min_slice, 100 * ctx.min_slice), ha="center", va="bottom",
                 color=st.MUTED, fontsize=st.FONT_SIZE - 1)
    fig.suptitle("%s: %s, %s input, %dx, chr1-22" % (
        title, st.PHASER_LABEL[ctx.phaser], st.CALLER_LABEL.get(ctx.caller, ctx.caller), ctx.depth), y=0.99)
    return save(fig, ctx, "%s_%s" % (name, ctx.phaser))


def fig_r5(ctx):
    return donut_figure(ctx, "class", CLASSES, "untagged", "Why reads stay untagged", "figR5_untagged_reasons")


# --------------------------------------------- R11: SNP votes behind a read

# the read's votes: at each phased het SNP of its block it carries the allele of
# haplotype 1 or 2 (plots/read_classes.py); a read without a tag is counted
# in the block that would count most
VOTES = [
    ("high", "≥2 votes, all for one haplotype", "#2F5D50", None),
    ("medium", "≥2 votes, majority with dissent", "#7FA693", None),
    ("single", "a single vote", "#D7B65D", None),
    ("tie", "as many votes for each", "#D9876A", None),
    ("none", "no vote", "#C9D1D8", None),
]


def fig_r11(ctx):
    return donut_figure(ctx, "confidence", VOTES, "reads", "SNP votes behind each read's haplotype",
                        "figR11_read_votes", share=False)


# ------------------------------------------- R6: phasing and haplotag rates

def rate_values(ctx):
    """metric -> {(ref, phaser, platform, depth): percent}"""
    auto = [r for r in table(path_of(ctx, "metrics", "summary_auto.tsv")) if r.get("caller") == ctx.caller]
    out = {"phased": {}, "tagged": {}}
    for r in auto:
        k = (r["reference"], r["phaser"], r["platform"], int(float(r["depth"])))
        tt = num(r["truth_total"])
        if tt:
            out["phased"][k] = 100 * num(r["covered"]) / tt
        out["tagged"][k] = 100 * num(r["tagging_completeness_min"])
    if ctx.tagged_unit == "bases":
        out["tagged"] = {}
        for r in table(path_of(ctx, "metrics", "summary_tagbases_auto_min2.tsv")):
            if r.get("caller") == ctx.caller:
                k = (r["reference"], r["phaser"], r["platform"], int(float(r["depth"])))
                out["tagged"][k] = 100 * num(r["bases_tagged_fraction"])
    return out


RATE_LABEL = {"phased": "Truth het SNPs phased (%)", "tagged": "Reads tagged (%)"}


def fig_r6(ctx):
    values = rate_values(ctx)
    if ctx.tagged_unit == "bases":
        RATE_LABEL["tagged"] = "Read bases tagged (%)"
    made = []
    for metric in ("phased", "tagged"):
        v = values[metric]
        if not v:
            continue
        finite = [x for x in v.values() if x == x]
        floor = max(0, 10 * math.floor((min(finite) - 3) / 10))
        for ref in st.REFERENCES:
            keys = [k for k in v if k[0] == ref]
            if not keys:
                continue
            phasers = [p for p in ctx.phasers if any(k[1] == p for k in keys)]
            depths = sorted({k[3] for k in keys})
            platforms = [p for p in st.PLATFORMS if any(k[2] == p for k in keys)]
            fig, axes = plt.subplots(len(phasers), 1, squeeze=False, sharex=True,
                                     figsize=(st.FULL_WIDTH * st.MM, (16 + 34 * len(phasers)) * st.MM),
                                     constrained_layout=True)
            nb = len(platforms)
            width = 0.8 / nb
            for ax, phaser in zip(axes[:, 0], phasers):
                for i, d in enumerate(depths):
                    for j, p in enumerate(platforms):
                        y = v.get((ref, phaser, p, d), math.nan)
                        if y != y:
                            continue
                        x = i + (j - (nb - 1) / 2) * width
                        ax.bar(x, y - floor, bottom=floor, width=width * 0.92,
                               color=st.PLATFORM_COLOR[p], linewidth=0, zorder=2)
                        dark = p in ("ont5k", "revio")
                        inside = y - floor > 0.35 * (100 - floor)
                        ax.text(x, (y - 1.2 * (100 - floor) / 100) if inside else y + 0.8 * (100 - floor) / 100,
                                "%.1f" % y, rotation=90, ha="center", va="top" if inside else "bottom",
                                fontsize=st.FONT_SIZE - 1,
                                color="white" if (inside and dark) else st.INK, zorder=3)
                ax.set_ylim(floor, 100)
                ax.set_title(st.PHASER_LABEL[phaser], loc="left", pad=2)
                ax.set_ylabel(RATE_LABEL[metric])
                y_grid(ax)
                if floor > 0:
                    axis_break(ax)
            axes[-1, 0].set_xticks(range(len(depths)))
            axes[-1, 0].set_xticklabels(["%d×" % d for d in depths])
            axes[-1, 0].set_xlabel("Coverage")
            axes[0, 0].legend(handles=[Patch(color=st.PLATFORM_COLOR[p], label=st.PLATFORM_LABEL[p])
                                       for p in platforms],
                              loc="lower center", bbox_to_anchor=(0.5, 1.12), ncol=len(platforms))
            fig.suptitle("%s, %s" % (RATE_LABEL[metric].replace(" (%)", ""), st.REFERENCE_LABEL[ref]))
            made.append(save(fig, ctx, "figR6_%s_%s" % (metric, ref)))
    if not made:
        raise NoData("no rows in summary_auto.tsv for %s" % ctx.caller)
    return ", ".join(made)


# ----------------------------------------- R7: genome-wide -> difficult regions

def fig_r7(ctx):
    auto = [r for r in table(path_of(ctx, "metrics", "summary_auto.tsv"))
            if r.get("caller") == ctx.caller and int(float(r.get("depth", 0))) == ctx.depth]
    strata = [r for r in table(path_of(ctx, "metrics_strata_auto", "summary_phase_strata.tsv"))
              if r.get("caller") == ctx.caller and int(float(r.get("depth", 0))) == ctx.depth
              and r.get("stratum") == ctx.stratum]
    if not auto or not strata:
        raise NoData("need summary_auto.tsv and summary_phase_strata.tsv (%s) at %dx" % (ctx.stratum, ctx.depth))

    def point(r):
        tt = num(r["truth_total"])
        return 100 * num(r["switch_rate"]), (100 * num(r["covered"]) / tt if tt else math.nan)

    start = {(r["reference"], r["phaser"], r["platform"]): point(r) for r in auto}
    end = {(r["reference"], r["phaser"], r["platform"]): point(r) for r in strata}
    made = []
    for phaser in ctx.phasers:
        keys = [k for k in start if k[1] == phaser and k in end]
        if not keys:
            continue
        refs = [r for r in st.REFERENCES if any(k[0] == r for k in keys)]
        xs = [v for k in keys for v in (start[k][0], end[k][0]) if v == v and v > 0]
        ys = [v for k in keys for v in (start[k][1], end[k][1]) if v == v]
        fig = plt.figure(figsize=(st.FULL_WIDTH * st.MM, 88 * st.MM), constrained_layout=True)
        grid = fig.add_gridspec(2, len(refs), height_ratios=[0.16, 1])
        legend_ax = fig.add_subplot(grid[0, :])
        legend_ax.axis("off")
        first = fig.add_subplot(grid[1, 0])
        axes = [[first] + [fig.add_subplot(grid[1, k], sharex=first, sharey=first) for k in range(1, len(refs))]]
        for ax in axes[0][1:]:
            ax.tick_params(labelleft=False)
        for ax, ref in zip(axes[0], refs):
            for p in st.PLATFORMS:
                k = (ref, phaser, p)
                if k not in start or k not in end:
                    continue
                (x0, y0), (x1, y1) = start[k], end[k]
                c, m = st.PLATFORM_COLOR[p], st.PLATFORM_MARKER[p]
                ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                            arrowprops=dict(arrowstyle="-|>", color=c, linestyle=(0, (3, 2)),
                                            linewidth=0.9, shrinkA=4, shrinkB=4, mutation_scale=8))
                ax.plot(x0, y0, marker=m, markersize=6, markerfacecolor="white", markeredgecolor=c,
                        markeredgewidth=1.2, linestyle="none", zorder=3)
                ax.plot(x1, y1, marker=m, markersize=6, color=c, markeredgecolor="white",
                        linestyle="none", zorder=3)
            ax.set_xscale("log")
            ax.set_title(st.REFERENCE_LABEL[ref])
            ax.set_xlabel("Switch error rate (%, lower to the right)")
            ax.grid(color=st.GRID, linewidth=0.5)
            ax.set_axisbelow(True)
        lo, hi = min(xs), max(xs)
        # reversed: fewer switches to the right, so the top right corner is best
        axes[0][0].set_xlim(hi * 1.6, lo / 1.6)
        log_ticks(axes[0][0].xaxis, lo / 1.6, hi * 1.6)
        ylo, yhi = min(ys), max(ys)
        pad = 0.08 * (yhi - ylo or 1)
        axes[0][0].set_ylim(ylo - pad, min(100.5, yhi + pad))
        axes[0][0].set_ylabel("Truth het SNPs phased (%)")
        handles = [series_handle(st.PLATFORM_COLOR[p], st.PLATFORM_MARKER[p], st.PLATFORM_LABEL[p], line=False)
                   for p in st.PLATFORMS if any(k[2] == p for k in keys)]
        handles += [Line2D([], [], marker="o", linestyle="none", markerfacecolor="white",
                           markeredgecolor=st.MUTED, markeredgewidth=1.2, label="Genome-wide (chr1-22)"),
                    Line2D([], [], marker="o", linestyle="none", color=st.MUTED, markeredgecolor="white",
                           label=stratum_label(ctx.stratum).replace("\n", " "))]
        top_legend(legend_ax, handles, ncol=3)
        fig.suptitle("%s: genome-wide to %s, %dx" % (st.PHASER_LABEL[phaser],
                                                     stratum_label(ctx.stratum).replace("\n", " ").lower(),
                                                     ctx.depth))
        made.append(save(fig, ctx, "figR7_difficult_drop_%s" % phaser))
    if not made:
        raise NoData("no phaser with both tables")
    return ", ".join(made)


# ------------------------------------------- R8: callers per genomic region

STRATA_R8 = [
    ("notinalldifficultregions", "Easy regions"),
    ("alldifficultregions", "All difficult"),
    ("SimpleRepeat_homopolymer_4to6_slop5", "Homopol 4–6 bp"),
    ("SimpleRepeat_homopolymer_7to11_slop5", "Homopol 7–11 bp"),
    ("SimpleRepeat_homopolymer_ge12_slop5", "Homopol ≥12 bp"),
    ("SimpleRepeat_imperfecthomopolge11_slop5", "Imperfect homopol ≥11 bp"),
    ("AllTandemRepeats_le50bp_slop5", "TR ≤50 bp"),
    ("AllTandemRepeats_51to200bp_slop5", "TR 51–200 bp"),
    ("AllTandemRepeats_201to10000bp_slop5", "TR 201 bp–10 kb"),
    ("AllTandemRepeats_ge10001bp_slop5", "TR >10 kb"),
    ("AllTandemRepeatsandHomopolymers_slop5", "TR and homopol"),
    ("segdups", "SegDups"),
    ("segdups_gt10kb", "SegDups >10 kb"),
    ("lowmappabilityall", "Low mappability"),
    ("satellites_slop5", "Satellites"),
    ("MHC", "MHC"),
    ("KIR", "KIR"),
    ("refseq_cds", "CDS"),
]
METRIC_LABEL = {"f1": "F1 score (%)", "precision": "Precision (%)", "recall": "Recall (%)"}


def fig_r8(ctx):
    rows = [r for r in table(path_of(ctx, "metrics_strata_auto", "summary_calls_strata.tsv"))
            if r.get("reference") == ctx.ref and int(float(r.get("depth", 0))) == ctx.depth]
    if not rows:
        raise NoData("no caller rows for %s at %dx in summary_calls_strata.tsv" % (ctx.ref, ctx.depth))
    val = {(r["platform"], r["caller"], r["stratum"]): 100 * num(r[ctx.metric]) for r in rows}
    callers = [c for c in st.CALLERS if c != "truth" and any(k[1] == c for k in val)]
    platforms = [p for p in st.PLATFORMS if any(k[0] == p for k in val)]
    strata = [(s, l) for s, l in STRATA_R8 if any(k[2] == s for k in val)]
    fig, axes = plt.subplots(len(platforms), 1, squeeze=False, sharex=True, sharey=True,
                             figsize=(st.FULL_WIDTH * st.MM, (30 + 34 * len(platforms)) * st.MM),
                             constrained_layout=True)
    width = 0.8 / len(callers)
    for ax, p in zip(axes[:, 0], platforms):
        for i, (s, _) in enumerate(strata):
            for j, c in enumerate(callers):
                y = val.get((p, c, s), math.nan)
                if y == y:
                    ax.bar(i + (j - (len(callers) - 1) / 2) * width, y, width=width * 0.92,
                           color=st.CALLER_COLOR[c], linewidth=0, zorder=2)
        ax.set_ylim(0, 100)
        ax.set_title(st.PLATFORM_LABEL[p], loc="left", pad=2)
        ax.set_ylabel(METRIC_LABEL[ctx.metric])
        y_grid(ax)
        ax.set_xlim(-0.6, len(strata) - 0.4)
    axes[-1, 0].set_xticks(range(len(strata)))
    axes[-1, 0].set_xticklabels([l for _, l in strata], rotation=40, ha="right", rotation_mode="anchor")
    axes[0, 0].legend(handles=[Patch(color=st.CALLER_COLOR[c], label=st.CALLER_LABEL[c]) for c in callers],
                      loc="lower center", bbox_to_anchor=(0.5, 1.15), ncol=len(callers))
    fig.suptitle("Het SNP calling across genomic regions, %s, %dx" % (st.REFERENCE_LABEL[ctx.ref], ctx.depth))
    return save(fig, ctx, "figR8_callers_strata_%s_%s" % (ctx.metric, ctx.ref))


# ----------------------------------------------- R9: reads rescued by CHM13

RESCUED = "#2F6F9F"
NOT_RESCUED = "#7A7A7A"


def fig_r9(ctx):
    made = []
    for platform in st.PLATFORMS:
        prefix = path_of(ctx, "metrics_reads", "rescue", "%s_%dx_%s" % (platform, ctx.depth, ctx.caller))
        combos = table(prefix + ".upset.tsv")
        if not combos:
            continue
        summary = table(prefix + ".summary.tsv")
        phasers = [c[len("untagged_GRCh38_"):] for c in combos[0] if c.startswith("untagged_GRCh38_")]
        sets = [("untagged_GRCh38_" + p, "%s untagged, GRCh38" % st.PHASER_LABEL.get(p, p), p, 0.45)
                for p in phasers]
        sets += [("tagged_T2T_" + p, "%s tagged, CHM13" % st.PHASER_LABEL.get(p, p), p, 1.0) for p in phasers]
        cols = []
        for r in combos:
            member = tuple(r[s[0]] == "1" for s in sets)
            rescued = any(r["untagged_GRCh38_" + p] == "1" and r["tagged_T2T_" + p] == "1" for p in phasers)
            cols.append((rescued, int(r["reads"]), member))
        groups = [(False, "Untagged on GRCh38, not rescued"), (True, "Rescued by CHM13")]
        shown = []
        for flag, _ in groups:
            g = sorted([c for c in cols if c[0] == flag], key=lambda c: -c[1])
            shown.append(g[:ctx.top])
        order = shown[0] + shown[1]
        if not order:
            continue
        set_size = [sum(c[1] for c in cols if c[2][i]) for i in range(len(sets))]

        fig = plt.figure(figsize=(st.FULL_WIDTH * st.MM, 100 * st.MM))
        # set sizes | set names | intersections over the membership matrix
        gs = fig.add_gridspec(2, 3, width_ratios=[0.17, 0.38, 1], height_ratios=[1.3, 0.75],
                              left=0.02, right=0.99, top=0.9, bottom=0.19, wspace=0.0, hspace=0.05)
        bar_ax = fig.add_subplot(gs[0, 2])
        mat_ax = fig.add_subplot(gs[1, 2], sharex=bar_ax)
        set_ax = fig.add_subplot(gs[1, 0], sharey=mat_ax)
        name_ax = fig.add_subplot(gs[1, 1], sharey=mat_ax)
        name_ax.axis("off")

        # the intersection bars, a gap between the groups
        xs, x = [], 0.0
        for gi, g in enumerate(shown):
            for _ in g:
                xs.append(x)
                x += 1
            x += 0.6
        heights = [c[1] for c in order]
        colors = [RESCUED if c[0] else NOT_RESCUED for c in order]
        bar_ax.bar(xs, heights, width=0.7, color=colors, linewidth=0, zorder=2)
        top = max(heights)
        for xi, h in zip(xs, heights):
            bar_ax.text(xi, h + 0.015 * top, thousands(h), ha="center", va="bottom",
                        fontsize=st.FONT_SIZE - 1, rotation=90 if len(order) > 12 else 0)
        bar_ax.set_ylim(0, top * 1.32)
        bar_ax.set_ylabel("Reads")
        bar_ax.tick_params(labelbottom=False, bottom=False)
        bar_ax.spines["bottom"].set_visible(False)
        y_grid(bar_ax)
        bar_ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: short(v) if v else "0"))
        # brackets over the groups
        i0 = 0
        for (flag, title), g in zip(groups, shown):
            if not g:
                continue
            a, b = xs[i0], xs[i0 + len(g) - 1]
            n = sum(c[1] for c in cols if c[0] == flag)
            yb = top * 1.14
            bar_ax.plot([a - 0.3, a - 0.3, b + 0.3, b + 0.3], [yb - 0.03 * top, yb, yb, yb - 0.03 * top],
                        color=st.INK, linewidth=0.6)
            bar_ax.text((a + b) / 2, yb + 0.02 * top, "%s\nn = %s" % (title, thousands(n)),
                        ha="center", va="bottom")
            i0 += len(g)

        # the membership matrix
        ny = len(sets)
        for i in range(ny):
            if i % 2 == 0:
                mat_ax.axhspan(i - 0.5, i + 0.5, color="#F2F2F2", zorder=0)
        for xi, c in zip(xs, order):
            on = [i for i in range(ny) if c[2][i]]
            mat_ax.scatter([xi] * ny, range(ny), s=22, color="#D6D6D6", zorder=2, linewidths=0)
            if len(on) > 1:
                mat_ax.plot([xi, xi], [min(on), max(on)], color=st.INK, linewidth=1.2, zorder=3)
            for i in on:
                col = st.PHASER_COLOR.get(sets[i][2], st.INK) if len(on) == 1 else st.INK
                mat_ax.scatter([xi], [i], s=22, color=col, zorder=4, linewidths=0)
        mat_ax.set_ylim(ny - 0.5, -0.5)
        mat_ax.set_xlim(xs[0] - 0.6, xs[-1] + 0.6)
        mat_ax.axis("off")

        # set sizes on the left, the labels between
        for i, (key, label, p, alpha) in enumerate(sets):
            set_ax.barh(i, set_size[i], height=0.55, color=st.PHASER_COLOR.get(p, st.INK), alpha=alpha,
                        linewidth=0)
            name_ax.text(0.97, i, label, transform=name_ax.get_yaxis_transform(), ha="right", va="center")
        set_ax.invert_xaxis()
        set_ax.set_yticks([])
        set_ax.spines["left"].set_visible(False)
        set_ax.spines["right"].set_visible(True)
        set_ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: short(v) if v else "0"))
        set_ax.tick_params(axis="x", labelsize=st.FONT_SIZE - 1)
        for lab in set_ax.get_xticklabels():
            lab.set_rotation(40)
            lab.set_ha("right")
        set_ax.set_xlabel("Reads in set")

        # what the rescued reads are, from the summary
        notes = []
        for p in phasers:
            def get(direction, cat, name, field="reads"):
                for r in summary:
                    if r["phaser"] == p and r["direction"] == direction and r["category"] == cat \
                            and r["name"] == name:
                        return num(r[field])
                return math.nan
            n = get("rescued", "total", "reads")
            hard = get("rescued", "stratum_on_T2T", "alldifficultregions", "share")
            bg = get("rescued", "stratum_on_T2T", "alldifficultregions", "background_share")
            ok = get("rescued", "tag_on_T2T", "correct", "share")
            lost = get("lost", "total", "reads")
            notes.append("%s: %s rescued, %.1f%% in difficult regions (all reads %.1f%%), %.1f%% tagged "
                         "correctly; %s tagged on GRCh38 only" % (
                             st.PHASER_LABEL.get(p, p), thousands(n), 100 * hard, 100 * bg, 100 * ok,
                             thousands(lost)))
        fig.text(0.5, 0.015, "\n".join(notes), ha="center", va="bottom", color=st.MUTED,
                 fontsize=st.FONT_SIZE - 1)
        fig.suptitle("Reads CHM13 tags that GRCh38 does not: %s, %dx, %s input" % (
            st.PLATFORM_LABEL[platform], ctx.depth, st.CALLER_LABEL.get(ctx.caller, ctx.caller)), y=0.985)
        made.append(save(fig, ctx, "figR9_rescued_%s" % platform))
    if not made:
        raise NoData("no metrics_reads/rescue/*.upset.tsv - run bin/rescue_reads.sh")
    return ", ".join(made)


# --------------------------------------- R10: difficult regions, GRCh38 vs T2T

STRATA_R10 = [
    ("notinalldifficultregions", "Easy regions"),
    ("alldifficultregions", "All difficult"),
    ("segdups", "SegDups"),
    ("lowmappabilityall", "Low mappability"),
    ("AllTandemRepeatsandHomopolymers_slop5", "TR and homopol"),
    ("satellites_slop5", "Satellites"),
    ("MHC", "MHC"),
]
R10_METRICS = [
    ("switch", "Switch error rate (%)"),
    ("n50", "Block N50 (Mb)"),
    ("accuracy", "Haplotag accuracy (%)"),
]


def zoomed(values, lo_cap=None, hi_cap=None):
    """y limits around the values, not from 0: the references differ by little"""
    v = [x for x in values if x == x]
    if not v:
        return 0, 1
    lo, hi = min(v), max(v)
    pad = 0.12 * (hi - lo or abs(hi) * 0.1 or 1)
    lo, hi = lo - pad, hi + pad
    if lo_cap is not None:
        lo = max(lo, lo_cap)
    if hi_cap is not None:
        hi = min(hi, hi_cap)
    return lo, hi


def fig_r10(ctx):
    strata_rows = [r for r in table(path_of(ctx, "metrics_strata_auto", "summary_phase_strata.tsv"))
                   if r.get("caller") == ctx.caller and int(float(r.get("depth", 0))) == ctx.depth
                   and r.get("phaser") == ctx.phaser]
    class_rows = [r for r in table(path_of(ctx, "metrics_reads", "summary_classes.tsv"))
                  if r.get("caller") == ctx.caller and r.get("depth") == str(ctx.depth)
                  and r.get("phaser") == ctx.phaser and r.get("group_type") == "all"]
    if not strata_rows:
        raise NoData("no %s %s %dx rows in summary_phase_strata.tsv" % (ctx.caller, ctx.phaser, ctx.depth))
    val = {}
    for r in strata_rows:
        k = (r["reference"], r["platform"], r["stratum"])
        val[("switch",) + k] = 100 * num(r["switch_rate"])
        val[("n50",) + k] = num(r["block_n50"]) / 1e6
    for r in class_rows:
        k = (r["reference"], r["platform"], r["stratum"])
        comp = num(r["compared"])
        if comp:
            val[("accuracy",) + k] = 100 * num(r["correct"]) / comp
    refs = [r for r in st.REFERENCES if any(k[1] == r for k in val)]
    platforms = [p for p in st.PLATFORMS if any(k[2] == p for k in val)]
    strata = [(s, l) for s, l in STRATA_R10 if any(k[3] == s for k in val)]
    metrics = [(m, l) for m, l in R10_METRICS if any(k[0] == m for k in val)]

    fig = plt.figure(figsize=(st.FULL_WIDTH * st.MM, (30 + 40 * len(metrics)) * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(len(metrics) + 1, len(platforms), height_ratios=[0.1] + [1] * len(metrics))
    legend_ax = fig.add_subplot(grid[0, :])
    legend_ax.axis("off")
    width = 0.8 / len(refs)
    for i, (m, label) in enumerate(metrics):
        lo, hi = zoomed([v for k, v in val.items() if k[0] == m],
                        lo_cap=0, hi_cap=100 if m != "n50" else None)
        first = None
        for j, p in enumerate(platforms):
            ax = fig.add_subplot(grid[1 + i, j], sharey=first)
            first = first or ax
            for si, (s, _) in enumerate(strata):
                for ri, ref in enumerate(refs):
                    y = val.get((m, ref, p, s), math.nan)
                    if y != y:
                        continue
                    x = si + (ri - (len(refs) - 1) / 2) * width
                    grch = ref == "GRCh38"
                    ax.bar(x, y - lo, bottom=lo, width=width * 0.9, linewidth=0.6, zorder=2,
                           color="white" if grch else st.PLATFORM_COLOR[p],
                           edgecolor=st.PLATFORM_COLOR[p], hatch="////" if grch else None)
            ax.set_ylim(lo, hi)
            if lo > 0:
                axis_break(ax)
            ax.set_xticks(range(len(strata)))
            if i == len(metrics) - 1:
                ax.set_xticklabels([l for _, l in strata], rotation=50, ha="right", rotation_mode="anchor")
            else:
                ax.set_xticklabels([])
            if j:
                ax.tick_params(labelleft=False)
            else:
                ax.set_ylabel(label)
            if i == 0:
                ax.set_title(st.PLATFORM_LABEL[p])
            y_grid(ax)
            ax.set_xlim(-0.6, len(strata) - 0.4)
    handles = [Patch(facecolor="white", edgecolor=st.INK, hatch="////", label="GRCh38"),
               Patch(facecolor=st.MUTED, edgecolor=st.MUTED, label="T2T-CHM13 (platform colour)")]
    top_legend(legend_ax, handles)
    fig.suptitle("%s, %s input, %dx: genomic regions on GRCh38 and T2T-CHM13" % (
        st.PHASER_LABEL[ctx.phaser], st.CALLER_LABEL.get(ctx.caller, ctx.caller), ctx.depth))
    return save(fig, ctx, "figR10_regions_refs_%s" % ctx.phaser)


# ------------------------------------- R12: phasing in easy and hard regions

R12_METRICS = [
    ("switch", "Switch error\nrate (%)"),
    ("hamming", "Hamming error\nrate (%)"),
    ("n50", "Block N50 (Mb)"),
    ("accuracy", "Haplotag\naccuracy (%)"),
]


def fig_r12(ctx):
    strata_rows = [r for r in table(path_of(ctx, "metrics_strata_auto", "summary_phase_strata.tsv"))
                   if r.get("caller") == ctx.caller and int(float(r.get("depth", 0))) == ctx.depth
                   and r.get("phaser") == ctx.phaser]
    class_rows = [r for r in table(path_of(ctx, "metrics_reads", "summary_classes.tsv"))
                  if r.get("caller") == ctx.caller and r.get("depth") == str(ctx.depth)
                  and r.get("phaser") == ctx.phaser and r.get("group_type") == "all"]
    if not strata_rows:
        raise NoData("no %s %s %dx rows in summary_phase_strata.tsv" % (ctx.caller, ctx.phaser, ctx.depth))
    val = {}
    for r in strata_rows:
        k = (r["reference"], r["platform"], r["stratum"])
        val[("switch",) + k] = 100 * num(r["switch_rate"])
        val[("hamming",) + k] = 100 * num(r.get("hamming_rate"))
        val[("n50",) + k] = num(r["block_n50"]) / 1e6
    for r in class_rows:
        comp = num(r["compared"])
        if comp:
            val[("accuracy", r["reference"], r["platform"], r["stratum"])] = 100 * num(r["correct"]) / comp
    made = []
    for ref in st.REFERENCES:
        keys = [k for k in val if k[1] == ref and val[k] == val[k]]
        if not keys:
            continue
        platforms = [p for p in st.PLATFORMS if any(k[2] == p for k in keys)]
        strata = [(s_, l) for s_, l in STRATA_R8 if any(k[3] == s_ for k in keys)]
        metrics = [(m, l) for m, l in R12_METRICS if any(k[0] == m for k in keys)]
        fig, axes = plt.subplots(len(metrics), 1, squeeze=False, sharex=True,
                                 figsize=(st.FULL_WIDTH * st.MM, (34 + 32 * len(metrics)) * st.MM),
                                 constrained_layout=True)
        width = 0.8 / len(platforms)
        for ax, (m, label) in zip(axes[:, 0], metrics):
            drawn = []
            for si, (s_, _) in enumerate(strata):
                for j, p in enumerate(platforms):
                    y = val.get((m, ref, p, s_), math.nan)
                    if y != y:
                        continue
                    drawn.append(y)
                    ax.bar(si + (j - (len(platforms) - 1) / 2) * width, y, width=width * 0.9,
                           color=st.PLATFORM_COLOR[p], linewidth=0, zorder=2)
            if m == "accuracy" and drawn:
                lo, hi = zoomed(drawn, lo_cap=0, hi_cap=100)
                ax.set_ylim(lo, hi)
                for si, (s_, _) in enumerate(strata):
                    if not any((m, ref, p, s_) in val for p in platforms):
                        ax.text(si, lo + 0.05 * (hi - lo), "n/a", ha="center", color=st.MUTED,
                                fontsize=st.FONT_SIZE - 1)
                if lo > 0:
                    axis_break(ax)
            ax.set_ylabel(label)
            y_grid(ax)
            # easy regions apart from the difficult ones
            ax.axvline(0.5, color=st.MUTED, linestyle=":", linewidth=0.6)
        ax = axes[-1, 0]
        ax.set_xticks(range(len(strata)))
        ax.set_xticklabels([l for _, l in strata], rotation=40, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.6, len(strata) - 0.4)
        axes[0, 0].legend(handles=[Patch(color=st.PLATFORM_COLOR[p], label=st.PLATFORM_LABEL[p])
                                   for p in platforms],
                          loc="lower center", bbox_to_anchor=(0.5, 1.08), ncol=len(platforms))
        fig.suptitle("Phasing across genomic regions: %s, %s, %s input, %dx" % (
            st.REFERENCE_LABEL[ref], st.PHASER_LABEL[ctx.phaser],
            st.CALLER_LABEL.get(ctx.caller, ctx.caller), ctx.depth))
        made.append(save(fig, ctx, "figR12_regions_phasing_%s_%s" % (ctx.phaser, ref)))
    return ", ".join(made)


# ------------------------------------------------ R13: what each region is

STRATUM_BAR = "#4A5C70"


def fig_r13(ctx):
    made = []
    for ref in st.REFERENCES:
        stats = {r["stratum"]: r for r in table(path_of(ctx, "metrics_strata_auto",
                                                         "regions_statistics_%s.tsv" % ref))}
        mech = {}
        for r in table(path_of(ctx, "metrics_strata_auto", "mechanism_%s.tsv" % ref)):
            if r.get("caller") == "-" and r.get("depth") == str(ctx.depth):
                mech[(r["platform"], r["stratum"], r["metric"])] = num(r["value"])
        if not stats:
            continue
        strata = [(s_, l) for s_, l in STRATA_R10 if s_ in stats]
        platforms = [p for p in st.PLATFORMS if any(k[0] == p for k in mech)]

        def stat(s_, col):
            return num(stats[s_].get(col))

        def share(s_, col):
            n = stat(s_, "src_variants")
            return 100 * stat(s_, col) / n if n else math.nan

        # (label, per stratum value or None, per platform metric or None, log, scale)
        rows = [
            ("Size in the\nbenchmark (Mb)", lambda s_: stat(s_, "bp") / 1e6, None, True),
            ("Truth het SNPs", lambda s_: stat(s_, "truth_het_snps"), None, True),
            ("Het SNPs per Mb", lambda s_: stat(s_, "het_snps_per_mb"), None, False),
            ("Median distance\nbetween het SNPs (bp)", lambda s_: stat(s_, "median_inter_snp_dist"), None, True),
            ("Indels (% of\ntruth variants)", lambda s_: share(s_, "src_indel"), None, False),
            ("1/2 genotypes (%\nof truth variants)", lambda s_: share(s_, "src_1_2"), None, False),
            ("Mean |VAF − 0.5|", None, ("vaf_dev", 1), False),
            ("Neighbouring SNP pairs\nwith >2 haplotypes (%)", None, ("paralog_link_frac", 100), False),
            ("Depth / genome", None, ("depth_ratio", 1), False),
            ("Depth with\nMAPQ < 20 (%)", None, ("lowmapq_frac", 100), False),
        ]
        rows = [r for r in rows if r[1] is not None or (platforms and any(
            (p, s_, r[2][0]) in mech for p in platforms for s_, _ in strata))]
        fig, axes = plt.subplots(len(rows), 1, squeeze=False, sharex=True,
                                 figsize=(st.FULL_WIDTH * st.MM, (30 + 21 * len(rows)) * st.MM),
                                 constrained_layout=True)
        width = 0.8 / max(1, len(platforms))
        for ax, (label, per_stratum, per_platform, log) in zip(axes[:, 0], rows):
            values = []
            if per_stratum is not None:
                for si, (s_, _) in enumerate(strata):
                    y = per_stratum(s_)
                    if y == y and (y > 0 or not log):
                        values.append(y)
                        ax.bar(si, y, width=0.55, color=STRATUM_BAR, linewidth=0, zorder=2)
                        ax.annotate(tick_label(y), (si, y), xytext=(0, 1.5), textcoords="offset points",
                                    ha="center", va="bottom", fontsize=st.FONT_SIZE - 2)
            else:
                metric, scale = per_platform
                for si, (s_, _) in enumerate(strata):
                    for j, p in enumerate(platforms):
                        y = mech.get((p, s_, metric), math.nan) * scale
                        if y == y:
                            ax.bar(si + (j - (len(platforms) - 1) / 2) * width, y, width=width * 0.9,
                                   color=st.PLATFORM_COLOR[p], linewidth=0, zorder=2)
                if metric == "depth_ratio":
                    ax.axhline(1, color=st.MUTED, linestyle="--", linewidth=0.6)
            if log and values:
                # plain numbers on the log axis, not 10^x
                ax.set_yscale("log")
                lo, hi = min(values) / 2, max(values) * 4
                ax.set_ylim(lo, hi)
                log_ticks(ax.yaxis, lo, hi)
            else:
                ax.margins(y=0.25)
            ax.set_ylabel(label, fontsize=st.FONT_SIZE - 1)
            y_grid(ax)
            ax.axvline(0.5, color=st.MUTED, linestyle=":", linewidth=0.6)
        ax = axes[-1, 0]
        ax.set_xticks(range(len(strata)))
        ax.set_xticklabels([l for _, l in strata], rotation=30, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.6, len(strata) - 0.4)
        if platforms:
            axes[0, 0].legend(handles=[Patch(color=STRATUM_BAR, label="Region (same for every platform)")] +
                              [Patch(color=st.PLATFORM_COLOR[p], label=st.PLATFORM_LABEL[p]) for p in platforms],
                              loc="lower center", bbox_to_anchor=(0.5, 1.15), ncol=3)
        fig.suptitle("The genomic regions on %s: size, variants and what the reads see (%dx)" % (
            st.REFERENCE_LABEL[ref], ctx.depth))
        made.append(save(fig, ctx, "figR13_regions_profile_%s" % ref))
    if not made:
        raise NoData("no regions_statistics_<ref>.tsv - run regions_statistics.sh")
    return ", ".join(made)


def tick_label(v):
    if v >= 1e6:
        return "%.1fM" % (v / 1e6)
    if v >= 1e4:
        return "%.0fk" % (v / 1e3)
    if v >= 100:
        return "%.0f" % v
    if v >= 10:
        return "%.1f" % v
    return "%.2g" % v


# -------------------------------------- R14: callers' input, by region

# the platform names on two lines: they label groups of three or four rows
PLATFORM_TWO_LINES = {"ont4k": "ONT\n4 kHz", "ont5k": "ONT\n5 kHz", "hifi": "PacBio\nSequel II",
                      "revio": "PacBio\nRevio"}

def fig_r14(ctx):
    """Heatmaps of phasing on each caller's het SNPs, per genomic region:
    the truth het SNPs of the region phased correctly (correct / truth_total,
    the same denominator for every input, so callers compare), and how many
    points of it each caller loses against the truth input. Rows: platform x
    input, columns: regions."""
    rows = [r for r in table(path_of(ctx, "metrics_strata_auto", "summary_phase_strata.tsv"))
            if int(float(r.get("depth", 0))) == ctx.depth and r.get("phaser") == ctx.phaser]
    if not rows:
        raise NoData("no %s %dx rows in summary_phase_strata.tsv" % (ctx.phaser, ctx.depth))
    val = {}
    for r in rows:
        tt = num(r.get("truth_total"))
        v = 100 * num(r.get("correct")) / tt if tt else 100 * num(r.get("recall_phased"))
        val[(r["reference"], r["caller"], r["platform"], r["stratum"])] = v
    inputs = [c for c in st.CALLERS if any(k[1] == c for k in val)]
    callers = [c for c in inputs if c != "truth"]
    if not callers:
        raise NoData("no caller input in summary_phase_strata.tsv")
    made = []
    for ref in st.REFERENCES:
        keys = [k for k in val if k[0] == ref]
        if not keys:
            continue
        platforms = [p for p in st.PLATFORMS if any(k[2] == p for k in keys)]
        strata = [(s_, l) for s_, l in STRATA_R8 if any(k[3] == s_ for k in keys)]
        # darker = more: more phased correctly / more lost to the caller
        for tag, row_inputs, value, bar_label in (
                ("phased", inputs, lambda p, c, s_: val.get((ref, c, p, s_), math.nan),
                 "Truth het SNPs phased correctly (%)"),
                ("loss", callers, lambda p, c, s_: val.get((ref, "truth", p, s_), math.nan)
                 - val.get((ref, c, p, s_), math.nan),
                 "Points lost against the truth input")):
            line = [(p, c) for p in platforms for c in row_inputs]
            matrix = [[value(p, c, s_) for s_, _ in strata] for p, c in line]
            shown = [v for row in matrix for v in row if v == v]
            if not shown:
                continue
            fig, ax = plt.subplots(figsize=(st.FULL_WIDTH * st.MM, (30 + 5.2 * len(line)) * st.MM),
                                   constrained_layout=True)
            im = heatmap(ax, matrix, [st.CALLER_LABEL[c] for _, c in line], [l for _, l in strata],
                         [(PLATFORM_TWO_LINES[p], len(row_inputs)) for p in platforms],
                         [("", len(strata))], label=lambda v: "%.1f" % v, cmap=st.SEQUENTIAL,
                         vmin=min(shown), vmax=max(shown))
            plt.setp(ax.get_xticklabels(), rotation=40, ha="right", rotation_mode="anchor")
            colorbar(fig, im, ax, bar_label)
            fig.suptitle("%s: %s, %s, %dx, on each input's het SNPs" % (
                bar_label, st.REFERENCE_LABEL[ref], st.PHASER_LABEL[ctx.phaser], ctx.depth))
            made.append(save(fig, ctx, "figR14_callers_%s_%s_%s" % (tag, ctx.phaser, ref)))
    if not made:
        raise NoData("nothing to draw")
    return ", ".join(made)


# --------------------------- R15: read accuracy, grid of input x phaser

def fig_r15(ctx):
    """Read-level haplotag accuracy with every block oriented by its SNPs
    (readacc_run.sh, chr1-22, blocks of >= 2 SNPs): precision = reads on the
    truth's haplotype / reads tagged by both, recall = the same / reads the
    truth tags. Rows: platform x input, columns: phasers; a PDF per metric and
    reference."""
    rows = [r for r in table(path_of(ctx, "metrics", "summary_readacc_auto_min2.tsv"))
            if int(float(r.get("depth", 0))) == ctx.depth]
    if not rows:
        raise NoData("no %dx rows in metrics/summary_readacc_auto_min2.tsv - run readacc_run.sh" % ctx.depth)
    val = {}
    for r in rows:
        k = (r["reference"], r["caller"], r["platform"], r["phaser"])
        val[("precision",) + k] = 100 * num(r["read_precision"])
        val[("recall",) + k] = 100 * num(r["read_recall"])
    inputs = [c for c in st.CALLERS if any(k[2] == c for k in val)]
    made = []
    for metric, label in (("precision", "Read haplotag precision (%)"), ("recall", "Read haplotag recall (%)")):
        for ref in st.REFERENCES:
            keys = [k for k in val if k[0] == metric and k[1] == ref]
            if not keys:
                continue
            platforms = [p for p in st.PLATFORMS if any(k[3] == p for k in keys)]
            phasers = [p for p in st.PHASERS if any(k[4] == p for k in keys)]
            line = [(p, c) for p in platforms for c in inputs]
            matrix = [[val.get((metric, ref, c, p, ph), math.nan) for ph in phasers] for p, c in line]
            shown = [v for row in matrix for v in row if v == v]
            if not shown:
                continue
            fig, ax = plt.subplots(figsize=(110 * st.MM, (30 + 5.2 * len(line)) * st.MM),
                                   constrained_layout=True)
            im = heatmap(ax, matrix, [st.CALLER_LABEL[c] for _, c in line],
                         [st.PHASER_LABEL[ph] for ph in phasers],
                         [(PLATFORM_TWO_LINES[p], len(inputs)) for p in platforms],
                         [("", len(phasers))], label=lambda v: "%.2f" % v, cmap=st.SEQUENTIAL,
                         vmin=min(shown), vmax=max(shown))
            plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")
            colorbar(fig, im, ax, label)
            fig.suptitle("%s, %s, %dx, blocks oriented by their SNPs" % (
                label.replace(" (%)", ""), st.REFERENCE_LABEL[ref], ctx.depth))
            made.append(save(fig, ctx, "figR15_read_%s_%dx_%s" % (metric, ctx.depth, ref)))
    if not made:
        raise NoData("nothing to draw")
    return ", ".join(made)


# ------------- R16: reads on the wrong haplotype per region, every input

ENRICH = LinearSegmentedColormap.from_list("enrich", ["#2E6F9E", "#F7F7F7", "#B5523B"])
MIN_WRONG = 10


def fig_r16(ctx):
    """Tagged reads on the wrong haplotype (1 - precision, blocks oriented by
    their SNPs, bin/read_classes.sh) per tier-1 region, for every phaser input:
    the truth and each caller. Rows: platform x phaser, columns: regions (a
    read is in a region if half its alignment or more lies there).
      figR16_wrong_<input>_<ref>   the percentage
      figR16_ratio_<input>         T2T-CHM13 / GRCh38; a cell with fewer than
                                   MIN_WRONG wrong reads on either reference is
                                   left out (a ratio of a few reads says nothing)"""
    rows = [r for r in table(path_of(ctx, "metrics_reads", "summary_classes.tsv"))
            if r.get("depth") == str(ctx.depth) and r.get("group_type") == "class"
            and r.get("group") == "tagged"]
    if not rows:
        raise NoData("no tagged-read rows at %dx in summary_classes.tsv" % ctx.depth)
    data = {}
    for r in rows:
        comp, corr = num(r["compared"]), num(r["correct"])
        data[(r["reference"], r["caller"], r["platform"], r["phaser"], r["stratum"])] = (comp, corr)
    order = ["genome"] + [s_ for s_, _ in STRATA_R10]
    strata = [s_ for s_ in order if any(k[4] == s_ for k in data)]
    labels = ["Whole genome\n(chr1-22)" if s_ == "genome" else stratum_label(s_) for s_ in strata]
    inputs = [c for c in st.CALLERS if any(k[1] == c for k in data)]

    def wrong(ref, c, p, ph, s_):
        v = data.get((ref, c, p, ph, s_))
        return 100 * (1 - v[1] / v[0]) if v and v[0] else math.nan

    def rows_of(c):
        return [(p, ph) for p in st.PLATFORMS for ph in st.PHASERS
                if any(k[1] == c and k[2] == p and k[3] == ph for k in data)]

    def draw(name, line, matrix, bar_label, title, **image):
        fig, ax = plt.subplots(figsize=(st.FULL_WIDTH * st.MM, (34 + 5 * len(line)) * st.MM),
                               constrained_layout=True)
        groups = [(PLATFORM_TWO_LINES[p], sum(1 for q, _ in line if q == p))
                  for p in st.PLATFORMS if any(q == p for q, _ in line)]
        im = heatmap(ax, matrix, [st.PHASER_LABEL[ph] for _, ph in line], labels, groups,
                     [("", len(strata))], label=lambda v: "%.2f" % v, **image)
        plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
        bar = colorbar(fig, im, ax, bar_label)
        if "norm" in image:
            # a ratio: a few plain ticks on both sides of 1
            lo, hi = image["norm"].vmin, image["norm"].vmax
            ticks = [t for t in (0.1, 0.2, 0.33, 0.5, 0.67, 0.8, 1, 1.25, 1.5, 2, 3, 5, 10) if lo <= t <= hi]
            bar.set_ticks(ticks)
            bar.set_ticklabels(["%g" % t for t in ticks])
        fig.suptitle(title)
        return save(fig, ctx, name)

    made = []
    # the percentage, one colour scale over every input of a reference
    for ref in st.REFERENCES:
        shown = [wrong(ref, c, p, ph, s_) for c in inputs for p, ph in rows_of(c) for s_ in strata]
        shown = [v for v in shown if v == v]
        if not shown:
            continue
        for c in inputs:
            line = rows_of(c)
            matrix = [[wrong(ref, c, p, ph, s_) for s_ in strata] for p, ph in line]
            if not any(v == v for row in matrix for v in row):
                continue
            made.append(draw("figR16_wrong_%s_%s" % (c, ref), line, matrix,
                             "Tagged reads on the wrong haplotype (%%), %d×" % ctx.depth,
                             "%s input, %s: tagged reads on the wrong haplotype" % (
                                 st.CALLER_LABEL[c], st.REFERENCE_LABEL[ref]),
                             cmap=st.SEQUENTIAL, vmin=min(shown), vmax=max(shown)))
    # T2T over GRCh38
    for c in inputs:
        line = rows_of(c)
        matrix = []
        for p, ph in line:
            row = []
            for s_ in strata:
                a, b = data.get(("T2T", c, p, ph, s_)), data.get(("GRCh38", c, p, ph, s_))
                ok = a and b and a[0] - a[1] >= MIN_WRONG and b[0] - b[1] >= MIN_WRONG
                ra, rb = wrong("T2T", c, p, ph, s_), wrong("GRCh38", c, p, ph, s_)
                row.append(ra / rb if ok and rb else math.nan)
            matrix.append(row)
        shown = [v for r_ in matrix for v in r_ if v == v and v > 0]
        if not shown:
            continue
        span = max(abs(math.log10(v)) for v in shown) or 1
        made.append(draw("figR16_ratio_%s" % c, line, matrix,
                         "Tagged reads on the wrong haplotype,\nT2T-CHM13 / GRCh38, %d×" % ctx.depth,
                         "%s input: T2T-CHM13 against GRCh38 (blank: under %d wrong reads)" % (
                             st.CALLER_LABEL[c], MIN_WRONG),
                         cmap=ENRICH, norm=LogNorm(vmin=10 ** -span, vmax=10 ** span)))
    if not made:
        raise NoData("nothing to draw")
    return ", ".join(made)



# ------------------- R17: callers per region, precision / recall / F1

R17_METRICS = [("precision", "Precision", 0.45), ("recall", "Recall", 0.7), ("f1", "F1 score", 1.0)]


def tint(color, f):
    """color mixed with white: f = 1 the colour itself, smaller is lighter"""
    rgb = matplotlib_colors.to_rgb(color)
    return tuple(1 - f * (1 - c) for c in rgb)


def fig_r17(ctx):
    """Het SNP calls of each caller per genomic region, at one coverage: a
    panel per region, three bars per caller - precision, recall, F1 - in light
    to full shades of the caller's colour, the values written in the bars
    (metrics_strata_auto/summary_calls_strata.tsv). A PDF per platform and
    reference; the y axis is shared and zoomed, the break marks it."""
    rows = [r for r in table(path_of(ctx, "metrics_strata_auto", "summary_calls_strata.tsv"))
            if int(float(r.get("depth", 0))) == ctx.depth]
    if not rows:
        raise NoData("no %dx rows in summary_calls_strata.tsv" % ctx.depth)
    val = {}
    for r in rows:
        for m, _, _ in R17_METRICS:
            val[(r["reference"], r["platform"], r["caller"], r["stratum"], m)] = num(r[m])
    made = []
    # one PDF per reference, then one with T2T as bars and GRCh38 as a tick on
    # each bar (other = the reference drawn as ticks)
    views = [(ref, None) for ref in st.REFERENCES] + [("T2T", "GRCh38")]
    for ref, other in views:
        for platform in st.PLATFORMS:
            keys = [k for k in val if k[0] == ref and k[1] == platform and val[k] == val[k]]
            if not keys:
                continue
            if other:
                okeys = [k for k in val if k[0] == other and k[1] == platform and val[k] == val[k]]
                if not okeys:
                    continue
                keys += okeys
            callers = [c for c in st.CALLERS if c != "truth" and any(k[2] == c for k in keys)]
            strata = [(s_, l) for s_, l in STRATA_R8 if any(k[3] == s_ for k in keys)]
            lo, hi = zoomed([val[k] for k in keys], lo_cap=0, hi_cap=1)
            lo = math.floor(lo * 50) / 50
            hi = 1.0
            ncols = 6
            nrows = math.ceil(len(strata) / ncols)
            fig, axes = plt.subplots(nrows, ncols, squeeze=False, sharey=True,
                                     figsize=(st.FULL_WIDTH * st.MM, (30 + 50 * nrows) * st.MM),
                                     layout="constrained")
            nbar = len(R17_METRICS)
            width = 0.9 / nbar
            for n, (s_, label) in enumerate(strata):
                ax = axes[n // ncols][n % ncols]
                for ci, c in enumerate(callers):
                    for mi, (m, mlabel, f) in enumerate(R17_METRICS):
                        y = val.get((ref, platform, c, s_, m), math.nan)
                        if y != y:
                            continue
                        x = ci + (mi - (nbar - 1) / 2) * width
                        color = tint(st.CALLER_COLOR[c], f)
                        ax.bar(x, y - lo, bottom=lo, width=width * 0.94, color=color, linewidth=0, zorder=2)
                        # the metric written up from the foot of the bar, where it
                        # clears both the value at the top and the GRCh38 line: left
                        # out of a bar too short for it
                        # the height the word takes, about 2.7% of the axis per letter
                        room = lo + (0.015 + 0.027 * len(mlabel)) * (hi - lo)
                        g = val.get((other, platform, c, s_, m), math.nan) if other else math.nan
                        if y > room + 0.01 * (hi - lo) and not (g == g and g < room + 0.01 * (hi - lo)):
                            ax.text(x, lo + 0.012 * (hi - lo), mlabel, rotation=90, ha="center", va="bottom",
                                    fontsize=st.FONT_SIZE - 2.5, zorder=3,
                                    color="white" if f > 0.6 and c != "longcalld" else st.INK)
                        # the value inside a tall bar, over a short one; not on the
                        # comparison, where the GRCh38 lines would cross it
                        if other:
                            g = val.get((other, platform, c, s_, m), math.nan)
                            if g == g:
                                ax.plot([x - width * 0.47, x + width * 0.47], [g, g], color=st.INK,
                                        linewidth=1.3, solid_capstyle="butt", zorder=4)
                            continue
                        inside = y - lo > 0.3 * (hi - lo)
                        dark = inside and f > 0.6 and c != "longcalld"
                        ax.text(x, y - 0.01 * (hi - lo) if inside else y + 0.01 * (hi - lo), "%.3f" % y,
                                rotation=90, ha="center", va="top" if inside else "bottom",
                                fontsize=st.FONT_SIZE - 2.5, color="white" if dark else st.INK, zorder=3)
                ax.set_ylim(lo, hi)
                ax.set_xlim(-0.55, len(callers) - 0.45)
                ax.set_xticks(range(len(callers)))
                ax.set_xticklabels([st.CALLER_LABEL[c] for c in callers], rotation=30, ha="right",
                                   rotation_mode="anchor", fontsize=st.FONT_SIZE - 1)
                ax.set_title(label.replace("\n", " ").replace("Imperfect homopol", "Imp. homopol"),
                             fontsize=st.FONT_SIZE)
                y_grid(ax)
                if n % ncols == 0:
                    ax.set_ylabel("Value")
                    if lo > 0:
                        axis_break(ax)
            for n in range(len(strata), nrows * ncols):
                axes[n // ncols][n % ncols].axis("off")
            # callers by colour, the metric by shade (grey key)
            handles = [Patch(color=st.CALLER_COLOR[c], label=st.CALLER_LABEL[c]) for c in callers]
            handles += [Patch(color=tint(st.MUTED, f), label=ml) for _, ml, f in R17_METRICS]
            if other:
                handles = [Patch(color=st.MUTED, label="%s (bars)" % st.REFERENCE_LABEL[ref]),
                           Line2D([], [], color=st.INK, linewidth=1.3, label="%s (line)" %
                                  st.REFERENCE_LABEL[other])] + handles
            fig.legend(handles=handles, loc="outside lower center", ncol=len(handles), frameon=False)
            title = st.REFERENCE_LABEL[ref] if not other else "%s against %s" % (
                st.REFERENCE_LABEL[ref], st.REFERENCE_LABEL[other])
            fig.suptitle("Het SNP calls across genomic regions: %s, %s, %dx" % (
                st.PLATFORM_LABEL[platform], title, ctx.depth))
            tag = ref if not other else "%s_vs_%s" % (ref, other)
            made.append(save(fig, ctx, "figR17_callers_prf_%s_%dx_%s" % (platform, ctx.depth, tag)))
    if not made:
        raise NoData("nothing to draw")
    return ", ".join(made)

FIGURES = [("R5", fig_r5), ("R6", fig_r6), ("R7", fig_r7), ("R8", fig_r8), ("R9", fig_r9), ("R10", fig_r10),
           ("R11", fig_r11), ("R12", fig_r12), ("R13", fig_r13),
           ("R14", fig_r14), ("R15", fig_r15), ("R16", fig_r16),
           ("R17", fig_r17)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=st.ROOT)
    ap.add_argument("--sample", default=st.SAMPLE)
    ap.add_argument("--out", default=None, help="default <root>/plots/pdf")
    ap.add_argument("--only", default="", help="R5,R6,R7,R8,R9")
    ap.add_argument("--caller", default="truth", help="phaser input")
    ap.add_argument("--depth", type=int, default=30)
    ap.add_argument("--phaser", default="whatshap", help="phaser of R5")
    ap.add_argument("--phasers", default="whatshap longphase margin hapcut2 hiphase", help="R6, R7")
    ap.add_argument("--tagged-unit", choices=["reads", "bases"], default="reads", help="R6")
    ap.add_argument("--stratum", default="alldifficultregions", help="end point of R7")
    ap.add_argument("--ref", default="T2T", help="reference of R8")
    ap.add_argument("--metric", choices=["f1", "precision", "recall"], default="f1", help="R8")
    ap.add_argument("--min-label", type=float, default=0.0, help="R5: smallest slice with a label")
    ap.add_argument("--min-slice", type=float, default=0.025,
                    help="R5: smallest drawn slice, as a fraction of the ring (0: true sizes)")
    ap.add_argument("--top", type=int, default=7, help="R9: columns per group")
    ap.add_argument("--png", action="store_true")
    ctx = ap.parse_args()
    ctx.out = ctx.out or os.path.join(ctx.root, "plots", "pdf")
    ctx.phasers = ctx.phasers.split()
    only = {s.strip() for s in ctx.only.split(",") if s.strip()}
    font = st.apply_style(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    if font != "Arial":
        print("note: Arial not found, using %s" % font, file=sys.stderr)
    for name, draw in FIGURES:
        if only and name not in only:
            continue
        try:
            print("%s: %s" % (name, draw(ctx)))
        except NoData as e:
            print("%s: skipped, %s" % (name, e))


if __name__ == "__main__":
    main()
