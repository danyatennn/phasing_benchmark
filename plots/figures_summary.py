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
from matplotlib.patches import Patch

import style as st
from data import sample_root
from figures import NoData, log_ticks, save, series_handle, stratum_label, top_legend, y_grid


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


def donut(ax, counts, center, min_label, min_slice):
    total = sum(counts.values())
    vals = [(k, col, hatch, counts.get(k, 0)) for k, _, col, hatch in CLASSES if counts.get(k, 0) > 0]
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
            ax.annotate(text, xy=(1.01 * x, 1.01 * y), xytext=(sign * 1.3, ty),
                        ha="left" if sign > 0 else "right", va="center", color=st.INK,
                        fontsize=st.FONT_SIZE - 1,
                        arrowprops=dict(arrowstyle="-", color=st.MUTED, linewidth=0.5, shrinkA=0, shrinkB=0))
    ax.set_xlim(-2.45, 2.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")


def fig_r5(ctx):
    rows = [r for r in table(path_of(ctx, "metrics_reads", "summary_classes.tsv"))
            if r.get("caller") == ctx.caller and r.get("depth") == str(ctx.depth)
            and r.get("phaser") == ctx.phaser and r.get("stratum") == "genome"]
    if not rows:
        raise NoData("no read classes for %s %s %dx" % (ctx.caller, ctx.phaser, ctx.depth))
    data = {}
    for r in rows:
        key = (r["reference"], r["platform"])
        if r["group_type"] in ("class", "all"):
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
            untagged = {k: counts.get(k, 0) for k, *_ in CLASSES}
            total_reads = counts.get("all", math.nan)
            n = sum(untagged.values())
            center = "%s\nuntagged\n(%.1f%%)" % (short(n), 100 * n / total_reads) \
                if total_reads == total_reads else "%s\nuntagged" % short(n)
            donut(ax, untagged, center, ctx.min_label, ctx.min_slice)
            if i == 0:
                ax.set_title(st.REFERENCE_LABEL[ref], pad=4)
        fig.text(0.02, top - (top - bottom) * (i + 0.5) / len(platforms), st.PLATFORM_LABEL[platform],
                 rotation=90, ha="center", va="center")
    handles = [Patch(facecolor=c, hatch=h, edgecolor="white", label=l) for _, l, c, h in CLASSES]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.55, 0.02))
    if ctx.min_slice > 0:
        fig.text(0.55, 0.005, "Slices under %g%% are drawn at %g%% so they can be seen; the labels give the "
                 "true numbers." % (100 * ctx.min_slice, 100 * ctx.min_slice), ha="center", va="bottom",
                 color=st.MUTED, fontsize=st.FONT_SIZE - 1)
    fig.suptitle("Why reads stay untagged: %s, %s input, %dx, chr1-22" % (
        st.PHASER_LABEL[ctx.phaser], st.CALLER_LABEL.get(ctx.caller, ctx.caller), ctx.depth), y=0.99)
    return save(fig, ctx, "figR5_untagged_reasons_%s" % ctx.phaser)


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


FIGURES = [("R5", fig_r5), ("R6", fig_r6), ("R7", fig_r7), ("R8", fig_r8), ("R9", fig_r9)]


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
