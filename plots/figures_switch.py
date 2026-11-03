#!/usr/bin/env python3
"""Where the phase switches along one chromosome, ONT against PacBio or GRCh38
against T2T-CHM13 (each reference in its own coordinates: no liftover):

  figR4a_switch_overview_<chrom>_<rows>   per row (a platform, or a platform on a
      reference): the switch rate per 1 Mb window (bars), the phase blocks
      (alternating shades), every switch (ticks), and the regions: inside the
      benchmark, segmental duplications, low mappability
  figR4b_switch_depth_<chrom>_<rows>      per row: mean depth of the tagged reads on
      the truth's haplotype (up) and on the other one (down), each block
      oriented against the truth, so a block boundary does not flip the colour
      and a switch inside a block does; the blocks as segments at 0 and the
      switches as ticks on top. --region zooms in.

From metrics_switchtracks/ (bin/switch_tracks.sh), one phaser (--phaser) and one
input (--caller, truth by default: the phaser is the only variable), 30x.

  python plots/figures_switch.py --chrom chr1 --png                    ONT vs PacBio on GRCh38
  python plots/figures_switch.py --chrom chr6 --region chr6:28000000-34000000
  python plots/figures_switch.py --compare references --chrom chr1 --platforms "ont5k revio"
  python plots/figures_switch.py --compare references --chrom chr6 \\
      --region chr6:28000000-34000000 --region-t2t chr6:28400000-34400000
"""

import argparse
import gzip
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

import style as st
from data import sample_root
from figures import NoData, same_decimals, save, top_legend

CONCORDANT = "#B5523B"
DISCORDANT = "#3D6FA3"
BLOCK_SHADES = ("#5B5B5B", "#A8A8A8")
BENCHMARK = "#E4E4E4"
SEGDUPS = "#4A5C70"
LOWMAP = "#9AA9B8"
MB = 1e6


def run_name(ctx, platform, ref):
    return "%s_%s_%s_%dx_wgs_%s_%s" % (platform, ref, ctx.sample, ctx.depth, ctx.caller, ctx.phaser)


def rows_of(ctx):
    """[(platform, reference)] in drawing order"""
    if ctx.compare == "platforms":
        return [(p, ctx.ref) for p in ctx.platforms]
    return [(p, r) for p in ctx.platforms for r in st.REFERENCES]


def row_label(ctx, platform, ref):
    if ctx.compare == "platforms":
        return st.PLATFORM_LABEL[platform]
    return "%s, %s" % (st.PLATFORM_LABEL[platform], st.REFERENCE_LABEL[ref])


def parse_region(text, chrom):
    """chr6:28000000-34000000 -> (chrom, start, end); None for the whole chromosome"""
    if not text:
        return None
    c, _, span = text.partition(":")
    if c != chrom:
        raise NoData("region %s is not on %s" % (text, chrom))
    if not span:
        return None
    a, _, b = span.replace(",", "").partition("-")
    return c, int(a), int(b)


def region_for(ctx, ref):
    return parse_region(ctx.region_t2t if ref == "T2T" else ctx.region, ctx.chrom)


def tsv(path, chrom):
    """rows of a TSV with a header, on one chromosome"""
    if not os.path.exists(path):
        return None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh
                if line.startswith(chrom + "\t")]


def bed(path, chrom):
    """[(start, end)] of a BED on one chromosome"""
    out = []
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.split("\t")
                if f[0] == chrom:
                    out.append((int(f[1]), int(f[2])))
    return out


def chrom_length(ctx, ref):
    path = os.path.join(sample_root(ctx.root, ctx.sample), "data", "chrom_lengths_%s.tsv" % ref)
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.split("\t")
                if f[0] == ctx.chrom:
                    return int(f[1])
    return None


def x_range(ctx, ref):
    region = region_for(ctx, ref)
    if region:
        return region[1], region[2]
    return 0, chrom_length(ctx, ref) or 0


def spans(ax, intervals, y0, h, color, lo, hi):
    for s, e in intervals:
        if e > lo and s < hi:
            ax.add_patch(Rectangle((s / MB, y0), (e - s) / MB, h, facecolor=color,
                                   edgecolor="none", linewidth=0))


# --------------------------------------------------------------- overview

def fig_overview(ctx):
    rows = rows_of(ctx)
    tracks = os.path.join(sample_root(ctx.root, ctx.sample), "metrics_switchtracks")
    data = []
    for platform, ref in rows:
        name = run_name(ctx, platform, ref)
        win = tsv(os.path.join(tracks, name + ".switch_windows.tsv"), ctx.chrom)
        if not win:
            continue
        data.append((platform, ref, win,
                     tsv(os.path.join(tracks, name + ".blocks.tsv"), ctx.chrom) or [],
                     tsv(os.path.join(tracks, name + ".switch_positions.tsv"), ctx.chrom) or []))
    if not data:
        raise NoData("no switch windows for %s in %s - run bin/switch_tracks.sh" % (ctx.chrom, tracks))

    # the rate axis is shared by every row
    top = 0.0
    for _, _, win, _, _ in data:
        for w in win:
            a = int(w["assessed"])
            if a >= ctx.min_assessed:
                top = max(top, 100 * int(w["switches"]) / a)
    top = top * 1.1 or 1.0

    n = len(data)
    fig = plt.figure(figsize=(st.FULL_WIDTH * st.MM, (14 + 30 * n) * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(1 + 2 * n, 1, height_ratios=[0.25] + [1.0, 0.55] * n)
    legend_ax = fig.add_subplot(grid[0])
    legend_ax.axis("off")
    first = None
    for i, (platform, ref, win, blocks, positions) in enumerate(data):
        lo, hi = x_range(ctx, ref)
        rate_ax = fig.add_subplot(grid[1 + 2 * i], sharex=first)
        first = first or rate_ax
        strip = fig.add_subplot(grid[2 + 2 * i], sharex=first)

        for w in win:
            a, s, e = int(w["assessed"]), int(w["start"]), int(w["end"])
            if a < ctx.min_assessed or e <= lo or s >= hi:
                continue
            rate_ax.bar(s / MB, 100 * int(w["switches"]) / a, width=(e - s) / MB, align="edge",
                        color=st.PLATFORM_COLOR[platform], linewidth=0)
        rate_ax.set_ylim(0, top)
        rate_ax.set_ylabel("Switch\nrate (%)")
        rate_ax.set_title(row_label(ctx, platform, ref), loc="left")
        rate_ax.tick_params(labelbottom=False)
        rate_ax.grid(axis="y", color=st.GRID, linewidth=0.5)
        rate_ax.set_axisbelow(True)

        # the strip: regions below, blocks above, switches over the blocks
        prepared = os.path.join(sample_root(ctx.root, ctx.sample), "regions", "prepared_auto")
        easy = bed(os.path.join(prepared, "%s_notinalldifficultregions.bed" % ref), ctx.chrom)
        hard = bed(os.path.join(prepared, "%s_alldifficultregions.bed" % ref), ctx.chrom)
        spans(strip, easy + hard, 0.0, 0.4, BENCHMARK, lo, hi)
        spans(strip, bed(os.path.join(prepared, "%s_lowmappabilityall.bed" % ref), ctx.chrom),
              0.0, 0.4, LOWMAP, lo, hi)
        spans(strip, bed(os.path.join(prepared, "%s_segdups.bed" % ref), ctx.chrom),
              0.0, 0.4, SEGDUPS, lo, hi)
        for k, b in enumerate(sorted(blocks, key=lambda b: int(b["start"]))):
            spans(strip, [(int(b["start"]), int(b["end"]) + 1)], 0.5, 0.3, BLOCK_SHADES[k % 2], lo, hi)
        xs = [(int(p["pos1"]) + int(p["pos2"])) / 2 / MB for p in positions]
        xs = [x for x in xs if lo / MB <= x <= hi / MB]
        if xs:
            strip.vlines(xs, 0.82, 1.0, color=st.INK, linewidth=0.4)
        strip.set_ylim(0, 1)
        strip.set_yticks([])
        strip.spines["left"].set_visible(False)
        if i < n - 1:
            strip.tick_params(labelbottom=False)
        else:
            strip.set_xlabel("%s position (Mb)" % ctx.chrom)
        print("  %s: %d blocks, %d switches on %s" % (row_label(ctx, platform, ref), len(blocks),
                                                      len(positions), ctx.chrom), file=sys.stderr)

    lo, hi = min(x_range(ctx, r)[0] for _, r, _, _, _ in data), max(x_range(ctx, r)[1] for _, r, _, _, _ in data)
    first.set_xlim(lo / MB, hi / MB)
    handles = [Patch(color=BLOCK_SHADES[0], label="Phase blocks"),
               Line2D([], [], color=st.INK, linewidth=0.8, label="Switch error"),
               Patch(color=BENCHMARK, label="Benchmark regions"),
               Patch(color=SEGDUPS, label="Segmental duplications"),
               Patch(color=LOWMAP, label="Low mappability")]
    top_legend(legend_ax, handles)
    return save(fig, ctx, "figR4a_switch_overview_%s_%s" % (ctx.chrom, ctx.tag))


# ------------------------------------------------------------------ depth

def fig_depth(ctx):
    rows = rows_of(ctx)
    tracks = os.path.join(sample_root(ctx.root, ctx.sample), "metrics_switchtracks")
    data = []
    for platform, ref in rows:
        name = run_name(ctx, platform, ref)
        depth = tsv(os.path.join(tracks, name + ".read_depth.tsv.gz"), ctx.chrom)
        if not depth:
            continue
        data.append((platform, ref, depth,
                     tsv(os.path.join(tracks, name + ".blocks.tsv"), ctx.chrom) or [],
                     tsv(os.path.join(tracks, name + ".switch_positions.tsv"), ctx.chrom) or []))
    if not data:
        raise NoData("no read depth for %s in %s - run bin/switch_tracks.sh with the read tables" % (
            ctx.chrom, tracks))

    n = len(data)
    fig, axes = plt.subplots(n, 1, figsize=(st.FULL_WIDTH * st.MM, (12 + 32 * n) * st.MM),
                             constrained_layout=True, squeeze=False)
    axes = axes[:, 0]
    ymax = 0.0
    for ax, (platform, ref, depth, blocks, positions) in zip(axes, data):
        lo, hi = x_range(ctx, ref)
        d = [r for r in depth if int(r["end"]) > lo and int(r["start"]) < hi]
        x = np.array([int(r["start"]) for r in d] + [int(d[-1]["end"])]) / MB if d else np.array([])
        conc = np.array([float(r["concordant"]) for r in d])
        disc = np.array([float(r["discordant"]) for r in d])
        if len(d):
            ax.stairs(conc, x, fill=True, color=CONCORDANT, linewidth=0)
            ax.stairs(-disc, x, fill=True, color=DISCORDANT, linewidth=0)
            ymax = max(ymax, conc.max(), disc.max())
        # blocks as segments at 0, alternating shades so every boundary shows
        for k, b in enumerate(sorted(blocks, key=lambda b: int(b["start"]))):
            s, e = int(b["start"]), int(b["end"])
            if e > lo and s < hi:
                ax.plot([s / MB, e / MB], [0, 0], color=BLOCK_SHADES[k % 2], linewidth=2.5,
                        solid_capstyle="butt", zorder=3)
        xs = [(int(p["pos1"]) + int(p["pos2"])) / 2 / MB for p in positions]
        xs = [v for v in xs if lo / MB <= v <= hi / MB]
        if xs:
            ax.plot(xs, [1.0] * len(xs), linestyle="none", marker="v", markersize=3, color=st.INK,
                    transform=ax.get_xaxis_transform(), clip_on=False, zorder=4)
        ax.set_xlim(lo / MB, hi / MB)
        ax.axhline(0, color=st.INK, linewidth=0.4)
        ax.text(1.005, 0.5, row_label(ctx, platform, ref), transform=ax.transAxes, rotation=270,
                ha="left", va="center")
        print("  %s: %d switches in view" % (row_label(ctx, platform, ref), len(xs)), file=sys.stderr)
    lim = ctx.depth_ymax or math.ceil(ymax / 10) * 10 or 10
    for ax in axes:
        ax.set_ylim(-lim, lim)
        ax.set_yticks([-lim, -lim / 2, 0, lim / 2, lim])
        ax.set_yticklabels(same_decimals([abs(v) for v in ax.get_yticks()]))
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    axes[-1].set_xlabel("%s position (Mb)" % ctx.chrom)
    fig.supylabel("Read depth: truth haplotype (up) / other haplotype (down)", fontsize=st.FONT_SIZE)
    return save(fig, ctx, "figR4b_switch_depth_%s_%s" % (ctx.chrom, ctx.tag))


FIGURES = [("overview", fig_overview), ("depth", fig_depth)]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--sample", default=st.SAMPLE)
    parser.add_argument("--out", default=None, help="default <root>/plots/pdf")
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--compare", choices=["platforms", "references"], default="platforms")
    parser.add_argument("--ref", default="GRCh38", help="reference of --compare platforms")
    parser.add_argument("--platforms", default=" ".join(st.PLATFORMS))
    parser.add_argument("--phaser", default="whatshap")
    parser.add_argument("--caller", default="truth")
    parser.add_argument("--depth", type=int, default=30)
    parser.add_argument("--region", default="", help="chr:start-end to zoom in (GRCh38, or --ref)")
    parser.add_argument("--region-t2t", default="", help="the same region in T2T coordinates")
    parser.add_argument("--min-assessed", type=int, default=20,
                        help="fewest assessed SNP pairs behind a window's switch rate")
    parser.add_argument("--depth-ymax", type=float, default=None, help="y limit of figR4b")
    parser.add_argument("--only", default="", help="overview, depth, comma-separated")
    parser.add_argument("--png", action="store_true", help="a PNG next to each PDF")
    ctx = parser.parse_args()
    ctx.out = ctx.out or os.path.join(ctx.root, "plots", "pdf")
    ctx.platforms = ctx.platforms.split()
    if ctx.compare == "references" and ctx.region and not ctx.region_t2t:
        sys.exit("--compare references with --region needs --region-t2t: the references have "
                 "different coordinates")
    if ctx.compare == "platforms" and ctx.ref == "T2T":
        ctx.region_t2t = ctx.region_t2t or ctx.region
    ctx.tag = (ctx.ref if ctx.compare == "platforms" else "refs") + "_" + ctx.phaser
    if ctx.region:
        ctx.tag += "_" + ctx.region.split(":")[-1].replace(",", "")
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
