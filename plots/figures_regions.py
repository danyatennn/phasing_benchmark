#!/usr/bin/env python3
"""Two figures of the genomic regions, one PDF per reference each:

  figR1_vaf_strata_<ref>      VAF of the truth het SNPs as each platform's reads
                              see them: a row per region, a column per platform,
                              the histogram as a density. A clean het sits at
                              0.5; reads of a paralogous copy pull it away.
                              From metrics_strata_auto/vaf_<ref>.tsv
                              (bin/vaf_strata.sh).
  figR2_tagged_strata_<ref>   reads tagged per region: a panel per region, a
                              bar per platform (the mean over the phasers run
                              on every platform), a point per phaser. A read is
                              tagged if whatshap haplotag put it on H1/H2 in a
                              block of two or more SNPs, and is in a region if
                              half its alignment or more lies there. From
                              metrics_reads/<run>.classes.tsv (bin/read_classes.sh).

Both are at one coverage (--depth, 30 by default: the allele matrices and the
read classes exist only there) and over chr1-22.

  python plots/figures_regions.py                       both figures, GRCh38 and T2T
  python plots/figures_regions.py --only vaf --refs T2T --png
  python plots/figures_regions.py --root /copy/of/the/project --out /tmp/pdf
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import style as st
from data import Metrics, sample_root
from figures import NoData, finite, save, series_handle, stratum_label, top_legend, y_grid

PLATFORM_SHORT = {"ont4k": "ONT\n4 kHz", "ont5k": "ONT\n5 kHz", "hifi": "PacBio\nSequel II",
                  "revio": "PacBio\nRevio"}


def region_order(db, present):
    """Whole genome first, then the tier-1 strata in the order of the other
    figures, then anything else; only regions with data."""
    order = ["genome"] + db.strata_names()
    order += sorted(set(present) - set(order))
    return [r for r in order if r in present]


# ------------------------------------------------------------ VAF histograms

def load_vaf(root, sample, ref, depth, source):
    """(platform, stratum) -> counts per bin, in the order of the bins."""
    path = os.path.join(sample_root(root, sample), "metrics_strata_auto", "vaf_%s.tsv" % ref)
    if not os.path.exists(path):
        raise NoData("no %s - run bin/vaf_strata.sh" % path)
    hist = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r.get("depth") != str(depth) or r.get("input") != source:
                continue
            k = (r["platform"], r["stratum"])
            hist.setdefault(k, []).append(int(r["sites"]))
    if not hist:
        raise NoData("no %s %dx rows in %s" % (source, depth, path))
    return hist


def rebin(counts, factor):
    n = len(counts) // factor * factor
    return np.asarray(counts[:n], dtype=float).reshape(-1, factor).sum(axis=1)


def fig_vaf(ctx, db, ref):
    hist = load_vaf(ctx.root, ctx.sample, ref, ctx.depth, ctx.input)
    platforms = [p for p in st.PLATFORMS if any(k[0] == p for k in hist)]
    regions = region_order(db, {k[1] for k in hist if sum(hist[k])})
    if not platforms or not regions:
        raise NoData("no VAF histogram with sites for %s" % ref)
    nbins = len(next(iter(hist.values())))
    factor = max(1, int(round(ctx.bin_width * nbins)))
    width = factor / nbins

    rows, cols = len(regions), len(platforms)
    fig = plt.figure(figsize=(st.FULL_WIDTH * st.MM, (12 + 21 * rows) * st.MM), constrained_layout=True)
    # a block per region, its name over the row: too long for a rotated label
    blocks = fig.subfigures(rows, 1, squeeze=False)[:, 0]
    for i, region in enumerate(regions):
        blocks[i].suptitle(stratum_label(region).replace("\n", " "), x=0.01, ha="left",
                           fontsize=st.FONT_SIZE)
        axes = blocks[i].subplots(1, cols, sharey=True, squeeze=False)[0]
        peak = 0.0
        for j, platform in enumerate(platforms):
            ax = axes[j]
            counts = hist.get((platform, region))
            total = sum(counts) if counts else 0
            if total:
                c = rebin(counts, factor)
                density = c / (c.sum() * width)
                peak = max(peak, density.max())
                lo = np.arange(len(c)) * width
                ax.bar(lo, density, width=width, align="edge", color=st.PLATFORM_COLOR[platform],
                       linewidth=0, zorder=2)
                mids = lo + width / 2
                median = mids[np.searchsorted(np.cumsum(c), c.sum() / 2)]
                ax.text(0.98, 0.95, "n = %s\nmedian %.2f" % (count_label(total), median),
                        transform=ax.transAxes, ha="right", va="top", color=st.MUTED)
            else:
                ax.text(0.5, 0.5, "no sites", transform=ax.transAxes, ha="center", va="center",
                        color=st.MUTED)
            ax.axvline(0.5, color=st.MUTED, linewidth=0.5, linestyle="--", zorder=1)
            ax.set_xlim(0, 1)
            ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
            ax.set_xticklabels(["0", "", "0.5", "", "1"])
            if i == 0:
                ax.set_title(st.PLATFORM_LABEL[platform])
            if j == 0:
                ax.set_ylabel("Density")
            if i == rows - 1:
                ax.set_xlabel("VAF")
            else:
                ax.tick_params(labelbottom=False)
        # room over the bars for the counts
        if peak:
            axes[0].set_ylim(0, 1.35 * peak)
    return save(fig, ctx, "figR1_vaf_strata_%s" % ref)


def count_label(n):
    if n >= 1e6:
        return "%.2fM" % (n / 1e6)
    if n >= 1e3:
        return "%.1fk" % (n / 1e3)
    return "%d" % n


# ---------------------------------------------------- reads tagged per region

def tagged_pct(table, region):
    tagged = table.get((region, "class", "tagged"), {}).get("reads", math.nan)
    total = table.get((region, "all", "all"), {}).get("reads", math.nan)
    return 100 * tagged / total if finite(tagged) and finite(total) and total else math.nan


def fig_tagged(ctx, db, ref):
    # phaser -> platform -> region -> % tagged
    values = {}
    present = set()
    for phaser in ctx.phasers:
        for platform in st.PLATFORMS:
            if not st.applies(platform, phaser):
                continue
            table = db.classes((ref, ctx.caller, platform, ctx.depth, phaser))
            for region in {k[0] for k in table}:
                v = tagged_pct(table, region)
                if finite(v):
                    values.setdefault(phaser, {}).setdefault(platform, {})[region] = v
                    present.add(region)
    if not values:
        raise NoData("no read classes for %s %s %dx" % (ref, ctx.caller, ctx.depth))
    regions = region_order(db, present)
    platforms = [p for p in st.PLATFORMS if any(p in v for v in values.values())]
    phasers = [p for p in ctx.phasers if p in values]
    # the bar is over the phasers run on every platform, so the platforms share it
    common = [p for p in phasers if all(st.applies(pl, p) for pl in platforms)]

    ncols = min(4, len(regions))
    nrows = math.ceil(len(regions) / ncols)
    fig = plt.figure(figsize=(st.FULL_WIDTH * st.MM, (10 + 42 * nrows) * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(nrows + 1, ncols, height_ratios=[0.12] + [1] * nrows)
    legend_ax = fig.add_subplot(grid[0, :])
    legend_ax.axis("off")

    x = np.arange(len(platforms))
    offsets = np.linspace(-0.27, 0.27, len(phasers)) if len(phasers) > 1 else [0.0]
    first = None
    for n, region in enumerate(regions):
        ax = fig.add_subplot(grid[1 + n // ncols, n % ncols], sharey=first)
        first = first or ax
        for i, platform in enumerate(platforms):
            color = st.PLATFORM_COLOR[platform]
            bar = [values[p][platform][region] for p in common
                   if region in values[p].get(platform, {})]
            if bar:
                ax.bar(i, np.mean(bar), width=0.75, color=color, alpha=0.3, linewidth=0, zorder=2)
            for off, phaser in zip(offsets, phasers):
                v = values[phaser].get(platform, {}).get(region)
                if v is not None:
                    ax.plot(i + off, v, marker=st.PHASER_MARKER[phaser], color=color,
                            markeredgecolor="white", linestyle="none", zorder=3)
        ax.set_title(stratum_label(region).replace("\n", " "))
        ax.set_xticks(x)
        ax.set_xticklabels([PLATFORM_SHORT[p] for p in platforms])
        ax.set_xlim(-0.6, len(platforms) - 0.4)
        ax.set_ylim(0, 100)
        y_grid(ax)
        if n % ncols:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("Reads tagged (%)")

    handles = [Patch(color=st.PLATFORM_COLOR[p], alpha=0.3, label=st.PLATFORM_LABEL[p]) for p in platforms]
    handles += [series_handle(st.INK, st.PHASER_MARKER[p], st.PHASER_LABEL[p], line=False) for p in phasers]
    top_legend(legend_ax, handles, ncol=min(len(handles), 5))
    return save(fig, ctx, "figR2_tagged_strata_%s" % ref)


# ------------------------------------------------------------------ main

FIGURES = [("vaf", fig_vaf), ("tagged", fig_tagged)]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--sample", default=st.SAMPLE)
    parser.add_argument("--out", default=None, help="default <root>/plots/pdf")
    parser.add_argument("--refs", default=" ".join(st.REFERENCES))
    parser.add_argument("--depth", type=int, default=30)
    parser.add_argument("--only", default="", help="vaf, tagged or both, comma-separated")
    parser.add_argument("--input", default="truth", help="sites of the VAF figure: truth or a caller")
    parser.add_argument("--caller", default="truth", help="phaser input of the tagged-reads figure")
    parser.add_argument("--phasers", default=" ".join(st.PHASERS))
    parser.add_argument("--bin-width", type=float, default=0.02)
    parser.add_argument("--png", action="store_true", help="a PNG next to each PDF")
    ctx = parser.parse_args()
    ctx.out = ctx.out or os.path.join(ctx.root, "plots", "pdf")
    ctx.phasers = ctx.phasers.split()
    only = {s.strip() for s in ctx.only.split(",") if s.strip()}

    font = st.apply_style(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    if font != "Arial":
        print("note: Arial not found, using %s" % font, file=sys.stderr)

    db = Metrics(ctx.root, ctx.sample, st.REGION)
    for ref in ctx.refs.split():
        for name, draw in FIGURES:
            if only and name not in only:
                continue
            try:
                print("%s %s: %s" % (name, ref, draw(ctx, db, ref)))
            except NoData as e:
                print("%s %s: skipped, %s" % (name, ref, e))


if __name__ == "__main__":
    main()
