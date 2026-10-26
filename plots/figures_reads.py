#!/usr/bin/env python3
"""Het SNPs hit by one read against its haplotag, one PDF per reference each:

  figR3a_accuracy_vs_snps_<ref>  haplotag accuracy of the tagged reads against
                                 the number of het SNPs the read hits: a panel
                                 per phaser, a line per platform
  figR3b_tagged_vs_snps_<ref>    share of the reads tagged, against the same
  figR3c_snps_per_read_<ref>     how many SNPs the reads hit: reads sorted by it,
                                 the count against the cumulative share of reads

From metrics_reads/summary_snps_accuracy.tsv (bin/snps_accuracy.sh), 30x,
chr1-22. "SNPs hit" (--count):
  snps   het SNPs of the input the read covers with REF or ALT
  votes  phased SNPs of the read's block that carry an allele of it - what
         whatshap haplotag decided on
Accuracy: tagged reads (H1/H2 in a block of two or more SNPs) with a truth
label, their block oriented against the truth by its SNPs; correct / compared.
A point with fewer than --min-compared reads is left out.

  python plots/figures_reads.py --png
  python plots/figures_reads.py --count votes --phasers whatshap --refs T2T
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

import style as st
from data import sample_root
from figures import NoData, panel_row, platform_handles, save, top_legend, y_grid

# groups of the SNP count, wider as it grows: (low, high, label)
GROUPS = [(1, 1, "1"), (2, 2, "2"), (3, 3, "3"), (4, 4, "4"), (5, 6, "5-6"), (7, 9, "7-9"),
          (10, 14, "10-14"), (15, 19, "15-19"), (20, 29, "20-29"), (30, 49, "30-49"),
          (50, 99, "50-99"), (100, 199, "100-199"), (200, 10 ** 9, "200+")]
COUNT_LABEL = {"snps": "Het SNPs hit by the read", "votes": "Phased SNPs the read votes with"}


def load(root, sample, ref, depth, caller, count):
    """(platform, phaser) -> {n: [reads, tagged, compared, correct]}"""
    path = os.path.join(sample_root(root, sample), "metrics_reads", "summary_snps_accuracy.tsv")
    if not os.path.exists(path):
        raise NoData("no %s - run bin/snps_accuracy.sh" % path)
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if (r.get("reference") != ref or r.get("caller") != caller
                    or r.get("depth") != str(depth) or r.get("count") != count):
                continue
            x = out.setdefault((r["platform"], r["phaser"]), {})
            x[int(r["n"])] = [int(r["reads"]), int(r["tagged"]), int(r["compared"]), int(r["correct"])]
    if not out:
        raise NoData("no %s %s %dx rows for %s in %s" % (caller, count, depth, ref, path))
    return out


def grouped(table, i, j, min_den):
    """x (geometric middle of each group), 100 * sum(col i) / sum(col j) per
    group; NaN where the denominator is below min_den."""
    xs, ys = [], []
    for lo, hi, _ in GROUPS:
        num = sum(v[i] for n, v in table.items() if lo <= n <= hi)
        den = sum(v[j] for n, v in table.items() if lo <= n <= hi)
        xs.append(math.sqrt(lo * min(hi, 2 * lo if hi > 10 ** 6 else hi)))
        ys.append(100 * num / den if den >= min_den else math.nan)
    return xs, ys


def group_axis(ax, last):
    """Log x axis labelled with the groups, up to the last one with a point."""
    ax.set_xscale("log")
    ticks = [math.sqrt(lo * min(hi, 2 * lo if hi > 10 ** 6 else hi)) for lo, hi, _ in GROUPS][:last + 1]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, i: GROUPS[i][2] if i is not None and i < len(ticks) else ""))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis="x", labelrotation=90)
    ax.set_xlim(ticks[0] / 1.3, ticks[-1] * 1.3)


def by_phaser(ctx, ref, name, i, j, min_den, ylabel, ylim=None):
    data = load(ctx.root, ctx.sample, ref, ctx.depth, ctx.caller, ctx.count)
    phasers = [p for p in ctx.phasers if any(k[1] == p for k in data)]
    if not phasers:
        raise NoData("none of %s in the table" % " ".join(ctx.phasers))
    fig, legend_ax, axes = panel_row(len(phasers), height=st.PANEL_HEIGHT + 12, sharex=True)
    drawn, last = [], 0
    for ax, phaser in zip(axes, phasers):
        for platform in st.PLATFORMS:
            table = data.get((platform, phaser))
            if not table:
                continue
            xs, ys = grouped(table, i, j, min_den)
            if not any(math.isfinite(y) for y in ys):
                continue
            ax.plot(xs, ys, color=st.PLATFORM_COLOR[platform], marker=st.PLATFORM_MARKER[platform],
                    markeredgecolor="white", zorder=3)
            drawn += [y for y in ys if math.isfinite(y)]
            last = max([last] + [g for g, y in enumerate(ys) if math.isfinite(y)])
        ax.set_title(st.PHASER_LABEL[phaser])
        ax.set_xlabel(COUNT_LABEL[ctx.count])
        y_grid(ax)
    if not drawn:
        plt.close(fig)
        raise NoData("no group with %d reads or more" % min_den)
    for ax in axes:
        group_axis(ax, last)
    if ylim:
        axes[0].set_ylim(*ylim)
    else:
        lo = min(drawn)
        axes[0].set_ylim(max(0, lo - 0.1 * (100 - lo)), 100.5)
    axes[0].set_ylabel(ylabel)
    present = {st.PLATFORM_LABEL[k[0]] for k in data}
    top_legend(legend_ax, [h for h in platform_handles() if h.get_label() in present])
    return save(fig, ctx, name)


def fig_accuracy(ctx, ref):
    return by_phaser(ctx, ref, "figR3a_accuracy_vs_snps_%s" % ref, 3, 2, ctx.min_compared,
                     "Haplotag accuracy (%)")


def fig_tagged(ctx, ref):
    return by_phaser(ctx, ref, "figR3b_tagged_vs_snps_%s" % ref, 1, 0, ctx.min_compared,
                     "Reads tagged (%)", ylim=(0, 100))


def fig_distribution(ctx, ref):
    """Reads sorted by the SNPs they hit, most first: the count against the
    share of reads with at least that many - the curve of a block-size plot.
    The SNP count does not depend on the phaser, so any one run will do."""
    data = load(ctx.root, ctx.sample, ref, ctx.depth, ctx.caller, "snps")
    fig, ax = plt.subplots(figsize=(90 * st.MM, 70 * st.MM), constrained_layout=True)
    drawn, top = [], 0
    for platform in st.PLATFORMS:
        tables = [v for k, v in data.items() if k[0] == platform]
        if not tables:
            continue
        table = tables[0]
        ns = sorted(table, reverse=True)
        reads = np.array([table[n][0] for n in ns], dtype=float)
        total = reads.sum()
        if not total:
            continue
        # a step per value: from the share before it to the share after it
        cum = np.concatenate([[0.0], np.cumsum(reads) / total])
        xs = np.repeat(cum, 2)[1:-1]
        ys = np.repeat(np.array(ns, dtype=float), 2)
        ax.plot(xs, ys, color=st.PLATFORM_COLOR[platform], label=st.PLATFORM_LABEL[platform])
        drawn.append(platform)
        median = ns[int(np.searchsorted(cum[1:], 0.5))]
        # the top 1% of the reads would squash the rest of the curve
        top = max(top, ns[int(np.searchsorted(cum[1:], 0.01))])
        print("  %s %s: median %d SNPs per read" % (ref, platform, median), file=sys.stderr)
    if not drawn:
        plt.close(fig)
        raise NoData("no read counts for %s" % ref)
    ax.axvline(0.5, color=st.INK, linestyle="-.", linewidth=0.8, zorder=1)
    ax.set_xlim(0, 1)
    # the curve starts at the top 1% of the reads; the rare longest ones run
    # off the axis (--ymax to show them)
    ax.set_ylim(0, ctx.ymax if ctx.ymax else 1.15 * top)
    ax.set_xlabel("Cumulative share of reads")
    ax.set_ylabel(COUNT_LABEL["snps"])
    ax.set_title("Het SNPs per read, %s, %dx" % (st.REFERENCE_LABEL[ref], ctx.depth))
    ax.legend(loc="upper right")
    y_grid(ax)
    return save(fig, ctx, "figR3c_snps_per_read_%s" % ref)


FIGURES = [("accuracy", fig_accuracy), ("tagged", fig_tagged), ("distribution", fig_distribution)]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--sample", default=st.SAMPLE)
    parser.add_argument("--out", default=None, help="default <root>/plots/pdf")
    parser.add_argument("--refs", default=" ".join(st.REFERENCES))
    parser.add_argument("--depth", type=int, default=30)
    parser.add_argument("--only", default="", help="accuracy, tagged, distribution, comma-separated")
    parser.add_argument("--caller", default="truth", help="phaser input")
    parser.add_argument("--count", choices=["snps", "votes"], default="snps")
    parser.add_argument("--phasers", default="whatshap longphase margin hapcut2 hiphase")
    parser.add_argument("--min-compared", type=int, default=200,
                        help="fewest reads behind a point")
    parser.add_argument("--ymax", type=float, default=None, help="top of the y axis of figR3c")
    parser.add_argument("--png", action="store_true", help="a PNG next to each PDF")
    ctx = parser.parse_args()
    ctx.out = ctx.out or os.path.join(ctx.root, "plots", "pdf")
    ctx.phasers = ctx.phasers.split()
    only = {s.strip() for s in ctx.only.split(",") if s.strip()}

    font = st.apply_style(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    if font != "Arial":
        print("note: Arial not found, using %s" % font, file=sys.stderr)

    for ref in ctx.refs.split():
        for name, draw in FIGURES:
            if only and name not in only:
                continue
            try:
                print("%s %s: %s" % (name, ref, draw(ctx, ref)))
            except NoData as e:
                print("%s %s: skipped, %s" % (name, ref, e))


if __name__ == "__main__":
    main()
