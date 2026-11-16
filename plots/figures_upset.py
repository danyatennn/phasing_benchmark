#!/usr/bin/env python3
"""Upset plots: a bar per combination of sets, the dots under it tell which
sets it is; groups of bars under a bracket with their total.

  figU1_crossref_<stratum>_<reads>   the reads of a stratum on one reference,
      on the other: aligned (primary, MAPQ >= 20, same chromosome), in its
      confident regions, in the same stratum; both directions; the bar is the
      mean share of the reads over the platforms, the line its range
      (metrics_reads/strata_crossref_sets_30x.tsv, bin/strata_crossref.sh)
  figU2_read_labels          the reads each of the three parental labelling
      methods labels: alignment to the parental assemblies, parental k-mers,
      informative Q100 variants; the mean share of the reads any method labels
      over the platforms, the line its range; the hatched top of a bar, the
      reads whose labels disagree
      (reads/<platform>_T2T_HG002_30x_wgs.consensus.tsv.gz, bin/consensus_truth.sh)

  python plots/figures_upset.py --out plots/pdf_upset
  python plots/figures_upset.py --only U1 --stratum segdups --reads labelled
"""

import argparse
import csv
import gzip
import itertools
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import style as st
from figures import NoData, paper_title, save, stratum_label

DIRECTIONS = [("T2T_to_GRCh38", "T2T-CHM13 → GRCh38"), ("GRCh38_to_T2T", "GRCh38 → T2T-CHM13")]
CROSS_SETS = ["aligned", "confident", "same"]
CROSS_LABEL = {"aligned": "Aligned (MAPQ ≥ 20)", "confident": "In confident regions", "same": "In the same stratum"}
CROSS_SET_COLOR = {"aligned": "#3E8DB3", "confident": "#D9B36A", "same": "#B05AA8"}
LABEL_SETS = ["matpat", "kmer", "q100snp"]
LABEL_NAME = {"matpat": "Parental assembly", "kmer": "Parental k-mers", "q100snp": "Informative Q100 variants"}
LABEL_SET_COLOR = {"matpat": "#3E8DB3", "kmer": "#D9B36A", "q100snp": "#B05AA8"}
DARK, LIGHT = "#7A7A7A", "#7FB6E0"          # the two groups of bars
DOT_ON, DOT_OFF, BAND = "#4D4D4D", "#DDDDDD", "#F1F1F1"


def pct(v):
    return "%.1f" % v if v >= 0.1 else "%.2f" % v


def draw_upset(fig, bar_cell, dot_cell, combos, values, ranges, colors, sets, set_label, set_color, groups,
               ylabel, show_sets=True, hatched=None, title=None, sharey=None):
    """bars in bar_cell over a dot matrix in dot_cell (cells of one grid, so
    that panels side by side line up). values/ranges/colors/hatched by combo;
    groups: [(first index, last index, text)] under brackets."""
    bar_ax, dot_ax = fig.add_subplot(bar_cell, sharey=sharey), fig.add_subplot(dot_cell)
    top = max(r[1] if r else v for v, r in zip(values, ranges))
    for x, (v, r, c) in enumerate(zip(values, ranges, colors)):
        part = hatched[x] if hatched else 0.0
        bar_ax.bar(x, v - part, width=0.7, color=c, zorder=3)
        if part > 0:
            bar_ax.bar(x, part, bottom=v - part, width=0.7, color=c, alpha=0.45, hatch="////",
                       edgecolor="white", linewidth=0, zorder=3)
        if r and r[1] > r[0]:
            bar_ax.plot([x, x], r, color=st.INK, linewidth=0.7, zorder=4)
        bar_ax.annotate(pct(v), (x, r[1] if r else v), xytext=(0, 1.5), textcoords="offset points",
                        ha="center", va="bottom", fontsize=st.SMALL_FONT, fontweight="bold", zorder=5)
    # brackets over the groups, their total on top
    y = top * 1.14
    for first, last, text in groups:
        bar_ax.plot([first - 0.3, first - 0.3, last + 0.3, last + 0.3], [y * 0.97, y, y, y * 0.97],
                    color=st.MUTED, linewidth=0.6, clip_on=False)
        bar_ax.text((first + last) / 2, y * 1.02, text, ha="center", va="bottom", fontsize=st.SMALL_FONT)
    if sharey is None:
        bar_ax.set_ylim(0, top * 1.3)
    bar_ax.set_xlim(-0.6, len(combos) - 0.4)
    bar_ax.set_xticks([])
    bar_ax.spines["bottom"].set_visible(False)
    bar_ax.grid(axis="y", color=st.GRID, linewidth=0.5, zorder=0)
    bar_ax.set_axisbelow(True)
    bar_ax.set_ylabel(ylabel if show_sets else "")
    if not show_sets:
        bar_ax.tick_params(labelleft=False)
    if title:
        bar_ax.set_title(title)
    # the dot matrix: a band under every other set
    for i in range(len(sets)):
        if i % 2 == 0:
            dot_ax.axhspan(i - 0.5, i + 0.5, color=BAND, zorder=0, linewidth=0)
    for x, combo in enumerate(combos):
        ys = [i for i, s in enumerate(sets) if s in combo]
        if len(ys) > 1:
            dot_ax.plot([x, x], [min(ys), max(ys)], color=DOT_ON, linewidth=1.0, zorder=2)
        for i, s in enumerate(sets):
            on = i in ys
            color = (set_color[s] if len(ys) == 1 else DOT_ON) if on else DOT_OFF
            dot_ax.plot(x, i, "o", markersize=3.2, color=color, markeredgewidth=0, zorder=3)
    dot_ax.set_xlim(-0.6, len(combos) - 0.4)
    dot_ax.set_ylim(len(sets) - 0.5, -0.5)
    dot_ax.set_yticks(range(len(sets)))
    dot_ax.set_yticklabels([set_label[s] for s in sets] if show_sets else [])
    dot_ax.set_xticks([])
    for side in ("top", "right", "bottom", "left"):
        dot_ax.spines[side].set_visible(False)
    dot_ax.tick_params(axis="y", length=0)
    return bar_ax


# ------------------------------------------------------------ U1: crossref

def fig_crossref(ctx):
    path = os.path.join(ctx.root, "metrics_reads", "strata_crossref_sets_30x.tsv")
    if not os.path.exists(path):
        raise NoData("no %s - run bin/strata_crossref.sh" % path)
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["stratum"] == ctx.stratum and r["reads"] == ctx.reads]
    if not rows:
        raise NoData("no rows of %s (%s reads) in %s" % (ctx.stratum, ctx.reads, path))
    totals = Counter()
    for r in rows:
        totals[(r["direction"], r["platform"])] += int(r["count"])
    share = defaultdict(lambda: defaultdict(dict))       # direction -> combo -> platform -> %
    for r in rows:
        combo = frozenset(s for s in CROSS_SETS if r[s] == "1")
        share[r["direction"]][combo][r["platform"]] = 100.0 * int(r["count"]) / totals[(r["direction"], r["platform"])]
    allc = [frozenset(c) for k in range(len(CROSS_SETS), -1, -1) for c in itertools.combinations(CROSS_SETS, k)]
    first = share[DIRECTIONS[0][0]]
    mean = lambda d, c: sum(d[c].values()) / len(d[c]) if d.get(c) else 0.0
    # the reads that stay in the stratum first, then the others; largest first in each
    stay = sorted([c for c in allc if "same" in c], key=lambda c: -mean(first, c))
    away = sorted([c for c in allc if "same" not in c], key=lambda c: -mean(first, c))
    combos = stay + away
    w, h = st.PAPER_PANELS["upset_crossref"]
    fig = plt.figure(figsize=(w * st.MM, h * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(2, len(DIRECTIONS), height_ratios=[2.2, 0.9])   # no legend: the brackets name the colours
    first_ax = None
    for i, (direction, title) in enumerate(DIRECTIONS):
        d = share.get(direction, {})
        values = [mean(d, c) for c in combos]
        ranges = [(min(d[c].values()), max(d[c].values())) if d.get(c) else None for c in combos]
        colors = [LIGHT if "same" in c else DARK for c in combos]
        groups = [(0, len(stay) - 1, "Same stratum: %s%%" % pct(sum(values[:len(stay)]))),
                  (len(stay), len(combos) - 1, "Elsewhere: %s%%" % pct(sum(values[len(stay):])))]
        ax = draw_upset(fig, grid[0, i], grid[1, i], combos, values, ranges, colors, CROSS_SETS, CROSS_LABEL, CROSS_SET_COLOR,
                        groups, "Reads of the stratum (%)", show_sets=(i == 0), title=title, sharey=first_ax)
        first_ax = first_ax or ax
    paper_title(fig, "Where the reads of %s land on the other reference, 30×"
                % stratum_label(ctx.stratum).replace("\n", " ").lower())
    return save(fig, ctx, "figU1_crossref_%s_%s" % (ctx.stratum, ctx.reads))


# ------------------------------------------------------ U2: read labels

def label_counts(path):
    """combination of the methods that label a read -> [reads, reads whose labels disagree]"""
    out = defaultdict(lambda: [0, 0])
    with gzip.open(path, "rt") as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        idx = [head.index(s) for s in LABEL_SETS]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            got = [f[i] for i in idx]
            combo = frozenset(s for s, v in zip(LABEL_SETS, got) if v in ("PAT", "MAT"))
            if not combo:
                continue
            out[combo][0] += 1
            out[combo][1] += int(len({v for v in got if v in ("PAT", "MAT")}) > 1)
    return out


def fig_labels(ctx):
    data = {}
    for p in st.PLATFORMS:
        path = os.path.join(ctx.root, "reads", "%s_T2T_HG002_30x_wgs.consensus.tsv.gz" % p)
        if os.path.exists(path):
            data[p] = label_counts(path)
            print("  read labels %s: %d reads" % (p, sum(v[0] for v in data[p].values())), file=sys.stderr)
    if not data:
        raise NoData("no reads/*_T2T_HG002_30x_wgs.consensus.tsv.gz - run bin/consensus_truth.sh")
    total = {p: sum(v[0] for v in d.values()) for p, d in data.items()}
    shares = lambda c, k: [100.0 * d[c][k] / total[p] for p, d in data.items()]
    mean = lambda c, k=0: sum(shares(c, k)) / len(data)
    by_n = [sorted([frozenset(c) for c in itertools.combinations(LABEL_SETS, k)], key=lambda c: -mean(c))
            for k in (3, 2, 1)]
    combos = [c for g in by_n for c in g]
    values = [mean(c) for c in combos]
    ranges = [(min(shares(c, 0)), max(shares(c, 0))) for c in combos]
    hatched = [mean(c, 1) for c in combos]
    palette = {3: "#5A5A5A", 2: "#9A9A9A", 1: LIGHT}
    colors = [palette[len(c)] for c in combos]
    groups, start = [], 0
    for g, name in zip(by_n, ("All three", "Two methods", "One method")):
        end = start + len(g) - 1
        groups.append((start, end, "%s%s%s%%" % (name, ": " if end > start else "\n", pct(sum(values[start:end + 1])))))
        start = end + 1
    w, h = st.PAPER_PANELS["upset_labels"]
    fig = plt.figure(figsize=(w * st.MM, h * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(3, 1, height_ratios=[2.2, 0.9, 0.2])
    draw_upset(fig, grid[0, 0], grid[1, 0], combos, values, ranges, colors, LABEL_SETS, LABEL_NAME, LABEL_SET_COLOR,
               groups, "Labelled reads (%)", hatched=hatched)
    legend_ax = fig.add_subplot(grid[2, 0])
    legend_ax.axis("off")
    legend_ax.legend(handles=[Patch(facecolor="#9A9A9A", alpha=0.45, hatch="////", edgecolor="white",
                                    label="Labels disagree")],
                     loc="center", ncol=1, borderaxespad=0, frameon=False)
    paper_title(fig, "Reads labelled by each parental method, T2T-CHM13, 30×")
    return save(fig, ctx, "figU2_read_labels")


FIGURES = [("U1", fig_crossref), ("U2", fig_labels)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=st.ROOT)
    ap.add_argument("--out", default=None, help="default <root>/plots/pdf")
    ap.add_argument("--only", default="", help="U1, U2, comma-separated")
    ap.add_argument("--stratum", default="alldifficultregions", help="U1: the stratum")
    ap.add_argument("--reads", choices=["labelled", "min2snps"], default="min2snps",
                    help="U1: every labelled read, or those of the recall denominator (>= 2 truth SNPs)")
    ap.add_argument("--png", action="store_true")
    ctx = ap.parse_args()
    ctx.out = ctx.out or os.path.join(ctx.root, "plots", "pdf")
    st.apply_style(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    only = {s.strip() for s in ctx.only.split(",") if s.strip()}
    for name, draw in FIGURES:
        if only and name not in only:
            continue
        try:
            print("%s: %s" % (name, draw(ctx)))
        except NoData as e:
            print("%s: skipped, %s" % (name, e))
        finally:
            plt.close("all")


if __name__ == "__main__":
    main()
