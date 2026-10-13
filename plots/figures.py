"""The figures. Every function draws one subfigure into its own PDF and raises
NoData when the metrics it needs are not there yet.

Which denominators are shared, and so what may be compared with what:
  phasers on the GIAB truth input   every phaser gets the same sites, so
                                    switch rate, N50 and phased fraction compare
                                    directly (fig1-fig4, fig7, figS1)
  callers                           each caller hands the phaser its own sites,
                                    so the switch rate is conditional on what was
                                    called (fig5a). The comparable numbers are
                                    the ones over truth_total, which is one
                                    constant per reference (fig5b)
  references                        truth_total differs between GRCh38 and T2T,
                                    so only rates are put side by side (fig6)
"""

import math
import os

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, LogLocator, NullLocator

import style as st
from data import NAN, ratio


class NoData(Exception):
    pass


def finite(v):
    return isinstance(v, (int, float)) and math.isfinite(v)


def positive(v):
    return v if finite(v) and v > 0 else NAN


# ------------------------------------------------------------------ plumbing

def save(fig, ctx, name):
    os.makedirs(ctx.out, exist_ok=True)
    fig.savefig(os.path.join(ctx.out, name + ".pdf"))
    if ctx.png:
        fig.savefig(os.path.join(ctx.out, name + ".png"))
    plt.close(fig)
    return name + ".pdf"


def panel_row(n, height=st.PANEL_HEIGHT, width=st.FULL_WIDTH, sharex=False):
    """A strip for the legend on top, n panels with one shared y axis below."""
    fig = plt.figure(figsize=(width * st.MM, height * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(2, n, height_ratios=[0.08, 1])
    legend_ax = fig.add_subplot(grid[0, :])
    legend_ax.axis("off")
    axes = []
    for i in range(n):
        ax = fig.add_subplot(
            grid[1, i],
            sharey=axes[0] if axes else None,
            sharex=axes[0] if (axes and sharex) else None)
        if i:
            ax.tick_params(labelleft=False)
        axes.append(ax)
    return fig, legend_ax, axes


def log_ticks(axis, lo, hi):
    """Plain numbers on a log axis, no 10^x - it would pull in the math font."""
    wide = finite(lo) and finite(hi) and lo > 0 and hi / lo > 300
    axis.set_major_locator(LogLocator(base=10, subs=(1,) if wide else (1, 2, 5), numticks=15))
    axis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    axis.set_minor_locator(NullLocator())


def y_grid(ax):
    ax.grid(axis="y", color=st.GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)


def empty_panel(ax):
    ax.text(0.5, 0.5, "not computed yet", transform=ax.transAxes,
            ha="center", va="center", color=st.MUTED)


def series_handle(color, marker, label, line=True):
    return Line2D([], [], color=color if line else "none", marker=marker,
                  markerfacecolor=color, markeredgecolor="white", label=label)


def phaser_handles(line=True):
    return [series_handle(st.PHASER_COLOR[p], st.PHASER_MARKER[p], st.PHASER_LABEL[p], line)
            for p in st.PHASERS]


def platform_handles(line=True):
    return [series_handle(st.PLATFORM_COLOR[p], st.PLATFORM_MARKER[p], st.PLATFORM_LABEL[p], line)
            for p in st.PLATFORMS]


def top_legend(legend_ax, handles, **kwargs):
    legend_ax.legend(handles=handles, loc="center", ncol=len(handles),
                     borderaxespad=0, **kwargs)


# ------------------------------------------------- lines over coverage: 1, 2, 7

def coverage_lines(db, ctx, name, ref, value, ylabel, by, log=False, caller="truth"):
    """by='platform': a panel per platform, a line per phaser.
       by='phaser':   a panel per phaser, a line per platform."""
    if by == "platform":
        panels, series = st.PLATFORMS, st.PHASERS
        panel_label, color, marker = st.PLATFORM_LABEL, st.PHASER_COLOR, st.PHASER_MARKER
        handles = phaser_handles()
        key = lambda panel, s, depth: (ref, caller, panel, depth, s)
    else:
        panels, series = st.PHASERS, st.PLATFORMS
        panel_label, color, marker = st.PHASER_LABEL, st.PLATFORM_COLOR, st.PLATFORM_MARKER
        handles = platform_handles()
        key = lambda panel, s, depth: (ref, caller, s, depth, panel)

    fig, legend_ax, axes = panel_row(len(panels))
    drawn = []
    for ax, panel in zip(axes, panels):
        lines = 0
        for s in series:
            ys = [value(db, key(panel, s, d)) for d in st.DEPTHS]
            if log:
                ys = [positive(y) for y in ys]
            if not any(finite(y) for y in ys):
                continue
            # a depth that is not computed breaks the line instead of being bridged
            ax.plot(st.DEPTHS, ys, color=color[s], marker=marker[s],
                    markeredgecolor="white", zorder=3)
            drawn.extend(y for y in ys if finite(y))
            lines += 1
        ax.set_title(panel_label[panel])
        ax.set_xticks(st.DEPTHS)
        ax.set_xlim(st.DEPTHS[0] - 5, st.DEPTHS[-1] + 5)
        ax.set_xlabel("Coverage (×)")
        y_grid(ax)
        if not lines:
            empty_panel(ax)
    if not drawn:
        plt.close(fig)
        raise NoData("no %s rows for %s" % (caller, ref))

    if log:
        axes[0].set_yscale("log")
        log_ticks(axes[0].yaxis, min(drawn), max(drawn))
    axes[0].set_ylabel(ylabel)
    top_legend(legend_ax, handles)
    return save(fig, ctx, name)


def switch_pct(db, key):
    return 100 * db.get(key, "switch_rate")


def n50_mb(db, key):
    return db.get(key, "block_n50") / 1e6


def phased_pct(db, key):
    return 100 * db.get(key, "phased_fraction")


def tagged_pct(db, key):
    return 100 * db.get(key, "tagging_completeness")


def haplotag_pct(db, key):
    return 100 * db.get(key, "haplotag_accuracy")


def fig1a(db, ctx, ref):
    return coverage_lines(db, ctx, "fig1a_switch_rate_by_platform_" + ref, ref,
                          switch_pct, "Switch error rate (%)", "platform", log=True)


def fig1b(db, ctx, ref):
    return coverage_lines(db, ctx, "fig1b_block_n50_by_platform_" + ref, ref,
                          n50_mb, "Phase block N50 (Mb)", "platform", log=True)


def fig1c(db, ctx, ref):
    return coverage_lines(db, ctx, "fig1c_phased_fraction_by_platform_" + ref, ref,
                          phased_pct, "Phased het SNPs (%)", "platform")


def fig2a(db, ctx, ref):
    return coverage_lines(db, ctx, "fig2a_switch_rate_by_phaser_" + ref, ref,
                          switch_pct, "Switch error rate (%)", "phaser", log=True)


def fig2b(db, ctx, ref):
    return coverage_lines(db, ctx, "fig2b_block_n50_by_phaser_" + ref, ref,
                          n50_mb, "Phase block N50 (Mb)", "phaser", log=True)


def fig7a(db, ctx, ref):
    return coverage_lines(db, ctx, "fig7a_reads_tagged_" + ref, ref,
                          tagged_pct, "Reads assigned a haplotype (%)", "platform")


def fig7b(db, ctx, ref):
    return coverage_lines(db, ctx, "fig7b_haplotag_accuracy_" + ref, ref,
                          haplotag_pct, "Read haplotag agreement\nwith truth (%)", "platform")


# ------------------------------------------------ accuracy against contiguity: 3

def fig3(db, ctx, ref, caller="truth"):
    sizes = [5 + 7 * i for i in range(len(st.DEPTHS))]
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), sharex=True)
    xs_all, ys_all = [], []
    for ax, platform in zip(axes, st.PLATFORMS):
        lines = 0
        for phaser in st.PHASERS:
            pts = []
            for depth, size in zip(st.DEPTHS, sizes):
                key = (ref, caller, platform, depth, phaser)
                x, y = positive(n50_mb(db, key)), positive(switch_pct(db, key))
                if finite(x) and finite(y):
                    pts.append((x, y, size))
            if not pts:
                continue
            xs, ys, ss = zip(*pts)
            ax.plot(xs, ys, color=st.PHASER_COLOR[phaser], linewidth=0.8, zorder=2)
            ax.scatter(xs, ys, s=ss, color=st.PHASER_COLOR[phaser],
                       marker=st.PHASER_MARKER[phaser], edgecolor="white",
                       linewidth=0.5, zorder=3)
            xs_all.extend(xs)
            ys_all.extend(ys)
            lines += 1
        ax.set_title(st.PLATFORM_LABEL[platform])
        ax.set_xlabel("Phase block N50 (Mb)")
        y_grid(ax)
        if not lines:
            empty_panel(ax)
    if not xs_all:
        plt.close(fig)
        raise NoData("no %s rows for %s" % (caller, ref))

    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    log_ticks(axes[0].xaxis, min(xs_all), max(xs_all))
    log_ticks(axes[0].yaxis, min(ys_all), max(ys_all))
    axes[0].set_ylabel("Switch error rate (%)")

    coverage = [Line2D([], [], color="none", marker="o", markerfacecolor=st.MUTED,
                       markeredgecolor="white", markersize=math.sqrt(s),
                       label="%d×" % d)
                for d, s in ((st.DEPTHS[0], sizes[0]), (st.DEPTHS[-1], sizes[-1]))]
    top_legend(legend_ax, phaser_handles() + coverage)
    return save(fig, ctx, "fig3_accuracy_vs_contiguity_" + ref)


# ------------------------------------------------------- kinds of error: 4, S1

def phaser_axis(ax):
    ax.set_xticks(range(len(st.PHASERS)))
    ax.set_xticklabels([st.PHASER_LABEL[p] for p in st.PHASERS],
                       rotation=35, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.6, len(st.PHASERS) - 0.4)
    ax.tick_params(axis="x", length=0)


def fig4(db, ctx, ref, caller="truth"):
    depth = ctx.depth
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), height=st.PANEL_HEIGHT + 6)
    bars = 0
    for ax, platform in zip(axes, st.PLATFORMS):
        here = 0
        for i, phaser in enumerate(st.PHASERS):
            rates = db.switch_flip((ref, caller, platform, depth, phaser))
            if rates is None:
                continue
            switches, flips = 100 * rates[0], 100 * rates[1]
            color = st.PHASER_COLOR[phaser]
            ax.bar(i, switches, width=0.62, color=color, edgecolor="white",
                   linewidth=0.8, zorder=3)
            ax.bar(i, flips, bottom=switches, width=0.62, color=color, alpha=0.45,
                   edgecolor="white", linewidth=0.8, zorder=3)
            here += 1
        ax.set_title(st.PLATFORM_LABEL[platform])
        phaser_axis(ax)
        y_grid(ax)
        if not here:
            empty_panel(ax)
        bars += here
    if not bars:
        plt.close(fig)
        raise NoData("no switch/flip counts for %s at %d×" % (ref, depth))

    axes[0].set_ylabel("Errors per 100 het SNP\npairs, %d×" % depth)
    axes[0].set_ylim(bottom=0)
    top_legend(legend_ax, [
        Patch(facecolor=st.MUTED, edgecolor="white", label="Switch error (wrong stretch)"),
        Patch(facecolor=st.MUTED, alpha=0.45, edgecolor="white", label="Flip error (single site)"),
    ])
    return save(fig, ctx, "fig4_switch_vs_flip_%dx_%s" % (depth, ref))


def figS1(db, ctx, ref, caller="truth"):
    """No replicate downsamples exist, so the spread over the 22 autosomes is
    the only variability these numbers come with."""
    depth = ctx.depth
    rng = np.random.RandomState(0)
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), height=st.PANEL_HEIGHT + 6)
    values = []
    for ax, platform in zip(axes, st.PLATFORMS):
        here = 0
        for i, phaser in enumerate(st.PHASERS):
            rates = [100 * r for r in
                     db.chromosome_switch_rates((ref, caller, platform, depth, phaser)) if r > 0]
            if not rates:
                continue
            color = st.PHASER_COLOR[phaser]
            ax.boxplot([rates], positions=[i], widths=0.56, patch_artist=True,
                       showfliers=False, zorder=3,
                       boxprops=dict(facecolor=color, alpha=0.35, edgecolor=color, linewidth=0.8),
                       medianprops=dict(color=st.INK, linewidth=1.0),
                       whiskerprops=dict(color=color, linewidth=0.8),
                       capprops=dict(color=color, linewidth=0.8))
            ax.scatter(i + rng.uniform(-0.17, 0.17, len(rates)), rates, s=5,
                       color=color, edgecolor="white", linewidth=0.3, zorder=4)
            values.extend(rates)
            here += 1
        ax.set_title(st.PLATFORM_LABEL[platform])
        phaser_axis(ax)
        y_grid(ax)
        if not here:
            empty_panel(ax)
    if not values:
        plt.close(fig)
        raise NoData("no per-chromosome compare for %s at %d×" % (ref, depth))

    axes[0].set_yscale("log")
    log_ticks(axes[0].yaxis, min(values), max(values))
    axes[0].set_ylabel("Switch error rate per\nautosome (%%), %d×" % depth)
    return save(fig, ctx, "figS1_switch_rate_per_chromosome_%dx_%s" % (depth, ref))


# ------------------------------------------------------------------- heatmaps

def heatmap(ax, matrix, row_labels, col_labels, row_groups, col_groups, **image):
    """Rows and columns come in labelled groups; a cell with no value is
    hatched, never painted as a low value."""
    arr = np.ma.masked_invalid(np.asarray(matrix, dtype=float))
    nrow, ncol = arr.shape
    ax.add_patch(Rectangle((-0.5, -0.5), ncol, nrow, facecolor="white",
                           edgecolor=st.MISSING, hatch="//////", linewidth=0, zorder=0))
    im = ax.imshow(arr, aspect="auto", interpolation="nearest", zorder=1, **image)

    row_edges = set(np.cumsum([n for _, n in row_groups])[:-1])
    col_edges = set(np.cumsum([n for _, n in col_groups])[:-1])
    for r in range(1, nrow):
        ax.axhline(r - 0.5, color="white", linewidth=2.5 if r in row_edges else 0.5, zorder=2)
    for c in range(1, ncol):
        ax.axvline(c - 0.5, color="white", linewidth=2.5 if c in col_edges else 0.5, zorder=2)

    ax.set_xticks(range(ncol))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(row_labels)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    start = 0
    for label, n in col_groups:
        ax.annotate(label, xy=(start + (n - 1) / 2, 1), xycoords=ax.get_xaxis_transform(),
                    xytext=(0, 3), textcoords="offset points", ha="center", va="bottom")
        start += n
    start = 0
    widest = max(len(str(l)) for l in row_labels)
    for label, n in row_groups:
        ax.annotate(label, xy=(0, start + (n - 1) / 2), xycoords=ax.get_yaxis_transform(),
                    xytext=(-(widest * 3.6 + 12), 0), textcoords="offset points",
                    rotation=90, ha="center", va="center")
        start += n
    return im


def log_norm(values):
    lo, hi = min(values), max(values)
    if hi <= lo:
        lo, hi = lo / 2, hi * 2
    return LogNorm(vmin=lo, vmax=hi)


def colorbar(fig, im, ax, label):
    bar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.015, aspect=22)
    bar.set_label(label)
    bar.outline.set_visible(False)
    bar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    bar.ax.yaxis.set_minor_locator(NullLocator())
    bar.ax.tick_params(length=2)
    return bar


def grid_cells(ref):
    """caller x phaser down, platform x depth across - the whole grid of one
    reference on one page."""
    rows = [(c, p) for c in st.CALLERS for p in st.PHASERS]
    cols = [(pl, d) for pl in st.PLATFORMS for d in st.DEPTHS]
    row_labels = [st.PHASER_LABEL[p] for _, p in rows]
    col_labels = [str(d) for _, d in cols]
    row_groups = [(st.CALLER_LABEL[c], len(st.PHASERS)) for c in st.CALLERS]
    col_groups = [(st.PLATFORM_LABEL[pl], len(st.DEPTHS)) for pl in st.PLATFORMS]
    keys = [[(ref, c, pl, d, p) for pl, d in cols] for c, p in rows]
    return keys, row_labels, col_labels, row_groups, col_groups


def fig5a(db, ctx, ref):
    keys, row_labels, col_labels, row_groups, col_groups = grid_cells(ref)
    matrix = [[positive(switch_pct(db, k)) for k in row] for row in keys]
    values = [v for row in matrix for v in row if finite(v)]
    if not values:
        raise NoData("no rows for %s" % ref)

    fig, ax = plt.subplots(figsize=(st.FULL_WIDTH * st.MM, 92 * st.MM), constrained_layout=True)
    im = heatmap(ax, matrix, row_labels, col_labels, row_groups, col_groups,
                 cmap=st.SEQUENTIAL, norm=log_norm(values))
    ax.set_xlabel("Coverage (×)")
    colorbar(fig, im, ax, "Switch error rate (%)")
    return save(fig, ctx, "fig5a_switch_rate_all_callers_" + ref)


def figS0(db, ctx, ref):
    keys, row_labels, col_labels, row_groups, col_groups = grid_cells(ref)
    code = {"done": 1.0, "locked": 0.0}
    matrix = [[code.get(db.status.get(k), NAN) for k in row] for row in keys]

    fig = plt.figure(figsize=(st.FULL_WIDTH * st.MM, 96 * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1, 0.05])
    ax = fig.add_subplot(grid[0])
    legend_ax = fig.add_subplot(grid[1])
    legend_ax.axis("off")
    heatmap(ax, matrix, row_labels, col_labels, row_groups, col_groups,
            cmap=ListedColormap(["#B9C4CF", "#4A5C70"]), vmin=0, vmax=1)
    ax.set_xlabel("Coverage (×)")
    legend_ax.legend(handles=[
        Patch(facecolor="#4A5C70", label="Computed (%d)" % db.count(ref, "done")),
        Patch(facecolor="#B9C4CF", label="Being written, locked (%d)" % db.count(ref, "locked")),
        Patch(facecolor="white", edgecolor=st.MISSING, hatch="//////", linewidth=0,
              label="Not computed (%d)" % db.count(ref, "missing")),
    ], loc="center", ncol=3, borderaxespad=0)
    return save(fig, ctx, "figS0_grid_completeness_" + ref)


def fig8(db, ctx, ref, caller="truth"):
    depth = ctx.depth
    strata = db.strata_names()
    cols = [(pl, p) for pl in st.PLATFORMS for p in st.PHASERS]
    matrix = []
    for stratum in strata:
        line = []
        for platform, phaser in cols:
            row = db.stratum((ref, caller, platform, depth, phaser), stratum)
            rate = ratio(row["switches"], row["assessed"]) if row else NAN
            line.append(positive(100 * rate))
        matrix.append(line)
    kept = [i for i, line in enumerate(matrix) if any(finite(v) for v in line)]
    if not kept:
        raise NoData("no metrics_strata rows for %s at %d× - regions.sh has not run" % (ref, depth))
    matrix = [matrix[i] for i in kept]
    strata = [strata[i] for i in kept]
    values = [v for line in matrix for v in line if finite(v)]

    fig, ax = plt.subplots(
        figsize=(st.FULL_WIDTH * st.MM, (34 + 5.5 * len(strata)) * st.MM), constrained_layout=True)
    im = heatmap(ax, matrix, [stratum_label(s) for s in strata],
                 [st.PHASER_LABEL[p] for _, p in cols],
                 [("", len(strata))],
                 [(st.PLATFORM_LABEL[pl], len(st.PHASERS)) for pl in st.PLATFORMS],
                 cmap=st.SEQUENTIAL, norm=log_norm(values))
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
    colorbar(fig, im, ax, "Switch error rate (%%), %d×" % depth)
    return save(fig, ctx, "fig8_switch_rate_by_stratum_%dx_%s" % (depth, ref))


STRATUM_LABEL = {
    "notinalldifficultregions": "Easy regions",
    "alldifficultregions": "All difficult regions",
    "segdups": "Segmental duplications",
    "lowmappabilityall": "Low mappability",
    "AllTandemRepeatsandHomopolymers_slop5": "Tandem repeats and homopolymers",
    "satellites_slop5": "Satellites",
    "MHC": "MHC",
    "chrX_nonPAR": "chrX non-PAR",
}


def stratum_label(name):
    return STRATUM_LABEL.get(name, name.replace("_", " "))


# ------------------------------------------- callers over one denominator: 5b

def fig5b(db, ctx, ref):
    """Both ends of a bar are fractions of the same truth het SNPs: the open end
    is what the caller found, the filled end what was also phased correctly. The
    gap is what phasing lost; the distance of the open end from 100 is what
    calling lost."""
    depth = ctx.depth
    rows = [(c, p) for c in st.CALLERS for p in st.PHASERS]
    gap = 0.9
    ypos = [i + gap * st.CALLERS.index(c) for i, (c, _) in enumerate(rows)]

    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), height=88, sharex=True)
    drawn = 0
    for ax, platform in zip(axes, st.PLATFORMS):
        here = 0
        for y, (caller, phaser) in zip(ypos, rows):
            key = (ref, caller, platform, depth, phaser)
            called = 100 * ratio(db.get(key, "shared"), db.get(key, "truth_total"))
            correct = 100 * db.get(key, "recall_phased")
            if not (finite(called) and finite(correct)):
                continue
            color = st.PHASER_COLOR[phaser]
            ax.plot([correct, called], [y, y], color=color, linewidth=1.2, zorder=2)
            ax.plot(called, y, marker="o", markerfacecolor="white", markeredgecolor=color,
                    markeredgewidth=0.9, zorder=3)
            ax.plot(correct, y, marker="o", markerfacecolor=color, markeredgecolor="white", zorder=4)
            here += 1
        ax.set_title(st.PLATFORM_LABEL[platform])
        ax.set_xlabel("Truth het SNPs (%)")
        ax.grid(axis="x", color=st.GRID, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
        if not here:
            empty_panel(ax)
        drawn += here
    if not drawn:
        plt.close(fig)
        raise NoData("no rows for %s at %d×" % (ref, depth))

    first = axes[0]
    first.set_yticks(ypos)
    first.set_yticklabels([st.PHASER_LABEL[p] for _, p in rows])
    first.set_ylim(ypos[-1] + 0.7, ypos[0] - 0.7)
    for i, caller in enumerate(st.CALLERS):
        block = ypos[i * len(st.PHASERS):(i + 1) * len(st.PHASERS)]
        first.annotate(st.CALLER_LABEL[caller], xy=(0, sum(block) / len(block)),
                       xycoords=first.get_yaxis_transform(), xytext=(-46, 0),
                       textcoords="offset points", rotation=90, ha="center", va="center")
    top_legend(legend_ax, [
        Line2D([], [], color="none", marker="o", markerfacecolor="white",
               markeredgecolor=st.MUTED, markeredgewidth=0.9, label="Called as het SNP"),
        Line2D([], [], color="none", marker="o", markerfacecolor=st.MUTED,
               markeredgecolor="white", label="Called and phased correctly, %d×" % depth),
    ])
    return save(fig, ctx, "fig5b_truth_recovered_%dx" % depth + "_" + ref)


# ------------------------------------------------------ GRCh38 against T2T: 6

def reference_scatter(db, ctx, name, value, label, log):
    a, b = st.REFERENCES[0], st.REFERENCES[1]
    fig = plt.figure(figsize=(118 * st.MM, 72 * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=[1, 0.42])
    ax = fig.add_subplot(grid[0, 0])
    legend_ax = fig.add_subplot(grid[0, 1])
    legend_ax.axis("off")

    seen = []
    for platform in st.PLATFORMS:
        for phaser in st.PHASERS:
            xs, ys = [], []
            for caller in st.CALLERS:
                for depth in st.DEPTHS:
                    x = value(db, (a, caller, platform, depth, phaser))
                    y = value(db, (b, caller, platform, depth, phaser))
                    if log:
                        x, y = positive(x), positive(y)
                    if finite(x) and finite(y):
                        xs.append(x)
                        ys.append(y)
            if xs:
                ax.scatter(xs, ys, s=13, color=st.PLATFORM_COLOR[platform],
                           marker=st.PHASER_MARKER[phaser], edgecolor="white",
                           linewidth=0.4, alpha=0.9, zorder=3)
                seen.extend(xs + ys)
    if not seen:
        plt.close(fig)
        raise NoData("no run is computed on both %s and %s" % (a, b))

    lo, hi = min(seen), max(seen)
    if log:
        lo, hi = lo / 1.4, hi * 1.4
        ax.set_xscale("log")
        ax.set_yscale("log")
        log_ticks(ax.xaxis, lo, hi)
        log_ticks(ax.yaxis, lo, hi)
    else:
        pad = 0.04 * (hi - lo) or 1
        lo, hi = lo - pad, hi + pad
    ax.plot([lo, hi], [lo, hi], color=st.MUTED, linewidth=0.6, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("%s, %s" % (label, st.REFERENCE_LABEL[a]))
    ax.set_ylabel("%s, %s" % (label, st.REFERENCE_LABEL[b]))
    ax.grid(color=st.GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

    first = legend_ax.legend(handles=platform_handles(line=False), title="Platform",
                             loc="upper left", borderaxespad=0)
    legend_ax.add_artist(first)
    second = legend_ax.legend(
        handles=[Line2D([], [], color="none", marker=st.PHASER_MARKER[p],
                        markerfacecolor=st.MUTED, markeredgecolor="white",
                        label=st.PHASER_LABEL[p]) for p in st.PHASERS],
        title="Phaser", loc="upper left", bbox_to_anchor=(0, 0.6), borderaxespad=0)
    for legend in (first, second):
        if hasattr(legend, "set_alignment"):
            legend.set_alignment("left")
    return save(fig, ctx, name)


def recall_pct(db, key):
    return 100 * db.get(key, "recall_phased")


def fig6a(db, ctx):
    return reference_scatter(db, ctx, "fig6a_reference_switch_rate",
                             switch_pct, "Switch error rate (%)", log=True)


def fig6b(db, ctx):
    return reference_scatter(db, ctx, "fig6b_reference_block_n50",
                             n50_mb, "Block N50 (Mb)", log=True)


def fig6c(db, ctx):
    return reference_scatter(db, ctx, "fig6c_reference_recall_phased",
                             recall_pct, "Correctly phased truth het SNPs (%)", log=False)


# (id, function, drawn once per reference?)
FIGURES = [
    ("figS0", figS0, True),
    ("fig1a", fig1a, True),
    ("fig1b", fig1b, True),
    ("fig1c", fig1c, True),
    ("fig2a", fig2a, True),
    ("fig2b", fig2b, True),
    ("fig3", fig3, True),
    ("fig4", fig4, True),
    ("fig5a", fig5a, True),
    ("fig5b", fig5b, True),
    ("fig6a", fig6a, False),
    ("fig6b", fig6b, False),
    ("fig6c", fig6c, False),
    ("fig7a", fig7a, True),
    ("fig7b", fig7b, True),
    ("fig8", fig8, True),
    ("figS1", figS1, True),
]
