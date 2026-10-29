"""The figures. Every function draws one subfigure into its own PDF and raises
NoData when the metrics it needs are not there yet.

Which denominators are shared, and so what may be compared with what:
  phasers on the GIAB truth input   every phaser gets the same sites, so
                                    switch rate, N50 and phased fraction compare
                                    directly (fig1-fig4, fig7, figS1)
  callers                           each caller hands the phaser its own sites,
                                    so its switch rate is over a different set
                                    of sites and is never shown caller against
                                    caller on its own. Callers are compared on
                                    what is divided by truth_total, one constant
                                    per reference: FNR, and truth het SNPs
                                    phased correctly (fig5a, fig5c, fig5d). FDR
                                    is over the caller's own calls, as precision
                                    always is (fig5b). The switch rate appears
                                    only against the share of the truth it was
                                    measured on (fig5e)
  references                        truth_total differs between GRCh38 and T2T,
                                    so only rates are put side by side (fig6)
"""

import math
import os
import textwrap

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
    """Plain numbers on a log axis, no 10^x - it would pull in the math font.
    The wider the range, the fewer labels per decade, so neighbours never touch."""
    span = hi / lo if finite(lo) and finite(hi) and lo > 0 else 1
    subs = (1,) if span > 300 else (1, 3) if span > 40 else (1, 2, 5)
    axis.set_major_locator(LogLocator(base=10, subs=subs, numticks=15))
    axis.set_major_formatter(FuncFormatter(lambda v, _: tick_number(v)))
    axis.set_minor_locator(NullLocator())


def tick_number(v):
    """0.005, 0.3, 20, 5k, 2M - never 1e+06."""
    if v >= 1e6:
        return "%gM" % (v / 1e6)
    if v >= 1e4:
        return "%gk" % (v / 1e3)
    return "%g" % v


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


def caller_handles(callers, line=True):
    return [series_handle(st.CALLER_COLOR[c], st.CALLER_MARKER[c], st.CALLER_LABEL[c], line)
            for c in callers]


def top_legend(legend_ax, handles, **kwargs):
    kwargs.setdefault("ncol", len(handles))
    legend_ax.legend(handles=handles, loc="center", borderaxespad=0, **kwargs)


# ------------------------------------------------- lines over coverage: 1, 2, 7

def coverage_lines(db, ctx, name, ref, value, ylabel, by, log=False, caller="truth"):
    """by='platform': a panel per platform, a line per phaser.
       by='phaser':   a panel per phaser, a line per platform.
       by='caller':   a panel per platform, a line per caller - for what the
                      caller alone decides, so the phaser is not a dimension."""
    if by == "caller":
        panels, series = st.PLATFORMS, [c for c in st.CALLERS if c != "truth"]
        panel_label, color, marker = st.PLATFORM_LABEL, st.CALLER_COLOR, st.CALLER_MARKER
        handles = caller_handles(series)
        # a caller's FNR and FDR from bin/calls_fdr.sh: inside the confident
        # regions on chr1-22, the same for every phaser of the cell
        key = lambda panel, s, depth: ("call", ref, s, panel, depth)
    elif by == "platform":
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
    """NG50 over the autosomes from plots/ng50_autosomes.py. whatshap's own
    block_n50 is not a fallback: on GRCh38 it counts chrX into the genome and
    on T2T it does not, so the two would not compare."""
    return db.get(key, "block_ng50") / 1e6


def truth_phased_pct(db, key):
    """Truth het SNPs phased in a block of two or more, out of all of them -
    one denominator for every run of the reference."""
    return 100 * ratio(db.get(key, "covered"), db.get(key, "truth_total"))


def tagged_pct(db, key):
    """Reads on chr1-22 tagged in a block of two SNPs or more (bin/autosomes.sh)."""
    return 100 * db.get(key, "tagging_completeness_min")


def tagged_bases_pct(db, key):
    return 100 * db.get(key, "bases_tagged_fraction")


def read_precision_pct(db, key):
    return 100 * db.get(key, "read_precision")


def read_recall_pct(db, key):
    return 100 * db.get(key, "read_recall")


def fnr_pct(db, key):
    if key and key[0] == "call":
        return 100 * db.call(key[1], key[2], key[3], key[4], "fnr")
    return 100 * db.get(key, "fnr")


def fdr_pct(db, key):
    if key and key[0] == "call":
        return 100 * db.call(key[1], key[2], key[3], key[4], "fdr")
    return 100 * db.get(key, "fdr")


def fig1a(db, ctx, ref):
    return coverage_lines(db, ctx, "fig1a_switch_rate_by_platform_" + ref, ref,
                          switch_pct, "Switch error rate (%)", "platform", log=True)


def fig1b(db, ctx, ref):
    return coverage_lines(db, ctx, "fig1b_block_n50_by_platform_" + ref, ref,
                          n50_mb, "Phase block NG50 (Mb)", "platform", log=True)


def fig1c(db, ctx, ref):
    return coverage_lines(db, ctx, "fig1c_truth_phased_by_platform_" + ref, ref,
                          truth_phased_pct, "Truth het SNPs phased (%)", "platform")


def fig2a(db, ctx, ref):
    return coverage_lines(db, ctx, "fig2a_switch_rate_by_phaser_" + ref, ref,
                          switch_pct, "Switch error rate (%)", "phaser", log=True)


def fig2b(db, ctx, ref):
    return coverage_lines(db, ctx, "fig2b_block_n50_by_phaser_" + ref, ref,
                          n50_mb, "Phase block NG50 (Mb)", "phaser", log=True)


def fig7a(db, ctx, ref):
    """By bases: a read without a het SNP on it cannot be tagged, and counting
    reads weighs the many short ONT reads like long ones (figS4 has the read
    count)."""
    return coverage_lines(db, ctx, "fig7a_bases_tagged_" + ref, ref,
                          tagged_bases_pct, "Read bases assigned\na haplotype (%)", "platform")


def figS4(db, ctx, ref):
    return coverage_lines(db, ctx, "figS4_reads_tagged_" + ref, ref,
                          tagged_pct, "Reads assigned a haplotype (%)", "platform")


def fig7b(db, ctx, ref):
    return coverage_lines(db, ctx, "fig7b_read_precision_" + ref, ref,
                          read_precision_pct,
                          "Reads on the truth haplotype,\nof reads tagged by both (%)", "platform")


def fig7c(db, ctx, ref):
    return coverage_lines(db, ctx, "fig7c_read_recall_" + ref, ref,
                          read_recall_pct,
                          "Reads on the truth haplotype,\nof reads tagged by the truth (%)", "platform")


# ------------------------------------------------ accuracy against contiguity: 3

def better_top_right(axes):
    """Fewer errors go up, so the best run sits in the top right corner of
    every panel; the axes share y, so inverting one inverts them all."""
    axes[0].invert_yaxis()


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
        ax.set_xlabel("Phase block NG50 (Mb)")
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
    better_top_right(axes)
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

    axes[0].set_ylabel("Errors per 100 pairs\nof het SNPs, %d×" % depth)
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
    every = [100 * r for platform in st.PLATFORMS for phaser in st.PHASERS
             for r in db.chromosome_switch_rates((ref, caller, platform, depth, phaser))]
    if not every:
        raise NoData("no per-chromosome compare for %s at %d×" % (ref, depth))
    # an autosome without a switch is 0, which a log axis cannot show: it is
    # drawn on a line of its own below the smallest rate, labelled 0, and kept
    # in the box statistics there instead of being dropped
    positive_rates = [r for r in every if r > 0]
    floor = (min(positive_rates) if positive_rates else 1e-4) / 3
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), height=st.PANEL_HEIGHT + 6)
    values = []
    for ax, platform in zip(axes, st.PLATFORMS):
        here = 0
        for i, phaser in enumerate(st.PHASERS):
            rates = [max(100 * r, floor) for r in
                     db.chromosome_switch_rates((ref, caller, platform, depth, phaser))]
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
    log_ticks(axes[0].yaxis, floor, max(values))
    ticks = [t for t in axes[0].get_yticks() if floor * 1.5 < t <= max(values) * 1.5]
    axes[0].set_yticks([floor] + ticks)
    axes[0].set_yticklabels(["0"] + ["%g" % t for t in ticks])
    axes[0].set_ylim(floor / 1.6, max(values) * 1.6)
    for ax in axes:
        ax.axhline(floor * 1.25, color=st.GRID, linewidth=0.6, linestyle=(0, (2, 2)), zorder=1)
    axes[0].set_ylabel("Switch error rate per\nautosome (%%), %d×" % depth)
    return save(fig, ctx, "figS1_switch_rate_per_chromosome_%dx_%s" % (depth, ref))


# ------------------------------------------------------------------- heatmaps

def heatmap(ax, matrix, row_labels, col_labels, row_groups, col_groups, label=None, **image):
    """Rows and columns come in labelled groups; a cell with no value is
    hatched, never painted as a low value. label(value) -> text writes the
    value into each cell, in ink or white, whichever reads on that colour."""
    arr = np.ma.masked_invalid(np.asarray(matrix, dtype=float))
    nrow, ncol = arr.shape
    ax.add_patch(Rectangle((-0.5, -0.5), ncol, nrow, facecolor="white",
                           edgecolor=st.MISSING, hatch="//////", linewidth=0, zorder=0))
    im = ax.imshow(arr, aspect="auto", interpolation="nearest", zorder=1, **image)

    if label is not None:
        for r in range(nrow):
            for c in range(ncol):
                if arr.mask[r, c] if np.ndim(arr.mask) else False:
                    continue
                red, green, blue, _ = im.cmap(im.norm(arr[r, c]))
                dark = 0.2126 * red + 0.7152 * green + 0.0722 * blue < 0.5
                ax.text(c, r, label(arr[r, c]), ha="center", va="center", zorder=3,
                        color="white" if dark else st.INK)

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
    group_labels_left(ax, row_groups)
    return im


def label_width_pt(ax, labels):
    """Width of the widest of these text objects in points, measured with the
    font in use - a guess from the number of letters breaks with a wider font."""
    renderer = ax.figure.canvas.get_renderer()
    widths = [t.get_window_extent(renderer).width for t in labels if t.get_text()]
    return max(widths) * 72 / ax.figure.dpi if widths else 0


def group_labels_left(ax, groups, positions=None):
    """A vertical label for each group of rows, left of the row labels.
    positions: the y of every row; default 0, 1, 2, ..."""
    if not any(label for label, _ in groups):
        return
    offset = label_width_pt(ax, ax.get_yticklabels()) + 3 + ax.yaxis.get_tick_padding() + 4
    start = 0
    for label, n in groups:
        ys = positions[start:start + n] if positions is not None else [start, start + n - 1]
        ax.annotate(label, xy=(0, (ys[0] + ys[-1]) / 2), xycoords=ax.get_yaxis_transform(),
                    xytext=(-offset, 0), textcoords="offset points",
                    rotation=90, ha="center", va="bottom", rotation_mode="anchor")
        start += n


def log_norm(values):
    lo, hi = min(values), max(values)
    if hi <= lo:
        lo, hi = lo / 2, hi * 2
    return LogNorm(vmin=lo, vmax=hi)


def colorbar(fig, im, ax, label):
    bar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.015, aspect=22)
    bar.set_label(label)
    bar.outline.set_visible(False)
    if isinstance(im.norm, LogNorm):
        bar.ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
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
    return coverage_lines(db, ctx, "fig5a_caller_fnr_" + ref, ref, fnr_pct,
                          "Truth het SNPs missed (FNR, %)", "caller", log=True, caller="called")


def fig5b(db, ctx, ref):
    return coverage_lines(db, ctx, "fig5b_caller_fdr_" + ref, ref, fdr_pct,
                          "Called het SNPs not in truth\n(FDR, %)", "caller", log=True, caller="called")


def short_number(v):
    """Two significant figures, so a value fits its cell and a small one is
    not printed as 0.00."""
    return "%.0f" % v if v >= 10 else "%.2g" % v


def grid_heatmap(db, ctx, ref, name, value, bar_label):
    """The whole grid of one reference - caller x phaser down, platform x depth
    across - on one log colour scale, every cell labelled with its value.
    Only for metrics whose denominator does not depend on the caller."""
    keys, row_labels, col_labels, row_groups, col_groups = grid_cells(ref)
    matrix = [[positive(value(db, k)) for k in row] for row in keys]
    values = [v for row in matrix for v in row if finite(v)]
    if not values:
        raise NoData("no rows for %s" % ref)

    fig, ax = plt.subplots(figsize=(st.FULL_WIDTH * st.MM, 92 * st.MM), constrained_layout=True)
    im = heatmap(ax, matrix, row_labels, col_labels, row_groups, col_groups,
                 label=lambda v: f"{v:.2f}", cmap=st.SEQUENTIAL, norm=log_norm(values))
    ax.set_xlabel("Coverage (×)")
    colorbar(fig, im, ax, bar_label)
    return save(fig, ctx, name + "_" + ref)


def truth_lost_pct(db, key):
    """Truth het SNPs not called, or called and phased wrong. Kept shares sit
    between 95 and 99.9%, and only their complement on a log scale spreads
    them apart."""
    return 100 - recall_pct(db, key)


def phasing_lost_pct(db, key):
    """Truth het SNPs the caller did find but that were left unphased or
    phased wrong - what the phaser lost, on the truth denominator."""
    return 100 * ratio(db.get(key, "shared") - db.get(key, "correct"), db.get(key, "truth_total"))


def reads_lost_pct(db, key):
    """Reads tagged by the truth that this run left untagged or put on the
    wrong haplotype; every run of a cell has the same reads."""
    return 100 - read_recall_pct(db, key)


def fig5c(db, ctx, ref):
    return grid_heatmap(db, ctx, ref, "fig5c_truth_lost_all_callers", truth_lost_pct,
                        "Truth het SNPs not called or phased wrongly (%)")


def figS2a(db, ctx, ref):
    return grid_heatmap(db, ctx, ref, "figS2a_phasing_lost_all_callers", phasing_lost_pct,
                        "Called truth het SNPs left unphased\nor phased wrongly (% of truth)")


def figS2b(db, ctx, ref):
    return grid_heatmap(db, ctx, ref, "figS2b_block_ng50_all_callers", n50_mb,
                        "Phase block NG50 (Mb)")


def figS2c(db, ctx, ref):
    return grid_heatmap(db, ctx, ref, "figS2c_reads_lost_all_callers", reads_lost_pct,
                        "Truth-tagged reads not assigned\nto their haplotype (%)")


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


def platform_phaser_rows():
    """platform x phaser, without the pairs a phaser does not run on (HiPhase on ONT)"""
    return [(pl, p) for pl in st.PLATFORMS for p in st.PHASERS if st.applies(pl, p)]


def platform_groups(rows):
    """(label, number of rows) per platform, for the left-hand group labels"""
    return [(st.PLATFORM_LABEL[pl], sum(1 for r in rows if r[0] == pl))
            for pl in st.PLATFORMS if any(r[0] == pl for r in rows)]


def strata_heatmap(ctx, name, strata, value, bar_label, note=None, fmt=None, **image):
    """Platform x phaser down, one stratum per column: 16 rows of wide cells, so
    every value fits its cell. value(platform, phaser, stratum) -> number; NaN is
    hatched as not computed, 0 drawn in the lightest colour and written as 0.
    note(platform, phaser, stratum) -> text or None: a cell with a note is
    dotted, left uncoloured and carries the note instead of a value - for a
    count of zero or of too few errors to compare, so that it does not read as
    the best result of the figure. fmt(value) -> text overrides the
    two significant figures written in a cell."""
    rows = platform_phaser_rows()
    matrix = [[value(pl, p, s) for s in strata] for pl, p in rows]
    notes = [[note(pl, p, s) if note else None for s in strata] for pl, p in rows]
    kept = [j for j in range(len(strata))
            if any(finite(line[j]) or text[j] for line, text in zip(matrix, notes))]
    if not kept:
        return None
    strata = [strata[j] for j in kept]
    matrix = [[line[j] for j in kept] for line in matrix]
    notes = [[line[j] for j in kept] for line in notes]
    noted = [(i, j, text) for i, line in enumerate(notes) for j, text in enumerate(line) if text]
    for i, j, _ in noted:
        matrix[i][j] = NAN
    shown = [v for line in matrix for v in line if finite(v) and v > 0]
    floor = min(shown) / 2 if shown else 1.0
    matrix = [[floor if finite(v) and v <= 0 else v for v in line] for line in matrix]
    if "norm" not in image:
        image["norm"] = log_norm(shown + [floor])
    fig, ax = plt.subplots(figsize=((62 + 13 * len(strata)) * st.MM, 100 * st.MM),
                           constrained_layout=True)
    im = heatmap(ax, matrix, [st.PHASER_LABEL[p] for _, p in rows],
                 [stratum_label(s) for s in strata], platform_groups(rows),
                 [("", len(strata))],
                 label=lambda v: "0.00" if v <= floor else (fmt or (lambda x: f"{x:.2f}"))(v), **image)
    for i, j, text in noted:
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="white", edgecolor=st.MISSING,
                               hatch="....", linewidth=0, zorder=1.5))
        ax.text(j, i, text, ha="center", va="center", color=st.MUTED, style="italic", zorder=3)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
    colorbar(fig, im, ax, bar_label)
    return save(fig, ctx, name)


def fig8(db, ctx, ref, caller="truth"):
    depth = ctx.depth

    def rate(platform, phaser, stratum):
        row = db.stratum((ref, caller, platform, depth, phaser), stratum)
        return 100 * ratio(row["switches"], row["assessed"]) if row else NAN

    def note(platform, phaser, stratum):
        """no switch at all: written as 0 in an uncoloured cell"""
        row = db.stratum((ref, caller, platform, depth, phaser), stratum)
        return "0.00" if row and row["assessed"] and not row["switches"] else None

    out = strata_heatmap(ctx, "fig8_switch_rate_by_stratum_%dx_%s" % (depth, ref),
                         db.strata_names(), rate, "Switch error rate (%%), %d×" % depth,
                         note=note, cmap=st.SEQUENTIAL)
    if out is None:
        raise NoData("no metrics_strata rows for %s at %d× - regions.sh has not run" % (ref, depth))
    return out


STRATUM_LABEL = {
    "genome": "Whole genome\n(chr1-22)",
    "notinalldifficultregions": "Easy regions",
    "alldifficultregions": "All difficult regions",
    "segdups": "Segmental duplications",
    "lowmappabilityall": "Low mappability",
    "AllTandemRepeatsandHomopolymers_slop5": "Tandem repeats and\nhomopolymers",
    "satellites_slop5": "Satellites",
    "MHC": "MHC",
    "chrX_nonPAR": "chrX non-PAR",
}


def stratum_label(name):
    return STRATUM_LABEL.get(name, name.replace("_", " "))


# ------------------------------------- what the strata are: S3, one file each
#
# metrics_strata/regions_statistics_<ref>.tsv from regions_statistics.sh. Every
# number is inside the confident BED unless it says raw, i.e. over what the
# per-stratum metrics are computed on.

STRATA_DARK = "#4A5C70"
STRATA_LIGHT = "#B9C4CF"


def strata_rows(db, ref):
    stats = db.strata_statistics(ref)
    names = [s for s in db.strata_names() if s in stats]
    if not names:
        raise NoData("no metrics_strata/regions_statistics_%s.tsv - regions_statistics.sh has not run" % ref)
    return names, stats


def strata_axes(n, legend=0):
    """One column wide; legend: how many rows of legend go above the bars."""
    height = 18 + 7 * n + 4.5 * legend
    fig, ax = plt.subplots(figsize=(92 * st.MM, height * st.MM), constrained_layout=True)
    return fig, ax


def legend_above(ax, handles, ncol):
    """A legend over the plotting area itself, not over the whole figure, so it
    stays above the bars however wide the stratum names are."""
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=ncol,
              borderaxespad=0, columnspacing=1.0)


def strata_y(ax, names):
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([stratum_label(s) for s in names])
    ax.set_ylim(len(names) - 0.4, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=st.GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)


def thousands(v):
    return "{:,}".format(int(round(v)))


def figS3a(db, ctx, ref):
    """Each stratum before and after the confident BED: the dark bar is what
    the benchmark sees, the number how much of the stratum that is."""
    names, stats = strata_rows(db, ref)
    raw = [positive(stats[s]["bp_raw"] / 1e6) for s in names]
    kept = [positive(stats[s]["bp"] / 1e6) for s in names]
    values = [v for v in raw + kept if finite(v)]
    if not values:
        raise NoData("no stratum sizes for %s" % ref)
    left = 10 ** math.floor(math.log10(min(values) / 2))
    fig, ax = strata_axes(len(names), legend=1)
    for i, (r, k) in enumerate(zip(raw, kept)):
        if finite(r):
            ax.barh(i, r - left, left=left, height=0.7, color=STRATA_LIGHT, zorder=2)
        if finite(k):
            ax.barh(i, k - left, left=left, height=0.7, color=STRATA_DARK, zorder=3)
        if finite(r) and finite(k):
            ax.text(max(r, k) * 1.15, i, "%.0f%%" % (100 * k / r), va="center", color=st.INK)
    ax.set_xscale("log")
    ax.set_xlim(left, max(values) * 4)
    log_ticks(ax.xaxis, left, max(values) * 4)
    ax.set_xlabel("Size (Mb)")
    strata_y(ax, names)
    legend_above(ax, [Patch(facecolor=STRATA_LIGHT, label="Whole stratum"),
                      Patch(facecolor=STRATA_DARK, label="In confident regions")], ncol=2)
    return save(fig, ctx, "figS3a_strata_size_" + ref)


def strata_bars(ctx, ref, name, names, values, xlabel, text, log=False):
    """One bar per stratum, the bar's own value written at its end."""
    if not any(finite(v) for v in values):
        return None
    fig, ax = strata_axes(len(names))
    top = max(v for v in values if finite(v))
    low = min(v for v in values if finite(v) and v > 0) if log else 0
    left = 10 ** math.floor(math.log10(low / 2)) if log else 0
    for i, v in enumerate(values):
        if not finite(v):
            continue
        ax.barh(i, v - left, left=left, height=0.7, color=STRATA_DARK, zorder=3)
        ax.annotate(text(i, v), xy=(v, i), xytext=(3, 0), textcoords="offset points",
                    va="center", color=st.INK)
    if log:
        ax.set_xscale("log")
        ax.set_xlim(left, top * 30)
        log_ticks(ax.xaxis, left, top * 30)
    else:
        ax.set_xlim(0, top * 1.35)
    ax.set_xlabel(xlabel)
    strata_y(ax, names)
    return save(fig, ctx, name + "_" + ref)


def figS3b(db, ctx, ref):
    """Density of the truth het SNPs: where SNPs are sparse, blocks break."""
    names, stats = strata_rows(db, ref)
    out = strata_bars(ctx, ref, "figS3b_strata_het_snp_density", names,
                      [stats[s]["het_snps_per_mb"] for s in names],
                      "Truth het SNPs per Mb", lambda i, v: thousands(v))
    if out is None:
        raise NoData("no het SNP density for %s" % ref)
    return out


def figS3g(db, ctx, ref):
    """How many truth het SNPs each stratum has, and their share of all of
    them: a rate over a few thousand sites is not as firm as one over a
    million."""
    names, stats = strata_rows(db, ref)
    count = [positive(stats[s]["truth_het_snps"]) for s in names]
    total = NAN
    if "alldifficultregions" in stats and "notinalldifficultregions" in stats:
        total = stats["alldifficultregions"]["truth_het_snps"] + \
            stats["notinalldifficultregions"]["truth_het_snps"]

    def text(i, v):
        share = " (%s%%)" % short_number(100 * v / total) if finite(total) and total else ""
        return thousands(v) + share

    out = strata_bars(ctx, ref, "figS3g_strata_het_snp_count", names, count,
                      "Truth het SNPs", text, log=True)
    if out is None:
        raise NoData("no het SNP counts for %s" % ref)
    return out


# the scored part dark, the rest in hues of their own so neighbours never blur
COMPOSITION = [
    ("Het SNP, benchmarked", "#2C3A4B"),
    ("Het indel or multiallelic", "#C8775A"),
    ("Hom alt", "#9DB7CC"),
    ("Hemizygous", "#E3CF9A"),
    ("Other", "#E6E6E6"),
]


def figS3c(db, ctx, ref):
    """Every variant of the GIAB truth in the stratum, split by what the
    benchmark does with it. Only the darkest part is scored; the rest is
    variation the per-stratum numbers say nothing about."""
    names, stats = strata_rows(db, ref)
    fig, ax = strata_axes(len(names), legend=3)
    drawn = 0
    for i, s in enumerate(names):
        r = stats[s]
        total = r["src_variants"]
        if not (finite(total) and total > 0):
            continue
        bench = r["truth_het_snps"]
        parts = [bench, max(r["src_0_1"] + r["src_1_2"] - bench, 0),
                 r["src_1_1"], r["src_hemi"], r["src_other"]]
        left = 0.0
        for v, (_, color) in zip(parts, COMPOSITION):
            share = 100 * v / total if finite(v) else 0
            ax.barh(i, share, left=left, height=0.7, color=color,
                    edgecolor="white", linewidth=0.4, zorder=3)
            left += share
        ax.text(101.5, i, "%.0f%%" % (100 * bench / total), va="center", color=st.INK)
        drawn += 1
    if not drawn:
        plt.close(fig)
        raise NoData("no truth composition for %s" % ref)
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("GIAB truth variants in the stratum (%)")
    strata_y(ax, names)
    legend_above(ax, [Patch(facecolor=c, edgecolor=st.GRID, linewidth=0.4, label=l)
                      for l, c in COMPOSITION], ncol=2)
    return save(fig, ctx, "figS3c_strata_composition_" + ref)


def figS3d(db, ctx, ref):
    """Phase blocks break where het SNPs are sparse: the typical gap."""
    names, stats = strata_rows(db, ref)
    out = strata_bars(ctx, ref, "figS3d_strata_snp_spacing", names,
                      [positive(stats[s]["median_inter_snp_dist"]) for s in names],
                      "Median distance to the next\ntruth het SNP (bp)",
                      lambda i, v: thousands(v), log=True)
    if out is None:
        raise NoData("no SNP spacing for %s" % ref)
    return out


def figS3e(db, ctx, ref):
    """SNPs this close share most reads, and a read error there touches both."""
    names, stats = strata_rows(db, ref)
    out = strata_bars(ctx, ref, "figS3e_strata_close_snps", names,
                      [100 * stats[s]["frac_snp_within_50bp"] for s in names],
                      "Truth het SNPs less than 50 bp\nfrom the next one (%)",
                      lambda i, v: short_number(v) + "%")
    if out is None:
        raise NoData("no close SNPs for %s" % ref)
    return out


def figS3f(db, ctx, ref):
    """How much the strata are the same sequence: Jaccard index of each pair
    inside the confident regions. Nested strata are not independent results."""
    names, _ = strata_rows(db, ref)
    pairs = db.strata_overlap(ref)
    names = [s for s in names if any((s, t) in pairs for t in names)]
    if not names:
        raise NoData("no metrics_strata/jaccard_%s.tsv" % ref)
    matrix = [[1.0 if a == b else pairs.get((a, b), NAN) for b in names] for a in names]
    n = len(names)
    fig, ax = plt.subplots(figsize=((60 + 9 * n) * st.MM, (34 + 7.5 * n) * st.MM),
                           constrained_layout=True)
    labels = [stratum_label(s) for s in names]
    im = heatmap(ax, matrix, labels, labels, [("", n)], [("", n)],
                 label=lambda v: "%.2f" % v, cmap=st.SEQUENTIAL, vmin=0, vmax=1)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
    colorbar(fig, im, ax, "Jaccard index")
    return save(fig, ctx, "figS3f_strata_overlap_" + ref)


# ------------------------------------------- callers over one denominator: 5b

def fig5d(db, ctx, ref):
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
    group_labels_left(first, [(st.CALLER_LABEL[c], len(st.PHASERS)) for c in st.CALLERS], ypos)
    top_legend(legend_ax, [
        Line2D([], [], color="none", marker="o", markerfacecolor="white",
               markeredgecolor=st.MUTED, markeredgewidth=0.9, label="Called as het SNP"),
        Line2D([], [], color="none", marker="o", markerfacecolor=st.MUTED,
               markeredgecolor="white", label="Called and phased correctly, %d×" % depth),
    ])
    return save(fig, ctx, "fig5d_truth_recovered_%dx" % depth + "_" + ref)


def fig5e(db, ctx, ref):
    """The switch rate of a caller is measured on the sites that caller found
    and the phaser phased, so it is drawn against how much of the truth that
    was. A caller that skips the hard sites moves down and to the left, and
    the move to the left is visible."""
    depth = ctx.depth
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), sharex=True)
    ys_all = []
    for ax, platform in zip(axes, st.PLATFORMS):
        here = 0
        for caller in st.CALLERS:
            for phaser in st.PHASERS:
                key = (ref, caller, platform, depth, phaser)
                x = 100 * ratio(db.get(key, "covered"), db.get(key, "truth_total"))
                y = positive(switch_pct(db, key))
                if not (finite(x) and finite(y)):
                    continue
                ax.scatter([x], [y], s=20, color=st.CALLER_COLOR[caller],
                           marker=st.PHASER_MARKER[phaser], edgecolor="white",
                           linewidth=0.5, zorder=3)
                ys_all.append(y)
                here += 1
        ax.set_title(st.PLATFORM_LABEL[platform])
        ax.set_xlabel("Truth het SNPs phased (%)")
        ax.grid(color=st.GRID, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        if not here:
            empty_panel(ax)
    if not ys_all:
        plt.close(fig)
        raise NoData("no rows for %s at %d×" % (ref, depth))

    axes[0].set_yscale("log")
    log_ticks(axes[0].yaxis, min(ys_all), max(ys_all))
    better_top_right(axes)
    axes[0].set_ylabel("Switch error rate on\nthose sites (%%), %d×" % depth)
    shapes = [Line2D([], [], color="none", marker=st.PHASER_MARKER[p], markerfacecolor="white",
                     markeredgecolor=st.INK, markeredgewidth=0.7, label=st.PHASER_LABEL[p])
              for p in st.PHASERS]
    # colour says caller here and shape says phaser, so the caller entries are
    # plain dots rather than the callers' own line markers
    colours = [Line2D([], [], color="none", marker="o", markerfacecolor=st.CALLER_COLOR[c],
                      markeredgecolor="white", label=st.CALLER_LABEL[c]) for c in st.CALLERS]
    top_legend(legend_ax, colours + shapes)
    return save(fig, ctx, "fig5e_switch_rate_vs_truth_phased_%dx" % depth + "_" + ref)


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
            for caller in ("truth",):     # see reference_pairs
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


# A difference of 10% between the references is invisible on a log-log plot
# that spans two decades, so the ratio T2T / GRCh38 of the same run gets plots
# of its own: a box per group, its runs as points over it. Every file is one
# panel; all files of one metric share their axis limits, so they can be put
# side by side.
#   scatter_all          every run, colour and shape are the phaser
#   scatter_ONT/PacBio   one technology, colour is the platform, shape the phaser
#   scatter_by_platform  every run, colour and shape are the platform
#   ratio_by_phaser      a box per phaser, all platforms together
#   ratio_by_technology  a box per phaser, ONT and PacBio apart
#   ratio_by_platform    a box per platform

RATIO_TICKS = [0.25, 0.5, 0.67, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 2, 4]
TECHNOLOGIES = [("ONT", ["ont4k", "ont5k"]), ("PacBio", ["hifi", "revio"])]
BAND = 1.1  # the shaded band: T2T within 10% of GRCh38


def reference_pairs(db, value):
    """(GRCh38, T2T) for every run computed on both, by (platform, phaser);
    the depths of a cell are its points. Truth input only: the two references
    were aligned differently (GRCh38 BAMs as their providers aligned them, T2T
    by us with minimap2), and Clair3 reacts to that; with the truth as input
    no caller stands between the reads and the comparison."""
    a, b = st.REFERENCES[0], st.REFERENCES[1]
    pairs = {}
    for platform in st.PLATFORMS:
        for phaser in st.PHASERS:
            pts = []
            for caller in ("truth",):
                for depth in st.DEPTHS:
                    x = positive(value(db, (a, caller, platform, depth, phaser)))
                    y = positive(value(db, (b, caller, platform, depth, phaser)))
                    if finite(x) and finite(y):
                        pts.append((x, y))
            if pts:
                pairs[(platform, phaser)] = pts
    if not pairs:
        raise NoData("no run is computed on both %s and %s" % tuple(st.REFERENCES))
    return pairs


def identity_axes(ax, lo, hi, label):
    a, b = st.REFERENCES[0], st.REFERENCES[1]
    ax.plot([lo, hi], [lo, hi], color=st.MUTED, linewidth=0.6, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    log_ticks(ax.xaxis, lo, hi)
    log_ticks(ax.yaxis, lo, hi)
    ax.set_xlabel("%s, %s" % (label, st.REFERENCE_LABEL[a]))
    ax.set_ylabel("%s, %s" % (label, st.REFERENCE_LABEL[b]))
    ax.grid(color=st.GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)


def reference_scatter_file(db, ctx, name, value, label, cells, style_of, handles):
    pairs = reference_pairs(db, value)
    seen = [v for pts in pairs.values() for p in pts for v in p]
    lo, hi = min(seen) / 1.3, max(seen) * 1.3
    fig, ax = plt.subplots(figsize=(88 * st.MM, 96 * st.MM), constrained_layout=True)
    drawn = 0
    for platform, phaser in cells:
        pts = pairs.get((platform, phaser))
        if not pts:
            continue
        color, marker = style_of(platform, phaser)
        xs, ys = zip(*pts)
        ax.scatter(xs, ys, s=12, color=color, marker=marker, edgecolor="white",
                   linewidth=0.4, alpha=0.9, zorder=3)
        drawn += 1
    if not drawn:
        plt.close(fig)
        raise NoData("none of these runs is computed on both references")
    identity_axes(ax, lo, hi, label)
    # over the square plot, two or three to a row: one row of four is wider
    # than the plot and would stretch the figure around it
    legend_above(ax, handles, ncol=3 if len(handles) in (3, 6) else 2)
    return save(fig, ctx, name)


FINE_RATIO_TICKS = [0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04,
                    1.05, 1.06, 1.08, 1.1]


# fig6 ratios: a linear axis with 1 in the middle, as far below as above
REFERENCE_RATIO_LIMITS = "symmetric"
# half spans whose half is round too: the ticks are 1, 1 +- half/2, 1 +- half
NICE_HALF_SPANS = [0.02, 0.04, 0.06, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.6, 2]


def ratio_axis(ax, ratios, label, limits=None):
    """Log ratio, symmetric about 1 and zoomed to the data: a spread of 2% gets an
    axis of about +-3%, not the +-25% it used to get, so the boxes stay readable.
    The band of +-10% is shaded only when it fits inside the axis; on a zoomed
    axis every point is within 10% anyway. With limits the axis is linear
    instead: (lo, hi) as given, or "symmetric" - 1 +- the smallest round step
    that holds every ratio, the same distance below 1 as above."""
    a, b = st.REFERENCES[0], st.REFERENCES[1]
    if limits:
        if limits == "symmetric":
            far = max(abs(r - 1) for r in ratios) * 1.05
            half = next((h for h in NICE_HALF_SPANS if h >= far), NICE_HALF_SPANS[-1])
            limits = (1 - half, 1 + half)
        lo, hi = limits
        ax.set_ylim(lo, hi)
        if BAND < hi:
            ax.axhspan(1 / BAND, BAND, color=st.GRID, linewidth=0, zorder=0)
        ax.axhline(1, color=st.MUTED, linewidth=0.6, linestyle=(0, (4, 3)), zorder=1)
        ax.set_yticks([lo + i * (hi - lo) / 4 for i in range(5)])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
        ax.yaxis.set_minor_locator(NullLocator())
        # a ratio has no unit: "Switch error rate (%)" becomes "Switch error rate"
        ax.set_ylabel("%s,\n%s / %s" % (label.replace(" (%)", "").replace(" (Mb)", ""),
                                        st.REFERENCE_LABEL[b], st.REFERENCE_LABEL[a]))
        return
    span = max(abs(math.log10(r)) for r in ratios)
    span = max(span * 1.25, math.log10(1.03))
    lo, hi = 10 ** -span, 10 ** span
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    if BAND < hi * 0.97:
        ax.axhspan(1 / BAND, BAND, color=st.GRID, linewidth=0, zorder=0)
    ax.axhline(1, color=st.MUTED, linewidth=0.6, linestyle=(0, (4, 3)), zorder=1)
    ticks = [t for t in RATIO_TICKS if lo <= t <= hi]
    if len(ticks) < 4:
        ticks = [t for t in FINE_RATIO_TICKS if lo <= t <= hi]
        while len(ticks) > 7:
            ticks = ticks[::2] if 1 in ticks[::2] else ticks[1::2]
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_ylabel("%s,\n%s / %s" % (label, st.REFERENCE_LABEL[b], st.REFERENCE_LABEL[a]))


def reference_box_file(db, ctx, name, value, label, groups, title=None, limits=None):
    """groups: [(group label, [(box label, colour, [(platform, phaser), ...])])].
    Every point is one run on both references (same platform, coverage and
    phaser); a box pools the coverages. The box is the middle half of the
    ratios, its whiskers reach the furthest ratio within 1.5 box lengths, the
    black bar is the median, written as the change on T2T in percent just
    over the highest point of the box."""
    pairs = reference_pairs(db, value)
    all_ratios = [y / x for pts in pairs.values() for x, y in pts]
    def change(m):
        """The median ratio as a change on T2T: 0.91 -> -9%, 1.013 -> +1.3%. A
        percent cannot be read as a coverage, as 0.91x could."""
        pct = 100 * (m - 1)
        text = ("%+.1f%%" if abs(pct) < 10 else "%+.0f%%") % pct
        return text.replace("-", "\u2212")
    boxes = sum(len(boxes) for _, boxes in groups)
    fig, ax = plt.subplots(figsize=((34 + 11 * boxes) * st.MM, 70 * st.MM),
                           constrained_layout=True)
    rng = np.random.RandomState(0)
    pos, ticks, tick_labels = 0, [], []
    for g, (group_label, members) in enumerate(groups):
        first = pos
        for box_label, color, cells in members:
            r = [y / x for c in cells for x, y in pairs.get(c, [])]
            if r:
                ax.boxplot([r], positions=[pos], widths=0.6, patch_artist=True,
                           showfliers=False, zorder=3,
                           boxprops=dict(facecolor=color, alpha=0.3, edgecolor=color, linewidth=0.8),
                           medianprops=dict(color=st.INK, linewidth=1.2),
                           whiskerprops=dict(color=color, linewidth=0.8),
                           capprops=dict(color=color, linewidth=0.8))
                ax.scatter(pos + rng.uniform(-0.2, 0.2, len(r)), r, s=5, color=color,
                           edgecolor="white", linewidth=0.3, alpha=0.8, zorder=4)
                # over the highest point, so no point sits under the number
                top = max(r)
                if isinstance(limits, tuple):
                    top = min(top, limits[1])
                ax.annotate(change(float(np.median(r))), xy=(pos, top), xytext=(0, 2),
                            textcoords="offset points", ha="center", va="bottom", annotation_clip=False)
            ticks.append(pos)
            tick_labels.append(box_label)
            pos += 1
        if group_label:
            ax.annotate(group_label, xy=((first + pos - 1) / 2, 1),
                        xycoords=ax.get_xaxis_transform(), xytext=(0, 13),
                        textcoords="offset points", ha="center", va="bottom")
            if g < len(groups) - 1:
                ax.axvline(pos - 0.5 + 0.35, color=st.GRID, linewidth=0.8, zorder=0)
        pos += 0.7 if g < len(groups) - 1 else 0
    ratio_axis(ax, all_ratios, label, limits)
    if title:
        fig.suptitle(textwrap.fill(title, max(30, int((34 + 11 * boxes) / 1.9))))
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.6, pos - 0.4)
    ax.tick_params(axis="x", length=0)
    return save(fig, ctx, name)


def platform_style(platform, phaser):
    return st.PLATFORM_COLOR[platform], st.PLATFORM_MARKER[platform]


def phaser_shapes():
    return [Line2D([], [], color="none", marker=st.PHASER_MARKER[p], markerfacecolor=st.MUTED,
                   markeredgecolor="white", label=st.PHASER_LABEL[p]) for p in st.PHASERS]


def reference_files(db, ctx, prefix, value, label, title=None, limits=REFERENCE_RATIO_LIMITS):
    """Every file of one metric; a file with nothing to show is skipped and
    the others are still drawn. title and limits go to the ratio files."""
    every = [(pl, ph) for pl in st.PLATFORMS for ph in st.PHASERS]
    jobs = [
        ("scatter_all", lambda: reference_scatter_file(
            db, ctx, prefix + "_scatter_all", value, label, every,
            lambda pl, ph: (st.PHASER_COLOR[ph], st.PHASER_MARKER[ph]),
            phaser_handles(line=False))),
        ("scatter_by_platform", lambda: reference_scatter_file(
            db, ctx, prefix + "_scatter_by_platform", value, label, every,
            platform_style, platform_handles(line=False))),
    ]
    for tech, platforms in TECHNOLOGIES:
        jobs.append(("scatter_" + tech, lambda tech=tech, platforms=platforms: reference_scatter_file(
            db, ctx, prefix + "_scatter_" + tech, value, label,
            [(pl, ph) for pl in platforms for ph in st.PHASERS],
            lambda pl, ph: (st.PLATFORM_COLOR[pl], st.PHASER_MARKER[ph]),
            [h for h, pl in zip(platform_handles(line=False), st.PLATFORMS) if pl in platforms]
            + phaser_shapes())))
    by_phaser = [(st.PHASER_LABEL[p], st.PHASER_COLOR[p], [(pl, p) for pl in st.PLATFORMS])
                 for p in st.PHASERS]
    jobs += [
        ("ratio_by_phaser", lambda: reference_box_file(
            db, ctx, prefix + "_ratio_by_phaser", value, label, [("", by_phaser)], title, limits)),
        ("ratio_by_technology", lambda: reference_box_file(
            db, ctx, prefix + "_ratio_by_technology", value, label,
            [(tech, [(st.PHASER_LABEL[p], st.PHASER_COLOR[p], [(pl, p) for pl in platforms])
                     for p in st.PHASERS])
             for tech, platforms in TECHNOLOGIES], title, limits)),
        ("ratio_by_platform", lambda: reference_box_file(
            db, ctx, prefix + "_ratio_by_platform", value, label,
            [(tech, [(st.PLATFORM_LABEL[pl], st.PLATFORM_COLOR[pl], [(pl, p) for p in st.PHASERS])
                     for pl in platforms])
             for tech, platforms in TECHNOLOGIES], title, limits)),
    ]
    written = []
    for _, job in jobs:
        try:
            written.append(job())
        except NoData:
            pass
        finally:
            plt.close("all")
    if not written:
        raise NoData("no run is computed on both %s and %s" % tuple(st.REFERENCES))
    return ", ".join(written)


def fig6a(db, ctx):
    return reference_files(db, ctx, "fig6a_reference_switch_rate", switch_pct,
                           "Switch error rate (%)",
                           "Comparison of switch error rates between T2T-CHM13 and GRCh38 "
                           "across sequencing platforms")


def fig6b(db, ctx):
    return reference_files(db, ctx, "fig6b_reference_block_ng50", n50_mb,
                           "Phase block NG50 (Mb)",
                           "Comparison of phase block NG50 between T2T-CHM13 and GRCh38 "
                           "across sequencing platforms")


def fig6c(db, ctx):
    return reference_scatter(db, ctx, "fig6c_reference_recall_phased",
                             recall_pct, "Correctly phased truth het SNPs (%)", log=False)


def matpat_recall(db):
    """read recall against the Q100 parental-copy labels, genome-wide, by
    (reference, caller, platform, depth, phaser): metrics_reads/
    summary_recall_strata_matpat.tsv of READ_TRUTH=matpat bin/read_recall_strata.sh,
    as figR18 ..._matpat. The same reads on both references - primary on
    GRCh38 and on T2T, a parental label, 2 truth het SNPs or more - so the
    denominator is the same for both; recall = correct / truth_tagged."""
    path = os.path.join(db.sample_root, "metrics_reads", "summary_recall_strata_matpat.tsv")
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r.get("stratum") != "genome":
                continue
            try:
                tagged, correct = float(r["truth_tagged"]), float(r["correct"])
                key = (r["reference"], r["caller"], r["platform"], int(float(r["depth"])), r["phaser"])
            except (KeyError, ValueError):
                continue
            if tagged > 0:
                out[key] = 100 * correct / tagged
    return out


def fig6d(db, ctx):
    recall = matpat_recall(db)
    if not recall:
        raise NoData("no metrics_reads/summary_recall_strata_matpat.tsv - run "
                     "READ_TRUTH=matpat bash bin/read_recall_strata.sh on both references")
    return reference_files(db, ctx, "fig6d_reference_read_recall_matpat",
                           lambda db, key: recall.get(key, NAN),
                           "Haplotagging read recall (%)",
                           "Comparison of haplotagging recall between T2T-CHM13 and GRCh38 "
                           "across sequencing platforms")


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
    ("fig5c", fig5c, True),
    ("figS2a", figS2a, True),
    ("figS2b", figS2b, True),
    ("figS2c", figS2c, True),
    ("fig5d", fig5d, True),
    ("fig5e", fig5e, True),
    ("fig6a", fig6a, False),
    ("fig6b", fig6b, False),
    ("fig6c", fig6c, False),
    ("fig6d", fig6d, False),
    ("fig7a", fig7a, True),
    ("fig7b", fig7b, True),
    ("fig7c", fig7c, True),
    ("fig8", fig8, True),
    ("figS3a", figS3a, True),
    ("figS3b", figS3b, True),
    ("figS3g", figS3g, True),
    ("figS3c", figS3c, True),
    ("figS3d", figS3d, True),
    ("figS3e", figS3e, True),
    ("figS3f", figS3f, True),
    ("figS1", figS1, True),
    ("figS4", figS4, True),
]
