"""Figures of the read-level and difficult-region analyses, and of the other
GIAB samples. Same rules as figures.py: one PDF per subfigure, a figure with
nothing to show raises NoData.

  fig7d   read precision, blocks oriented by an anchor read against by SNPs
  fig7e   truth-labelled reads not tagged to their haplotype (1 - recall)
          against three read truths (whatshap on the truth VCF, alignment to
          the parental Q100 copies, parental k-mers), blocks oriented by their
          SNPs, reads labelled by both the truth VCF and that truth; bars
  fig7f   tagged reads on the wrong haplotype (1 - precision), the same truths
  fig7g   truth-tagged reads missed (1 - recall): anchor over SNP orientation
          against the SNP one, every run (Bland-Altman / MA plot)
  fig10a  why reads are not tagged: share of reads per class
  fig10b  the same in bases
  fig10c  confidence of the tagged reads
  fig10d  accuracy of the tagged reads by confidence
  fig10e  what the reads of each class look like (one file per platform)
  fig11a  where the switch errors gather: enrichment per stratum
  fig11b  the same for phase block ends
  fig11c  why a stratum is hard: allele balance, paralog-like links, depth,
          low MAPQ, read errors (one file per platform)
  fig12a-d T2T-CHM13 against GRCh38 on the whole genome and in each stratum:
          switch rate, phase block ends, reads tagged, reads on the wrong
          haplotype (the ratio of the two)
  fig9a-d the GIAB samples side by side at 30x, bars: switch rate, NG50,
          truth-tagged reads missed, FNR
"""

import math

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

import style as st
from bars import bar, group_ticks, label_bar, offsets, scale_axes
from data import NAN, ratio
from figures import (NoData, colorbar, empty_panel, finite, heatmap, log_norm, log_ticks, panel_row,
                     panel_size, paper_title,
                     platform_groups, platform_phaser_rows,
                     phaser_axis, phaser_handles, platform_handles, positive, save,
                     short_number, strata_heatmap, stratum_label, top_legend, y_grid)

CLASS_LABEL = {
    "tagged": "Tagged",
    "lowmapq": "MAPQ < 20",
    "no_snp_outside_bed": "Outside confident regions",
    "no_snp": "No het SNP on the read",
    "snps_unphased": "No usable phased SNP",
    "phased_unassigned": "SNPs phased, read not assigned",
}
CLASS_ORDER = list(CLASS_LABEL)
CLASS_COLOR = {
    "tagged": "#2C3A4B",
    "lowmapq": "#E39B7B",
    "no_snp_outside_bed": "#D9D9D9",
    "no_snp": "#A9B8C6",
    "snps_unphased": "#6D7F92",
    "phased_unassigned": "#E3A034",
}
# read_classes.py splits MAPQ < 20 by segdups / low mappability cut to the
# confident regions, where few low-MAPQ reads are; the figures show one class
CLASS_PARTS = {"lowmapq": ("lowmapq_paralog", "lowmapq_other")}
SUMMED = ("reads", "bases", "compared", "correct")
# the weights of the per-read means and per-base rates when parts are merged;
# vaf_dev and paralog_link_frac are averaged over reads, not over their SNPs
WEIGHTED = {"mean_length": "reads", "mean_mapq": "reads", "mean_bq": "reads",
            "vaf_dev": "reads", "paralog_link_frac": "reads", "softclip_frac": "bases",
            "nm_rate": "bases", "mismatch_rate": "bases", "indel_rate": "bases"}


def group_row(table, stratum, group_type, group):
    """One row of a read_classes.py table; a merged class sums the counts of
    its parts and averages their features, weighted as in WEIGHTED."""
    parts = CLASS_PARTS.get(group) if group_type == "class" else None
    if not parts:
        return table.get((stratum, group_type, group))
    rows = [table[(stratum, group_type, p)] for p in parts if (stratum, group_type, p) in table]
    if not rows:
        return None
    out = {c: sum(r.get(c) or 0 for r in rows) for c in SUMMED}
    for col, weight in WEIGHTED.items():
        pairs = [(r[col], r.get(weight) or 0) for r in rows if finite(r.get(col))]
        total = sum(w for _, w in pairs)
        out[col] = sum(v * w for v, w in pairs) / total if total else NAN
    return out
CONF_LABEL = {"high": "High", "medium": "Medium", "single": "One SNP", "tie": "Tie", "none": "No SNP"}
CONF_ORDER = list(CONF_LABEL)
CONF_COLOR = {"high": "#2C3A4B", "medium": "#6D7F92", "single": "#B9C4CF", "tie": "#E3A034",
              "none": "#E6E6E6"}

FEATURES = [
    ("mean_length", "Read length (kb)", lambda v: v / 1000),
    ("mean_mapq", "MAPQ", lambda v: v),
    ("mismatch_rate", "Mismatches (%)", lambda v: 100 * v),
    ("indel_rate", "Indels (%)", lambda v: 100 * v),
    ("softclip_frac", "Soft-clipped (%)", lambda v: 100 * v),
    ("mean_bq", "Base quality", lambda v: v),
    ("vaf_dev", "|VAF - 0.5|", lambda v: v),
    ("paralog_link_frac", "Links > 2 patterns (%)", lambda v: 100 * v),
]

DIVERGING = LinearSegmentedColormap.from_list("enrich", ["#2E6F9E", "#F7F7F7", "#B5523B"])


def platform_bars(ax, groups, order, colors, value):
    """Stacked 100% bars, one per phaser."""
    for i, phaser in enumerate(st.PHASERS):
        shares = groups.get(phaser)
        if not shares:
            continue
        bottom = 0.0
        for g in order:
            v = shares.get(g, 0.0)
            ax.bar(i, v, bottom=bottom, width=0.7, color=colors[g], edgecolor="white",
                   linewidth=0.4, zorder=3)
            bottom += v


def class_shares(db, ref, depth, caller, column, group_type, order, tagged_only=False):
    """platform -> phaser -> group -> % of the reads (or bases) of the run on chr1-22"""
    out = {}
    for platform in st.PLATFORMS:
        for phaser in st.PHASERS:
            table = db.classes((ref, caller, platform, depth, phaser))
            if not table:
                continue
            gtype = "tagged_confidence" if tagged_only else group_type
            counts = {g: (group_row(table, "genome", gtype, g) or {}).get(column, 0) or 0 for g in order}
            total = sum(counts.values())
            if not total:
                continue
            out.setdefault(platform, {})[phaser] = {g: 100 * counts[g] / total for g in order}
    return out


def stacked_figure(db, ctx, ref, name, column, group_type, order, labels, colors, ylabel,
                   tagged_only=False):
    depth = ctx.depth
    shares = class_shares(db, ref, depth, "truth", column, group_type, order, tagged_only)
    if not shares:
        raise NoData("no metrics_reads classes for %s at %dx - run bin/read_classes.sh" % (ref, depth))
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), height=st.PANEL_HEIGHT + 16)
    for ax, platform in zip(axes, st.PLATFORMS):
        if platform in shares:
            platform_bars(ax, shares[platform], order, colors, column)
        else:
            empty_panel(ax)
        ax.set_title(st.PLATFORM_LABEL[platform])
        phaser_axis(ax)
        y_grid(ax)
    axes[0].set_ylim(0, 100)
    axes[0].set_ylabel(ylabel + ", %d×" % depth)
    top_legend(legend_ax, [Patch(facecolor=colors[g], label=labels[g]) for g in order], ncol=3)
    return save(fig, ctx, "%s_%dx_%s" % (name, depth, ref))


CLASS_TICK = {
    "tagged": "Tagged",
    "lowmapq": "MAPQ < 20",
    "no_snp_outside_bed": "Outside\nconfident\nregions",
    "no_snp": "No het SNP\non the read",
    "snps_unphased": "No usable\nphased SNP",
    "phased_unassigned": "SNPs phased,\nread not\nassigned",
}
UNTAGGED = CLASS_ORDER[1:]


def pct_text(v):
    return "%.1f" % v if v >= 1 else "%.2f" % v


def grouped_bars(ax, groups, shares, width=0.2):
    """A group per read class, a bar per platform, each value written above
    its bar so a 0.1% class is read as easily as a 70% one."""
    pad = 0.015 * ax.get_ylim()[1]
    for i, g in enumerate(groups):
        for j, platform in enumerate(st.PLATFORMS):
            v = shares.get(platform, {}).get(g)
            if v is None:
                continue
            x = i + (j - 1.5) * width
            ax.bar(x, v, width=width * 0.9, color=st.PLATFORM_COLOR[platform], zorder=3)
            ax.text(x, v + pad, pct_text(v), rotation=90, ha="center", va="bottom",
                    fontsize=st.FONT_SIZE - 1, color=st.INK, zorder=4)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([CLASS_TICK[g] for g in groups])
    ax.set_xlim(-0.55, len(groups) - 0.45)
    ax.tick_params(axis="x", length=0)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    y_grid(ax)


def class_bars(db, ctx, ref, name, column, unit, ylabel, phaser="whatshap"):
    """Tagged on its own scale, the reasons a read is not tagged on another,
    a bar per platform; one phaser (the classes beside 'tagged' barely depend
    on it, the tables have all four)."""
    depth = ctx.depth
    shares = class_shares(db, ref, depth, "truth", column, "class", CLASS_ORDER)
    shares = {p: v[phaser] for p, v in shares.items() if phaser in v}
    if not shares:
        raise NoData("no metrics_reads classes for %s at %dx - run bin/read_classes.sh" % (ref, depth))
    fig = plt.figure(figsize=(st.FULL_WIDTH * st.MM, 72 * st.MM), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[0.08, 1], width_ratios=[1.4, len(UNTAGGED)])
    legend_ax = fig.add_subplot(grid[0, :])
    legend_ax.axis("off")
    left, right = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    rest = [v for p in shares.values() for g, v in p.items() if g != "tagged"]
    left.set_ylim(0, 100)
    right.set_ylim(0, max(rest) * 1.3 if rest else 1)
    grouped_bars(left, ["tagged"], shares)
    grouped_bars(right, UNTAGGED, shares)
    left.set_ylabel(ylabel + ", %d×" % depth)
    left.set_title("%s tagged" % unit)
    right.set_title("%s not tagged, by reason (%s)" % (unit, st.PHASER_LABEL[phaser]))
    top_legend(legend_ax, [Patch(facecolor=st.PLATFORM_COLOR[p], label=st.PLATFORM_LABEL[p])
                           for p in st.PLATFORMS])
    return save(fig, ctx, "%s_%dx_%s" % (name, depth, ref))


def fig10a(db, ctx, ref):
    return class_bars(db, ctx, ref, "fig10a_read_classes_reads", "reads", "Reads",
                      "Reads on chr1-22 (%)")


def fig10b(db, ctx, ref):
    return class_bars(db, ctx, ref, "fig10b_read_classes_bases", "bases", "Read bases",
                      "Read bases on chr1-22 (%)")


def fig10c(db, ctx, ref):
    return stacked_figure(db, ctx, ref, "fig10c_tagged_confidence", "reads", "confidence",
                          CONF_ORDER, CONF_LABEL, CONF_COLOR, "Tagged reads (%)", tagged_only=True)


CONF_TICK = {
    "high": "High\n(2+ SNPs,\nall agree)",
    "medium": "Medium\n(majority,\nsome against)",
    "single": "One SNP",
    "tie": "Tie\n(as many\nfor each)",
}


def fig10d(db, ctx, ref):
    """Tagged reads on the wrong haplotype, by how the read's SNPs vote for
    its haplotype, blocks oriented by their SNPs; platform x phaser down, one
    confidence class per column. Fewer than 20 compared reads: hatched."""
    depth = ctx.depth
    confs = CONF_ORDER[:4]
    rows = platform_phaser_rows()

    def wrong(platform, phaser, conf):
        row = db.classes((ref, "truth", platform, depth, phaser)).get(("genome", "tagged_confidence", conf))
        if not row or not row.get("compared") or row["compared"] < 20:
            return NAN
        return 100 * (1 - row["correct"] / row["compared"])

    matrix = [[wrong(pl, p, c) for c in confs] for pl, p in rows]
    shown = [v for line in matrix for v in line if finite(v) and v > 0]
    if not shown:
        raise NoData("no accuracy by confidence for %s at %dx" % (ref, depth))
    floor = min(shown) / 2
    matrix = [[floor if finite(v) and v <= 0 else v for v in line] for line in matrix]
    fig, ax = plt.subplots(figsize=((62 + 18 * len(confs)) * st.MM, 100 * st.MM),
                           constrained_layout=True)
    im = heatmap(ax, matrix, [st.PHASER_LABEL[p] for _, p in rows], [CONF_TICK[c] for c in confs],
                 platform_groups(rows),
                 [("", len(confs))], label=lambda v: "0" if v <= floor else short_number(v),
                 cmap=st.SEQUENTIAL, norm=log_norm(shown + [floor]))
    colorbar(fig, im, ax, "Tagged reads on the wrong haplotype (%%), %d×" % depth)
    return save(fig, ctx, "fig10d_accuracy_by_confidence_%dx_%s" % (depth, ref))


def feature_heatmap(ax, fig, rows, row_labels, values, title, features=None):
    """rows x features; colour scaled per column, raw values written in."""
    features = features or FEATURES
    matrix = np.full((len(rows), len(features)), np.nan)
    text = [["" for _ in features] for _ in rows]
    for i, r in enumerate(rows):
        for j, (col, _, f) in enumerate(features):
            v = values.get(r, {}).get(col, NAN)
            if finite(v):
                matrix[i, j] = f(v)
                text[i][j] = short_number(abs(f(v))) if f(v) >= 0 else "-" + short_number(-f(v))
    scaled = np.full_like(matrix, np.nan)
    for j in range(matrix.shape[1]):
        col = matrix[:, j]
        ok = np.isfinite(col)
        if ok.sum() == 0:
            continue
        lo, hi = np.nanmin(col), np.nanmax(col)
        scaled[ok, j] = 0.5 if hi == lo else (col[ok] - lo) / (hi - lo)
    im = heatmap(ax, scaled, row_labels, [label for _, label, _ in features],
                 [("", len(rows))], [("", len(features))], cmap=st.SEQUENTIAL, vmin=0, vmax=1)
    for i in range(len(rows)):
        for j in range(len(features)):
            if text[i][j]:
                ax.text(j, i, text[i][j], ha="center", va="center", zorder=3,
                        color="white" if scaled[i, j] > 0.55 else st.INK)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
    ax.set_title(title)
    return im


def fig10e(db, ctx, ref, phaser="whatshap"):
    """Per platform: the read classes against their features, truth input,
    one phaser (the classes beside 'tagged' barely depend on it)."""
    depth = ctx.depth
    written = []
    for platform in st.PLATFORMS:
        table = db.classes((ref, "truth", platform, depth, phaser))
        values = {c: group_row(table, "genome", "class", c) for c in CLASS_ORDER}
        rows = [c for c in CLASS_ORDER if values[c]]
        if not rows:
            continue
        fig, ax = plt.subplots(figsize=(st.FULL_WIDTH * 0.75 * st.MM, (30 + 7 * len(rows)) * st.MM),
                               constrained_layout=True)
        feature_heatmap(ax, fig, rows, [CLASS_LABEL[c] for c in rows], values,
                        "%s, %s, %d×" % (st.PLATFORM_LABEL[platform], st.PHASER_LABEL[phaser], depth))
        written.append(save(fig, ctx, "fig10e_class_features_%s_%dx_%s" % (platform, depth, ref)))
    if not written:
        raise NoData("no read classes for %s at %dx" % (ref, depth))
    return ", ".join(written)


def enrichment_heatmap(ctx, ref, name, strata, value, label, note=None):
    depth = ctx.depth
    values = [value(pl, p, s) for pl in st.PLATFORMS for p in st.PHASERS for s in strata]
    values = [v for v in values if finite(v) and v > 0]
    if not values:
        raise NoData("no %s for %s at %dx" % (name, ref, depth))
    # symmetric about 1, so depletion and enrichment get the same colour depth
    span = max(abs(math.log10(v)) for v in values) or 1
    return strata_heatmap(ctx, "%s_%dx_%s" % (name, depth, ref), strata, value,
                          "%s, %d×" % (label, depth), note=note, cmap=DIVERGING,
                          norm=LogNorm(vmin=10 ** -span, vmax=10 ** span))


def switch_enrichment(db, key, stratum):
    """Switch rate among the SNPs of the stratum (regions.sh, the numbers of
    fig8) over the run's switch rate on chr1-22: a flip counts in the stratum
    of the flipped SNP, and a mis-phased stretch in the stratum it lies in."""
    row = db.stratum(key, stratum)
    genome = db.get(key, "switch_rate") if db.has(key) else NAN
    if not row or not finite(genome) or genome <= 0:
        return NAN
    return ratio(row["switches"], row["assessed"]) / genome


def zero_note(row):
    """'0' for a stratum row without a switch, else None"""
    return "0" if row and row["assessed"] and not row["switches"] else None


def fig11a(db, ctx, ref):
    depth = ctx.depth
    return enrichment_heatmap(
        ctx, ref, "fig11a_switch_enrichment", db.strata_names(),
        lambda pl, p, s: switch_enrichment(db, (ref, "truth", pl, depth, p), s),
        "Switch error rate in the stratum,\nrelative to the genome",
        note=lambda pl, p, s: zero_note(db.stratum((ref, "truth", pl, depth, p), s)))


def fig11b(db, ctx, ref):
    depth = ctx.depth
    mech = db.mechanism(ref)
    if not mech:
        raise NoData("no mechanism_%s.tsv - run bin/strata_mechanism.sh" % ref)
    strata = [s for s in db.strata_names()
              if any(k[4] == s and k[5] == "break_enrichment" for k in mech)]
    return enrichment_heatmap(
        ctx, ref, "fig11b_block_end_enrichment", strata,
        lambda pl, p, s: mech.get((pl, depth, "truth", p, s, "break_enrichment"), NAN),
        "Phase block ends per het SNP,\nrelative to the genome")


STRATUM_FEATURES = [
    ("vaf_dev", "|VAF - 0.5|", lambda v: v),
    ("paralog_link_frac", "Links > 2 patterns (%)", lambda v: 100 * v),
    ("lowmapq_frac", "Depth at MAPQ < 20 (%)", lambda v: 100 * v),
    ("mismatch_rate", "Read mismatches (%)", lambda v: 100 * v),
    ("indel_rate", "Read indels (%)", lambda v: 100 * v),
    ("tagged_reads", "Reads tagged (%)", lambda v: 100 * v),
]


def fig11c(db, ctx, ref, phaser="whatshap"):
    """Per platform: each stratum against what makes it hard, truth input."""
    depth = ctx.depth
    mech = db.mechanism(ref)
    written = []
    for platform in st.PLATFORMS:
        table = db.classes((ref, "truth", platform, depth, phaser))
        strata = [s for s in db.strata_names()
                  if (platform, depth, "-", "-", s, "truth_snps") in mech or ("%s" % s, "all", "all") in table]
        if not strata:
            continue
        values = {}
        for s in strata:
            v = {k: mech.get((platform, depth, "-", "-", s, k), NAN)
                 for k in ("vaf_dev", "paralog_link_frac", "lowmapq_frac")}
            allr = table.get((s, "all", "all"), {})
            tagged = table.get((s, "class", "tagged"), {})
            v["mismatch_rate"] = allr.get("mismatch_rate", NAN)
            v["indel_rate"] = allr.get("indel_rate", NAN)
            v["tagged_reads"] = (tagged.get("reads", 0) / allr["reads"]) if allr.get("reads") else NAN
            values[s] = v
        fig, ax = plt.subplots(figsize=(st.FULL_WIDTH * 0.75 * st.MM, (30 + 7 * len(strata)) * st.MM),
                               constrained_layout=True)
        feature_heatmap(ax, fig, strata, [stratum_label(s) for s in strata], values,
                        "%s, %d×" % (st.PLATFORM_LABEL[platform], depth), features=STRATUM_FEATURES)
        written.append(save(fig, ctx, "fig11c_stratum_features_%s_%dx_%s" % (platform, depth, ref)))
    if not written:
        raise NoData("no stratum features for %s at %dx" % (ref, depth))
    return ", ".join(written)


def fig7d(db, ctx, ref):
    """Every run, reads on the wrong haplotype (1 - read precision)
    with blocks oriented by an anchor read against by their SNPs (both with
    blocks of two SNPs or more, chr1-22). Colour is the platform, shape the
    phaser. Each platform's legend entry carries the median over its runs of
    anchor / SNP: how many times more reads the anchor read puts on the wrong
    haplotype. Every run sits over the diagonal, so the platforms go into the
    empty corner under it and the phasers under the plot."""
    fig = plt.figure(figsize=panel_size("fig1d"), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1, 0.13])
    ax = fig.add_subplot(grid[0])
    legend_ax = fig.add_subplot(grid[1])
    legend_ax.axis("off")
    seen, ratios, phasers = [], {}, []
    for platform in st.PLATFORMS:
        for phaser in st.PHASERS:
            xs, ys = [], []
            for caller in st.CALLERS:
                for depth in st.DEPTHS:
                    key = (ref, caller, platform, depth, phaser)
                    x = 100 * (1 - db.get(key, "read_precision"))
                    y = 100 * (1 - db.get(key, "read_precision_anchor_truthtags"))
                    if finite(x) and finite(y):
                        xs.append(max(x, 1e-3))
                        ys.append(max(y, 1e-3))
            if xs:
                ax.scatter(xs, ys, s=8, color=st.PLATFORM_COLOR[platform], marker=st.PHASER_MARKER[phaser],
                           edgecolor="white", linewidth=0.3, alpha=0.9, zorder=3)
                seen += xs + ys
                ratios.setdefault(platform, []).extend(y / x for x, y in zip(xs, ys))
                if phaser not in phasers:
                    phasers.append(phaser)
    if not seen:
        plt.close(fig)
        raise NoData("no anchor-oriented accuracy for %s - run readacc_anchor.sh" % ref)
    lo, hi = min(seen) / 1.5, max(seen) * 1.5
    # x reaches further right and y further down than the runs, so the corner
    # under the diagonal is free for the platform legend
    ax.plot([lo / 6, hi * 6], [lo / 6, hi * 6], color=st.MUTED, linewidth=0.6, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi * 6)
    ax.set_ylim(lo / 6, hi)
    log_ticks(ax.xaxis, lo, hi * 6)
    log_ticks(ax.yaxis, lo / 6, hi)
    ax.set_xlabel("Reads on the wrong haplotype (%),\nblocks oriented by their SNPs")
    ax.set_ylabel("Reads on the wrong haplotype (%),\nblocks oriented by an anchor read")
    ax.grid(color=st.GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(handles=[Line2D([], [], color="none", marker="o", markerfacecolor=st.PLATFORM_COLOR[p],
                              markeredgecolor="white",
                              label="%s  ×%.2f" % (st.PLATFORM_LABEL[p], float(np.median(ratios[p]))))
                       for p in st.PLATFORMS if p in ratios],
              title="Anchor / SNP (median)", loc="lower right", borderaxespad=0.2,
              handletextpad=0.1, labelspacing=0.25)
    legend_ax.legend(handles=[Line2D([], [], color="none", marker=st.PHASER_MARKER[p], markerfacecolor=st.MUTED,
                                     markeredgecolor="white", label=st.PHASER_LABEL[p]) for p in phasers],
                     loc="center", ncol=3, borderaxespad=0, handletextpad=0.2, columnspacing=0.8)
    paper_title(fig, "Comparison of haplotagging error rates using different phase-block "
                     "orientation strategies")
    return save(fig, ctx, "fig7d_anchor_vs_variant_orientation_" + ref)


def fig7g(db, ctx, ref):
    """Every run as in a Bland-Altman / MA plot: x the share of truth-tagged
    reads not assigned to their haplotype (1 - read recall) with blocks
    oriented by their SNPs, y how many times more there are with blocks
    oriented by an anchor read. No difference is the dashed line at 1; over it
    the SNP orientation is better. y runs from 0.8 to 2, so the gap from 1
    shows; runs outside are left out. Fig. 1d, as large as fig6 by platform
    beside it; no legend of its own: colour (platform) and shape (phaser) are
    those of Fig. 1c, whose legend serves both."""
    fig, ax = plt.subplots(figsize=panel_size("small"), constrained_layout=True)
    xs_all, ys_all = [], []
    for platform in st.PLATFORMS:
        for phaser in st.PHASERS:
            if not st.applies(platform, phaser):
                continue
            xs, ys = [], []
            for caller in st.CALLERS:
                for depth in st.DEPTHS:
                    key = (ref, caller, platform, depth, phaser)
                    snp = 1 - db.get(key, "read_recall")
                    anchor = 1 - db.get(key, "read_recall_anchor_truthtags")
                    if finite(snp) and finite(anchor) and snp > 0 and anchor > 0:
                        xs.append(100 * snp)
                        ys.append(anchor / snp)
            ax.scatter(xs, ys, s=8, color=st.PLATFORM_COLOR[platform], marker=st.PHASER_MARKER[phaser],
                       alpha=0.8, edgecolor="white", linewidth=0.3, zorder=3)
            xs_all += xs
            ys_all += ys
    if not xs_all:
        plt.close(fig)
        raise NoData("no anchor-oriented read recall for %s - run readacc_anchor.sh" % ref)
    xlo, xhi = min(xs_all) / 1.2, max(xs_all) * 1.2
    ylo, yhi = 0.8, 2.0
    ax.set_xscale("log")
    ax.set_yscale("log", base=2)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    log_ticks(ax.xaxis, xlo, xhi)
    ax.yaxis.set_major_locator(FixedLocator([0.8, 1, 1.25, 1.5, 2]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.axhline(1, color=st.MUTED, linewidth=0.6, linestyle=(0, (4, 3)), zorder=2)
    paper_title(fig, "Comparison of haplotagging error rates using different phase-block "
                     "orientation strategies")
    ax.set_xlabel("Reads not on their haplotype (%),\nSNP orientation")
    ax.set_ylabel("Fold increase,\nanchor read / SNP orientation")
    ax.grid(color=st.GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    return save(fig, ctx, "fig7g_anchor_vs_variant_recall_ratio_" + ref)


TRUTH_LABEL = {"truthtags": "GIAB truth VCF (whatshap haplotag)",
               "matpat": "Q100 parental copies (alignment)",
               "yak": "Parental k-mers (yak)"}
TRUTH_COLOR = {"truthtags": "#B9C4CF", "matpat": "#2C3A4B", "yak": "#6D7F92"}


def read_truth_bars(db, ctx, ref, column, ylabel, name):
    """Truth-input runs: a panel per platform, a group per phaser, a bar per
    read truth - which haplotype a read comes from, by the truth VCF, by
    alignment to the two parental copies of the Q100 assembly, or by parental
    k-mers (bin/readacc_truths.sh). Every block is oriented by its SNPs, the
    same for all three truths, and only reads with a label of both this truth
    and the truth VCF count, so the bars differ by the labels alone. The bars
    are the error, not the accuracy: 97% against 98% is invisible on bars, 3%
    against 2% is not."""
    depth = ctx.depth
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), height=st.PANEL_HEIGHT + 12)
    xs, w = offsets(len(TRUTH_LABEL))
    drawn = []
    for ax, platform in zip(axes, st.PLATFORMS):
        phasers = [p for p in st.PHASERS if st.applies(platform, p)]
        for x, phaser in enumerate(phasers):
            key = (ref, "truth", platform, depth, phaser)
            for dx, truth in zip(xs, TRUTH_LABEL):
                v = bar(ax, x + dx, 100 * (1 - db.get(key, "%s_common_var_%s" % (column, truth))),
                        w, TRUTH_COLOR[truth])
                if v is not None:
                    drawn.append(v)
                    label_bar(ax, x + dx, v, w)
        group_ticks(ax, [st.PHASER_LABEL[p] for p in phasers])
        ax.set_title(st.PLATFORM_LABEL[platform])
    if not drawn:
        plt.close(fig)
        raise NoData("no read accuracy against the read truths for %s at %dx" % (ref, depth))
    scale_axes(axes, drawn, False, ylabel)
    top_legend(legend_ax, [Patch(color=TRUTH_COLOR[t], label=label) for t, label in TRUTH_LABEL.items()])
    return save(fig, ctx, name)


def fig7e(db, ctx, ref):
    """Reads with a label that the run does not tag to that haplotype: 1 -
    read recall (untagged reads and reads on the wrong haplotype)."""
    return read_truth_bars(db, ctx, ref, "read_recall",
                           "Truth-labelled reads not tagged\nto their haplotype (%%), %d×" % ctx.depth,
                           "fig7e_read_missed_by_truth_%dx_%s" % (ctx.depth, ref))


def fig7f(db, ctx, ref):
    """Tagged reads on the wrong haplotype by each truth: 1 - read precision."""
    return read_truth_bars(db, ctx, ref, "read_precision",
                           "Tagged reads on the wrong\nhaplotype (%%), %d×" % ctx.depth,
                           "fig7f_read_wrong_by_truth_%dx_%s" % (ctx.depth, ref))


# ------------------------------------------------------------- other samples

# one neutral shade per sample, darkest first, so a sample never looks like a
# platform or a phaser
SAMPLE_SHADES = ["#2C3A4B", "#6D7F92", "#B9C4CF", "#8C8C8C", "#D4D4D4"]


def phaser_value(column, scale, complement=False):
    """value(db, ref, platform, depth, phaser) of the truth-input run"""
    def value(db, ref, platform, depth, phaser):
        v = db.get((ref, "truth", platform, depth, phaser), column)
        return scale * (1 - v) if complement else scale * v
    return value


def caller_value(column, scale):
    return lambda db, ref, platform, depth, caller: scale * db.call(ref, caller, platform, depth, column)


SAMPLE_METRICS = [
    ("fig9a_samples_switch_rate", phaser_value("switch_rate", 100), "Switch error rate (%)", True, "phaser"),
    ("fig9b_samples_block_ng50", phaser_value("block_ng50", 1e-6), "Phase block NG50 (Mb)", True, "phaser"),
    ("fig9c_samples_reads_missed", phaser_value("read_recall", 100, complement=True),
     "Truth-tagged reads not tagged\nto their haplotype (%)", False, "phaser"),
    ("fig9d_samples_caller_fnr", caller_value("fnr", 100), "Truth het SNPs missed (FNR, %)", True, "caller"),
]


def samples_figure(dbs, ctx, name, value, ylabel, log, by, depth=30, ref="GRCh38"):
    """The GIAB samples side by side at 30x on GRCh38: a panel per platform that
    two samples or more were sequenced on, a group per phaser on the truth
    input (or per caller), a bar per sample. Each sample is scored against its
    own truth (HG002 v5.0q, HG005 v4.2.1, HG001 Platinum Pedigree v1.2), so a
    difference is the sample and its truth together."""
    series_all = st.PHASERS if by == "phaser" else [c for c in st.CALLERS if c != "truth"]

    def val(sample, platform, s):
        v = value(dbs[sample], ref, platform, depth, s)
        return positive(v) if log else v

    has = lambda sample, platform: any(finite(val(sample, platform, s)) for s in series_all)
    samples = [s for s in dbs if any(has(s, p) for p in st.PLATFORMS)]
    platforms = [p for p in st.PLATFORMS if sum(has(s, p) for s in samples) >= 2]
    if not platforms:
        raise NoData("no platform with two samples at %dx" % depth)
    color = {s: SAMPLE_SHADES[i % len(SAMPLE_SHADES)] for i, s in enumerate(samples)}
    fig, legend_ax, axes = panel_row(len(platforms), height=st.PANEL_HEIGHT + 12)
    xs, w = offsets(len(samples))
    drawn = []
    for ax, platform in zip(axes, platforms):
        series = [s for s in series_all if by != "phaser" or st.applies(platform, s)]
        for x, s in enumerate(series):
            for dx, sample in zip(xs, samples):
                v = bar(ax, x + dx, val(sample, platform, s), w, color[sample])
                if v is not None:
                    drawn.append(v)
                    label_bar(ax, x + dx, v, w)
        group_ticks(ax, [st.PHASER_LABEL[s] if by == "phaser" else st.CALLER_LABEL[s] for s in series])
        ax.set_title(st.PLATFORM_LABEL[platform])
    if not drawn:
        plt.close(fig)
        raise NoData("no values for any sample at %dx" % depth)
    scale_axes(axes, drawn, log, "%s, %d×" % (ylabel, depth))
    top_legend(legend_ax, [Patch(color=color[s], label=s) for s in samples])
    return save(fig, ctx, "%s_%dx_%s" % (name, depth, ref))


def sample_figures(dbs, ctx):
    out = []
    for name, value, ylabel, log, by in SAMPLE_METRICS:
        try:
            out.append(samples_figure(dbs, ctx, name, value, ylabel, log, by))
        except NoData:
            pass
        finally:
            plt.close("all")
    if not out:
        raise NoData("no sample has results at 30x")
    return ", ".join(out)


# (id, function, drawn once per reference?)
# ------------------------------------- T2T against GRCh38 in the strata: 12
#
# Each metric(db, ref, platform, phaser, stratum, depth) -> (value, errors): the
# rate in the stratum ("genome" for the whole of chr1-22) and the number of
# errors behind it (None where the rate is not a count of errors); truth input.
# Under MIN_ERRORS errors on either reference the two rates are not compared.

MIN_ERRORS = 10


def stratum_switch(db, ref, platform, phaser, stratum, depth):
    """switch errors / assessed pairs among the SNPs of the stratum (regions.sh)"""
    key = (ref, "truth", platform, depth, phaser)
    if stratum == "genome":
        if not db.has(key):
            return NAN, NAN
        n, k = db.get(key, "assessed"), db.get(key, "switches")
    else:
        row = db.stratum(key, stratum)
        if not row:
            return NAN, NAN
        n, k = row["assessed"], row["switches"]
    if not (finite(n) and n > 0 and finite(k)):
        return NAN, NAN
    return 100 * k / n, k


def stratum_block_ends(db, ref, platform, phaser, stratum, depth):
    """phase block ends per 1,000 truth het SNPs of the stratum, a block end
    placed at the middle of the gap between two blocks (strata_mechanism.py);
    NG50 has no per-stratum form, this is what takes its place"""
    mech = db.mechanism(ref)
    snps_of = lambda s: mech.get((platform, depth, "-", "-", s, "truth_snps"), NAN)
    snps = (snps_of("notinalldifficultregions") + snps_of("alldifficultregions")
            if stratum == "genome" else snps_of(stratum))
    k = mech.get((platform, depth, "truth", phaser, stratum, "breaks"), NAN)
    if not (finite(snps) and snps > 0 and finite(k)):
        return NAN, NAN
    return 1000 * k / snps, k


def stratum_tagged(db, ref, platform, phaser, stratum, depth):
    """reads tagged in a block of 2+ SNPs, of the reads with half their
    alignment or more in the stratum (read_classes.py)"""
    table = db.classes((ref, "truth", platform, depth, phaser))
    every = table.get((stratum, "all", "all"))
    if not every or not every.get("reads"):
        return NAN, NAN
    tagged = table.get((stratum, "class", "tagged")) or {}
    return 100 * (tagged.get("reads") or 0) / every["reads"], None


def stratum_wrong(db, ref, platform, phaser, stratum, depth):
    """tagged reads on the wrong haplotype, of those with a truth label"""
    tagged = db.classes((ref, "truth", platform, depth, phaser)).get((stratum, "class", "tagged"))
    if not tagged or not tagged.get("compared"):
        return NAN, NAN
    n = tagged["compared"]
    return 100 * (n - tagged["correct"]) / n, n - tagged["correct"]


def comparable(g, g_errors, t, t_errors):
    """two rates are put side by side only on MIN_ERRORS errors or more on each
    reference (and a rate that is not an error count only when both are > 0)"""
    if g_errors is None or t_errors is None:
        return g > 0 and t > 0
    return g_errors >= MIN_ERRORS and t_errors >= MIN_ERRORS


def reference_strata(db, ctx, name, metric, label, higher_is_better=False):
    """T2T over GRCh38 for one rate: platform x phaser down, the whole genome
    and each stratum across. A stratum of one name is cut to each reference's
    own benchmark regions, so the two columns are the same kind of sequence,
    not the same bases. Blue: better on T2T. Fewer than MIN_ERRORS errors on
    either reference: dotted and '-', not comparable (both rates in t12g)."""
    depth = ctx.depth
    strata = ["genome"] + db.strata_names()

    def pair(pl, p, s):
        return (metric(db, "GRCh38", pl, p, s, depth), metric(db, "T2T", pl, p, s, depth))

    def value(pl, p, s):
        (g, _), (t, _) = pair(pl, p, s)
        return t / g if finite(g) and finite(t) and g > 0 and t > 0 else NAN

    def note(pl, p, s):
        (g, g_errors), (t, t_errors) = pair(pl, p, s)
        if not (finite(g) and finite(t)):
            return None
        return "–" if not comparable(g, g_errors, t, t_errors) else None

    values = [value(pl, p, s) for pl in st.PLATFORMS for p in st.PHASERS for s in strata]
    values = [v for v in values if finite(v) and v > 0]
    if not values:
        raise NoData("no %s on both references at %dx" % (name, depth))
    # symmetric about 1 and never narrower than 0.5-2, so a small difference
    # stays pale instead of filling the scale
    span = max(max(abs(math.log10(v)) for v in values), math.log10(2))
    out = strata_heatmap(ctx, "%s_%dx" % (name, depth), strata, value,
                         "%s,\nT2T-CHM13 / GRCh38, %d×" % (label, depth), note=note,
                         fmt=lambda v: "%.2f" % v,
                         cmap=DIVERGING.reversed() if higher_is_better else DIVERGING,
                         norm=LogNorm(vmin=10 ** -span, vmax=10 ** span))
    if out is None:
        raise NoData("no %s on both references at %dx" % (name, depth))
    return out


def fig12a(db, ctx):
    return reference_strata(db, ctx, "fig12a_reference_strata_switch_rate", stratum_switch,
                            "Switch error rate")


def fig12b(db, ctx):
    return reference_strata(db, ctx, "fig12b_reference_strata_block_ends", stratum_block_ends,
                            "Phase block ends per het SNP")


def fig12c(db, ctx):
    return reference_strata(db, ctx, "fig12c_reference_strata_reads_tagged", stratum_tagged,
                            "Reads tagged", higher_is_better=True)


def fig12d(db, ctx):
    return reference_strata(db, ctx, "fig12d_reference_strata_wrong_haplotype", stratum_wrong,
                            "Tagged reads on the wrong haplotype")


FIGURES_EXTRA = [
    ("fig7d", fig7d, True),
    ("fig7e", fig7e, True),
    ("fig7f", fig7f, True),
    ("fig7g", fig7g, True),
    ("fig10a", fig10a, True),
    ("fig10b", fig10b, True),
    ("fig10c", fig10c, True),
    ("fig10d", fig10d, True),
    ("fig10e", fig10e, True),
    ("fig11a", fig11a, True),
    ("fig11b", fig11b, True),
    ("fig11c", fig11c, True),
    ("fig12a", fig12a, False),
    ("fig12b", fig12b, False),
    ("fig12c", fig12c, False),
    ("fig12d", fig12d, False),
]
