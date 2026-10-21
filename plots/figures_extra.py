"""Figures of the read-level and difficult-region analyses, and of the other
GIAB samples. Same rules as figures.py: one PDF per subfigure, a figure with
nothing to show raises NoData.

  fig7d   read precision, blocks oriented by an anchor read against by SNPs
  fig7e   read recall against three truths (whatshap on the truth VCF,
          alignment to the parental Q100 copies, parental k-mers)
  fig10a  why reads are not tagged: share of reads per class
  fig10b  the same in bases
  fig10c  confidence of the tagged reads
  fig10d  accuracy of the tagged reads by confidence
  fig10e  what the reads of each class look like (one file per platform)
  fig11a  where the switch errors gather: enrichment per stratum
  fig11b  the same for phase block ends
  fig11c  why a stratum is hard: allele balance, paralog-like links, depth,
          low MAPQ, read errors (one file per platform)
  fig9a-d the other GIAB samples at 30x: switch rate, NG50, read recall, FNR
"""

import math

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import style as st
from data import NAN
from figures import (NoData, colorbar, empty_panel, finite, heatmap, log_ticks, panel_row,
                     phaser_axis, phaser_handles, platform_handles, positive, save,
                     short_number, strata_heatmap, stratum_label, top_legend, y_grid)

CLASS_LABEL = {
    "tagged": "Tagged",
    "lowmapq_paralog": "MAPQ < 20, paralogous region",
    "lowmapq_other": "MAPQ < 20, elsewhere",
    "no_snp_outside_bed": "Outside confident regions",
    "no_snp": "No het SNP on the read",
    "snps_unphased": "SNPs not phased (or single-SNP block)",
    "phased_unassigned": "SNPs phased, read not assigned",
}
CLASS_ORDER = list(CLASS_LABEL)
CLASS_COLOR = {
    "tagged": "#2C3A4B",
    "lowmapq_paralog": "#B5523B",
    "lowmapq_other": "#E39B7B",
    "no_snp_outside_bed": "#D9D9D9",
    "no_snp": "#A9B8C6",
    "snps_unphased": "#6D7F92",
    "phased_unassigned": "#E3A034",
}
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
            total = sum(table.get(("genome", gtype, g), {}).get(column, 0) or 0 for g in order)
            if not total:
                continue
            out.setdefault(platform, {})[phaser] = {
                g: 100 * (table.get(("genome", gtype, g), {}).get(column, 0) or 0) / total for g in order}
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


def fig10a(db, ctx, ref):
    return stacked_figure(db, ctx, ref, "fig10a_read_classes_reads", "reads", "class",
                          CLASS_ORDER, CLASS_LABEL, CLASS_COLOR, "Reads on chr1-22 (%)")


def fig10b(db, ctx, ref):
    return stacked_figure(db, ctx, ref, "fig10b_read_classes_bases", "bases", "class",
                          CLASS_ORDER, CLASS_LABEL, CLASS_COLOR, "Read bases on chr1-22 (%)")


def fig10c(db, ctx, ref):
    return stacked_figure(db, ctx, ref, "fig10c_tagged_confidence", "reads", "confidence",
                          CONF_ORDER, CONF_LABEL, CONF_COLOR, "Tagged reads (%)", tagged_only=True)


def fig10d(db, ctx, ref):
    """Accuracy of tagged reads in each confidence class, blocks oriented by
    their SNPs; a class with fewer than 20 compared reads is not drawn."""
    depth = ctx.depth
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), height=st.PANEL_HEIGHT + 6)
    drawn = []
    for ax, platform in zip(axes, st.PLATFORMS):
        here = 0
        for j, phaser in enumerate(st.PHASERS):
            table = db.classes((ref, "truth", platform, depth, phaser))
            for i, conf in enumerate(CONF_ORDER[:4]):
                row = table.get(("genome", "tagged_confidence", conf))
                if not row or not row.get("compared") or row["compared"] < 20:
                    continue
                err = 100 * (1 - row["correct"] / row["compared"])
                err = max(err, 0.01)
                ax.plot(i + (j - 1.5) * 0.15, err, marker=st.PHASER_MARKER[phaser],
                        color=st.PHASER_COLOR[phaser], markeredgecolor="white", linestyle="none", zorder=3)
                drawn.append(err)
                here += 1
        ax.set_xticks(range(4))
        ax.set_xticklabels([CONF_LABEL[c] for c in CONF_ORDER[:4]], rotation=35, ha="right",
                           rotation_mode="anchor")
        ax.set_xlim(-0.6, 3.6)
        ax.set_title(st.PLATFORM_LABEL[platform])
        y_grid(ax)
        if not here:
            empty_panel(ax)
    if not drawn:
        plt.close(fig)
        raise NoData("no accuracy by confidence for %s at %dx" % (ref, depth))
    axes[0].set_yscale("log")
    log_ticks(axes[0].yaxis, min(drawn), max(drawn))
    axes[0].set_ylabel("Tagged reads on the wrong\nhaplotype (%%), %d×" % depth)
    top_legend(legend_ax, phaser_handles(line=False))
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
        rows = [c for c in CLASS_ORDER if ("genome", "class", c) in table]
        if not rows:
            continue
        values = {c: table[("genome", "class", c)] for c in rows}
        fig, ax = plt.subplots(figsize=(st.FULL_WIDTH * 0.75 * st.MM, (30 + 7 * len(rows)) * st.MM),
                               constrained_layout=True)
        feature_heatmap(ax, fig, rows, [CLASS_LABEL[c] for c in rows], values,
                        "%s, %s, %d×" % (st.PLATFORM_LABEL[platform], st.PHASER_LABEL[phaser], depth))
        written.append(save(fig, ctx, "fig10e_class_features_%s_%dx_%s" % (platform, depth, ref)))
    if not written:
        raise NoData("no read classes for %s at %dx" % (ref, depth))
    return ", ".join(written)


def enrichment_heatmap(db, ctx, ref, name, metric, label):
    depth = ctx.depth
    mech = db.mechanism(ref)
    if not mech:
        raise NoData("no mechanism_%s.tsv - run bin/strata_mechanism.sh" % ref)
    strata = [s for s in db.strata_names()
              if any(k[4] == s and k[5] == metric for k in mech)]

    def value(platform, phaser, stratum):
        return mech.get((platform, depth, "truth", phaser, stratum, metric), NAN)

    values = [value(pl, p, s) for pl in st.PLATFORMS for p in st.PHASERS for s in strata]
    values = [v for v in values if finite(v) and v > 0]
    if not values:
        raise NoData("no %s for %s at %dx" % (metric, ref, depth))
    # symmetric about 1, so depletion and enrichment get the same colour depth
    span = max(abs(math.log10(v)) for v in values) or 1
    return strata_heatmap(ctx, "%s_%dx_%s" % (name, depth, ref), strata, value,
                          "%s, %d×" % (label, depth), cmap=DIVERGING,
                          norm=LogNorm(vmin=10 ** -span, vmax=10 ** span))


def fig11a(db, ctx, ref):
    return enrichment_heatmap(db, ctx, ref, "fig11a_switch_enrichment", "switch_enrichment",
                              "Switch errors per het SNP,\nrelative to the genome")


def fig11b(db, ctx, ref):
    return enrichment_heatmap(db, ctx, ref, "fig11b_block_end_enrichment", "break_enrichment",
                              "Phase block ends per het SNP,\nrelative to the genome")


STRATUM_FEATURES = [
    ("vaf_dev", "|VAF - 0.5|", lambda v: v),
    ("paralog_link_frac", "Links > 2 patterns (%)", lambda v: 100 * v),
    ("depth_ratio", "Depth / genome", lambda v: v),
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
                 for k in ("vaf_dev", "paralog_link_frac", "depth_ratio", "lowmapq_frac")}
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
    """Every run, read precision with blocks oriented by an anchor read against
    by their SNPs (both with blocks of two SNPs or more, chr1-22)."""
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
                    key = (ref, caller, platform, depth, phaser)
                    x = 100 * (1 - db.get(key, "read_precision"))
                    y = 100 * (1 - db.get(key, "read_precision_anchor_truthtags"))
                    if finite(x) and finite(y):
                        xs.append(max(x, 1e-3))
                        ys.append(max(y, 1e-3))
            if xs:
                ax.scatter(xs, ys, s=12, color=st.PLATFORM_COLOR[platform], marker=st.PHASER_MARKER[phaser],
                           edgecolor="white", linewidth=0.4, alpha=0.9, zorder=3)
                seen += xs + ys
    if not seen:
        plt.close(fig)
        raise NoData("no anchor-oriented accuracy for %s - run readacc_anchor.sh" % ref)
    lo, hi = min(seen) / 1.5, max(seen) * 1.5
    ax.plot([lo, hi], [lo, hi], color=st.MUTED, linewidth=0.6, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    log_ticks(ax.xaxis, lo, hi)
    log_ticks(ax.yaxis, lo, hi)
    ax.set_xlabel("Reads on the wrong haplotype (%),\nblocks oriented by their SNPs")
    ax.set_ylabel("Reads on the wrong haplotype (%),\nblocks oriented by an anchor read")
    ax.grid(color=st.GRID, linewidth=0.5, zorder=0)
    first = legend_ax.legend(handles=platform_handles(line=False), title="Platform", loc="upper left",
                             borderaxespad=0)
    legend_ax.add_artist(first)
    legend_ax.legend(handles=[Line2D([], [], color="none", marker=st.PHASER_MARKER[p], markerfacecolor=st.MUTED,
                                     markeredgecolor="white", label=st.PHASER_LABEL[p]) for p in st.PHASERS],
                     title="Phaser", loc="upper left", bbox_to_anchor=(0, 0.6), borderaxespad=0)
    return save(fig, ctx, "fig7d_anchor_vs_variant_orientation_" + ref)


TRUTH_LABEL = {"truthtags": "Truth VCF", "matpat": "Q100 copies", "yak": "Parental k-mers"}


def fig7e(db, ctx, ref):
    """Read recall of the truth-input runs against the three read-level truths,
    blocks oriented by an anchor read."""
    depth = ctx.depth
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), height=st.PANEL_HEIGHT + 10)
    drawn = []
    for ax, platform in zip(axes, st.PLATFORMS):
        for j, phaser in enumerate(st.PHASERS):
            key = (ref, "truth", platform, depth, phaser)
            for i, truth in enumerate(TRUTH_LABEL):
                v = 100 * db.get(key, "read_recall_anchor_%s" % truth)
                if finite(v):
                    ax.plot(i + (j - 1.5) * 0.15, v, marker=st.PHASER_MARKER[phaser],
                            color=st.PHASER_COLOR[phaser], markeredgecolor="white", linestyle="none", zorder=3)
                    drawn.append(v)
        ax.set_xticks(range(len(TRUTH_LABEL)))
        ax.set_xticklabels(list(TRUTH_LABEL.values()), rotation=35, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.6, len(TRUTH_LABEL) - 0.4)
        ax.set_title(st.PLATFORM_LABEL[platform])
        y_grid(ax)
    if not drawn:
        plt.close(fig)
        raise NoData("no anchor accuracy for %s at %dx" % (ref, depth))
    axes[0].set_ylabel("Labelled reads on their\ntruth haplotype (%%), %d×" % depth)
    top_legend(legend_ax, phaser_handles(line=False))
    return save(fig, ctx, "fig7e_read_recall_by_truth_%dx_%s" % (depth, ref))


# ------------------------------------------------------------- other samples

SAMPLE_METRICS = [
    ("fig9a_samples_switch_rate", "switch_rate", 100, "Switch error rate (%)", True, "phaser"),
    ("fig9b_samples_block_ng50", "block_ng50", 1e-6, "Phase block NG50 (Mb)", True, "phaser"),
    ("fig9c_samples_read_recall", "read_recall", 100, "Reads on the truth haplotype,\nof reads tagged by the truth (%)", False, "phaser"),
    ("fig9d_samples_caller_fnr", "fnr", 100, "Truth het SNPs missed (FNR, %)", True, "caller"),
]


def samples_figure(dbs, ctx, name, column, scale, ylabel, log, by, depth=30, ref="GRCh38"):
    samples = list(dbs)
    fig, legend_ax, axes = panel_row(len(st.PLATFORMS), height=st.PANEL_HEIGHT + 8)
    series = st.PHASERS if by == "phaser" else [c for c in st.CALLERS if c != "truth"]
    color = st.PHASER_COLOR if by == "phaser" else st.CALLER_COLOR
    marker = st.PHASER_MARKER if by == "phaser" else st.CALLER_MARKER
    drawn = []
    for ax, platform in zip(axes, st.PLATFORMS):
        for j, s in enumerate(series):
            xs, ys = [], []
            for i, sample in enumerate(samples):
                db = dbs[sample]
                if by == "phaser":
                    v = scale * db.get((ref, "truth", platform, depth, s), column)
                else:
                    v = scale * db.call(ref, s, platform, depth, column)
                if log:
                    v = positive(v)
                if finite(v):
                    xs.append(i + (j - (len(series) - 1) / 2) * 0.14)
                    ys.append(v)
            if xs:
                ax.plot(xs, ys, marker=marker[s], color=color[s], markeredgecolor="white",
                        linestyle="none", zorder=3)
                drawn += ys
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=35, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.6, len(samples) - 0.4)
        ax.set_title(st.PLATFORM_LABEL[platform])
        y_grid(ax)
    if not drawn:
        plt.close(fig)
        raise NoData("no %s for any sample at %dx" % (column, depth))
    if log:
        axes[0].set_yscale("log")
        log_ticks(axes[0].yaxis, min(drawn), max(drawn))
    axes[0].set_ylabel("%s, %d×" % (ylabel, depth))
    from figures import caller_handles
    top_legend(legend_ax, phaser_handles(line=False) if by == "phaser" else caller_handles(series, line=False))
    return save(fig, ctx, "%s_%dx_%s" % (name, depth, ref))


def sample_figures(dbs, ctx):
    out = []
    for name, column, scale, ylabel, log, by in SAMPLE_METRICS:
        try:
            out.append(samples_figure(dbs, ctx, name, column, scale, ylabel, log, by))
        except NoData:
            pass
        finally:
            plt.close("all")
    if not out:
        raise NoData("no sample has results at 30x")
    return ", ".join(out)


# (id, function, drawn once per reference?)
FIGURES_EXTRA = [
    ("fig7d", fig7d, True),
    ("fig7e", fig7e, True),
    ("fig10a", fig10a, True),
    ("fig10b", fig10b, True),
    ("fig10c", fig10c, True),
    ("fig10d", fig10d, True),
    ("fig10e", fig10e, True),
    ("fig11a", fig11a, True),
    ("fig11b", fig11b, True),
    ("fig11c", fig11c, True),
]
