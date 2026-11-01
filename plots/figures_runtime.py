#!/usr/bin/env python3
"""Run time and peak memory of the phasers at 30x (bin/runtime_run.sh), one
heatmap per metric and reference: a row per phaser, grouped by its input
(GIAB truth, then each caller), a column per platform, the value in each cell.
A cell with no measured run (HiPhase on ONT, a failed run) is hatched.

  figT1_runtime_30x_<ref>   wall-clock time of phase_one (min)
  figT2_memory_30x_<ref>    peak resident memory (GB)

From <root>/runtime/runtime_30x_<ref>.tsv; only runs with exit 0.

  python plots/figures_runtime.py
  python plots/figures_runtime.py --refs GRCh38 --png
"""

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt

import style as st
from figures import NoData, colorbar, heatmap, log_norm, paper_title, save

DEPTH = 30

# (file tag, title, value from a table row, colour bar label, cell text)
METRICS = [
    ("figT1_runtime", "Phasing run time at 30×", lambda r: float(r["wall_s"]) / 60, "min",
     lambda v: "%.0f" % v if v >= 10 else "%.1f" % v),
    ("figT2_memory", "Phasing peak memory at 30×", lambda r: float(r["max_rss_mb"]) / 1024, "GB",
     lambda v: "%.0f" % v if v >= 10 else "%.1f" % v),
]


def load(root, ref):
    """(caller, platform, phaser) -> table row of the last successful run"""
    path = os.path.join(root, "runtime", "runtime_%dx_%s.tsv" % (DEPTH, ref))
    rows = {}
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("exit") == "0" and r.get("depth") == str(DEPTH):
                rows[(r["caller"], r["platform"], r["phaser"])] = r
    return rows


def runtime_heatmap(ctx, ref, rows, tag, title, value, unit, text):
    callers = [c for c in st.CALLERS if any(k[0] == c for k in rows)]
    platforms = [p for p in st.PLATFORMS if any(k[1] == p for k in rows)]
    if not callers or not platforms:
        raise NoData("no measured run on %s - run bin/runtime_run.sh" % ref)
    matrix = [[value(rows[(c, p, ph)]) if (c, p, ph) in rows else math.nan for p in platforms]
              for c in callers for ph in st.PHASERS]
    values = [v for row in matrix for v in row if math.isfinite(v) and v > 0]
    if not values:
        raise NoData("no positive %s values on %s" % (tag, ref))

    nrows = len(callers) * len(st.PHASERS)
    fig, ax = plt.subplots(figsize=(88 * st.MM, (24 + 3.2 * nrows) * st.MM), constrained_layout=True)
    im = heatmap(ax, matrix,
                 [st.PHASER_LABEL[ph] for _ in callers for ph in st.PHASERS],
                 [st.PLATFORM_LABEL[p] for p in platforms],
                 [(st.CALLER_LABEL[c], len(st.PHASERS)) for c in callers],
                 [("", len(platforms))],
                 label=text, cmap=st.SEQUENTIAL, norm=log_norm(values))
    ax.set_xticklabels([st.PLATFORM_LABEL[p] for p in platforms], rotation=30, ha="right",
                       rotation_mode="anchor")
    colorbar(fig, im, ax, unit, unit=unit)
    paper_title(fig, "%s, %s" % (title, st.REFERENCE_LABEL[ref]))
    return save(fig, ctx, "%s_%dx_%s" % (tag, DEPTH, ref))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--out", default=None, help="default <root>/plots/pdf")
    parser.add_argument("--refs", default=" ".join(st.REFERENCES))
    parser.add_argument("--png", action="store_true", help="a PNG next to each PDF")
    ctx = parser.parse_args()
    ctx.out = ctx.out or os.path.join(ctx.root, "plots", "pdf")

    font = st.apply_style(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    if font != "Arial":
        print("note: Arial not found, using %s" % font, file=sys.stderr)

    for ref in ctx.refs.split():
        rows = load(ctx.root, ref)
        for tag, title, value, unit, text in METRICS:
            try:
                print("%s %s: %s" % (tag, ref, runtime_heatmap(ctx, ref, rows, tag, title, value, unit, text)))
            except NoData as e:
                print("%s %s: skipped, %s" % (tag, ref, e))
            finally:
                plt.close("all")


if __name__ == "__main__":
    main()
