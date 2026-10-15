#!/usr/bin/env python3
"""Draft figures from whatever is in metrics/ right now.

  python plots/make_plots.py                 every figure, PDFs in plots/pdf/
  python plots/make_plots.py --only fig1,fig5
  python plots/make_plots.py --depth 20      the coverage of the single-depth figures
  python plots/make_plots.py --png           a PNG next to each PDF, for a quick look
  python plots/make_plots.py --refs GRCh38   one reference only (GRCh38 against T2T needs both)

One PDF per subfigure, text kept as editable type. A run that is not computed
yet, or whose row is locked by a host, is left out; a figure with nothing to
show is skipped and the rest are still drawn.
"""

import argparse
import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt

import style as st
from data import Metrics
from figures import FIGURES, NoData


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--out", default=None, help="default <root>/plots/pdf")
    parser.add_argument("--only", default="", help="comma-separated prefixes, e.g. fig1,fig5a")
    parser.add_argument("--depth", type=int, default=30, choices=st.DEPTHS)
    parser.add_argument("--png", action="store_true")
    parser.add_argument("--refs", default=" ".join(st.REFERENCES),
                        help="references to draw, default: %(default)s")
    ctx = parser.parse_args()
    ctx.refs = [r for r in ctx.refs.split() if r in st.REFERENCES]
    both = all(r in ctx.refs for r in st.REFERENCES)
    ctx.out = ctx.out or os.path.join(ctx.root, "plots", "pdf")

    font = st.apply_style(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    print("font      %s" % font)
    if font != st.FONT_STACK[0]:
        print("          %s is not installed here - put Arial.ttf into plots/fonts/ and rerun"
              % st.FONT_STACK[0])

    db = Metrics(ctx.root, st.SAMPLE, st.REGION).load(
        st.REFERENCES, st.CALLERS, st.PLATFORMS, st.DEPTHS, st.PHASERS)
    print("metrics   %s" % db.metrics_dir)
    for ref in st.REFERENCES:
        print("  %-7s computed %3d   locked %3d   not computed %3d" % (
            ref, db.count(ref, "done"), db.count(ref, "locked"), db.count(ref, "missing")))

    wanted = [w for w in ctx.only.split(",") if w]
    written = skipped = failed = 0
    print("figures   %s" % ctx.out)
    for fig_id, draw, per_reference in FIGURES:
        if wanted and not any(fig_id.startswith(w) for w in wanted):
            continue
        if not per_reference and not both:
            continue
        for ref in (ctx.refs if per_reference else [None]):
            label = fig_id + (" " + ref if ref else "")
            try:
                name = draw(db, ctx, ref) if per_reference else draw(db, ctx)
                print("  ok      %-14s %s" % (label, name))
                written += 1
            except NoData as why:
                print("  skip    %-14s %s" % (label, why))
                skipped += 1
            except Exception:
                # one broken figure must not take the others down with it
                print("  FAILED  %-14s" % label)
                traceback.print_exc()
                failed += 1
            finally:
                plt.close("all")

    print("%d written, %d skipped for lack of data, %d failed" % (written, skipped, failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
