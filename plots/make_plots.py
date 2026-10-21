#!/usr/bin/env python3
"""Draft figures from whatever is in metrics/ right now.

  python plots/make_plots.py                 every figure, PDFs in plots/pdf/
  python plots/make_plots.py --only fig1,fig5
  python plots/make_plots.py --depth 20      the coverage of the single-depth figures
  python plots/make_plots.py --png           a PNG next to each PDF, for a quick look
  python plots/make_plots.py --refs GRCh38   one reference only (GRCh38 against T2T needs both)
  python plots/make_plots.py --phasers "whatshap longphase hapcut2"
                                             these phasers only, e.g. while one is re-run
  python plots/make_plots.py --samples "HG002 HG005"   samples of the fig9 comparison
                                                       (default: HG002 and samples/*)

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
from data import sample_root
from figures import FIGURES, NoData
from figures_extra import FIGURES_EXTRA, sample_figures


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--out", default=None, help="default <root>/plots/pdf")
    parser.add_argument("--only", default="", help="comma-separated prefixes, e.g. fig1,fig5a")
    parser.add_argument("--depth", type=int, default=30, choices=st.DEPTHS)
    parser.add_argument("--png", action="store_true")
    parser.add_argument("--refs", default=" ".join(st.REFERENCES),
                        help="references to draw, default: %(default)s")
    parser.add_argument("--samples", default="",
                        help="samples of the fig9 comparison at 30x (default: HG002 and samples/*)")
    parser.add_argument("--phasers", default=" ".join(st.PHASERS),
                        help="phasers to draw, default: %(default)s")
    parser.add_argument("--rules", default="auto_min2",
                        help="read-accuracy rules of the tables to read (readacc_<rules>, ...)")
    ctx = parser.parse_args()
    ctx.refs = [r for r in ctx.refs.split() if r in st.REFERENCES]
    # in place: the figures read st.PHASERS when they are drawn
    st.PHASERS[:] = [p for p in st.PHASERS if p in ctx.phasers.split()]
    if not st.PHASERS:
        parser.error("--phasers names none of %s" % " ".join(st.PHASER_LABEL))
    both = all(r in ctx.refs for r in st.REFERENCES)
    ctx.out = ctx.out or os.path.join(ctx.root, "plots", "pdf")

    font = st.apply_style(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    print("font      %s" % font)
    if font != st.FONT_STACK[0]:
        print("          %s is not installed here - put Arial.ttf into plots/fonts/ and rerun"
              % st.FONT_STACK[0])

    db = Metrics(ctx.root, st.SAMPLE, st.REGION, readacc_rules=ctx.rules).load(
        st.REFERENCES, st.CALLERS, st.PLATFORMS, st.DEPTHS, st.PHASERS)
    print("metrics   %s" % db.metrics_dir)
    for ref in st.REFERENCES:
        print("  %-7s computed %3d   no autosome table %3d   locked %3d   not computed %3d" % (
            ref, db.count(ref, "done"), db.count(ref, "no_auto"), db.count(ref, "locked"),
            db.count(ref, "missing")))
        print("          NG50 over autosomes %3d   tagged bases %3d   read accuracy %3d" % (
            db.coverage(ref, "block_ng50"), db.coverage(ref, "bases_tagged_fraction"),
            db.coverage(ref, "read_precision")))

    wanted = [w for w in ctx.only.split(",") if w]
    written = skipped = failed = 0
    print("figures   %s" % ctx.out)
    for fig_id, draw, per_reference in FIGURES + FIGURES_EXTRA:
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

    # the other samples, side by side with HG002 at 30x on GRCh38
    if not wanted or any("fig9".startswith(w) or w.startswith("fig9") for w in wanted):
        names = ctx.samples.split() or (["HG002"] + sorted(
            d for d in os.listdir(os.path.join(ctx.root, "samples"))
            if os.path.isdir(os.path.join(ctx.root, "samples", d)))
            if os.path.isdir(os.path.join(ctx.root, "samples")) else ["HG002"])
        dbs = {}
        for name in names:
            if not os.path.isdir(sample_root(ctx.root, name)):
                continue
            dbs[name] = Metrics(ctx.root, name, st.REGION, readacc_rules=ctx.rules).load(
                ["GRCh38"], st.CALLERS, st.PLATFORMS, [30], st.PHASERS)
            print("  %-7s GRCh38 30x computed %3d" % (name, dbs[name].count("GRCh38", "done")))
        try:
            out = sample_figures(dbs, ctx)
            print("  ok      %-14s %s" % ("fig9", out))
            written += 1
        except NoData as why:
            print("  skip    %-14s %s" % ("fig9", why))
            skipped += 1
        except Exception:
            print("  FAILED  %-14s" % "fig9")
            traceback.print_exc()
            failed += 1
        finally:
            plt.close("all")

    print("%d written, %d skipped for lack of data, %d failed" % (written, skipped, failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
