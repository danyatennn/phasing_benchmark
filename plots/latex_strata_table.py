#!/usr/bin/env python3
"""LaTeX tables (booktabs) of precision, recall and F1 per genomic region at one
coverage, for the variant callers and for one phaser on every input, one table
per reference. Every number is checked against the counts it comes from before
it is written; a mismatch stops the script (--force writes anyway and marks it).

  Variant calling  metrics_strata_auto/summary_calls_strata.tsv
                   precision = shared / emitted, recall = shared / truth_total
                   (shared: the caller's PASS het SNPs that are truth het SNPs)
  Phasing          metrics_strata_auto/summary_phase_strata.tsv
                   precision = correct / emitted, recall = correct / truth_total
                   (correct: truth het SNPs phased right, covered - hamming)
  F1 = 2PR / (P + R) in both, recomputed here.

  python plots/latex_strata_table.py --phaser whatshap --out tables/strata_prf.tex
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import style as st
from data import sample_root

REGIONS = [
    ("notinalldifficultregions", "Easy regions"),
    ("alldifficultregions", "All difficult"),
    ("SimpleRepeat_homopolymer_4to6_slop5", "Homopolymer 4--6 bp"),
    ("SimpleRepeat_homopolymer_7to11_slop5", "Homopolymer 7--11 bp"),
    ("SimpleRepeat_homopolymer_ge12_slop5", "Homopolymer $\\geq$12 bp"),
    ("SimpleRepeat_imperfecthomopolge11_slop5", "Imperfect homopolymer $\\geq$11 bp"),
    ("AllTandemRepeats_le50bp_slop5", "TR $\\leq$50 bp"),
    ("AllTandemRepeats_51to200bp_slop5", "TR 51--200 bp"),
    ("AllTandemRepeats_201to10000bp_slop5", "TR 201 bp--10 kb"),
    ("AllTandemRepeats_ge10001bp_slop5", "TR $>$10 kb"),
    ("AllTandemRepeatsandHomopolymers_slop5", "TR and homopolymers"),
    ("segdups", "Segmental duplications"),
    ("segdups_gt10kb", "SegDups $>$10 kb"),
    ("lowmappabilityall", "Low mappability"),
    ("satellites_slop5", "Satellites"),
    ("MHC", "MHC"),
    ("KIR", "KIR"),
    ("refseq_cds", "CDS"),
]
PLATFORM_TEX = {"ont4k": "ONT 4 kHz", "ont5k": "ONT 5 kHz", "hifi": "PacBio Sequel II",
                "revio": "PacBio Revio"}
TOL = 6e-4          # the TSVs carry 6 decimals; a stored rate may differ by rounding


def read(path):
    if not os.path.exists(path):
        sys.exit("missing %s" % path)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return math.nan
    return x


def f1(p, r):
    return 2 * p * r / (p + r) if p + r > 0 else 0.0


def collect(rows, kind, depth, phaser, problems):
    """(ref, input, platform, stratum) -> (P, R, F1, n_truth); checks every row"""
    out = {}
    for r in rows:
        if int(float(r.get("depth", 0))) != depth:
            continue
        if kind == "phase" and r.get("phaser") != phaser:
            continue
        key = (r["reference"], r["caller"], r["platform"], r["stratum"])
        tt, em = num(r["truth_total"]), num(r["emitted"])
        if kind == "calls":
            hit = num(r["shared"])
            stored = (num(r["precision"]), num(r["recall"]), num(r["f1"]))
            # the counts must add up too
            if num(r["fp"]) != em - hit or num(r["fn"]) != tt - hit:
                problems.append("%s %s: fp/fn do not match emitted, shared, truth_total" % (kind, key))
        else:
            hit = num(r["correct"])
            stored = (num(r["precision_phased"]), num(r["recall_phased"]), math.nan)
        p = hit / em if em else 0.0
        rc = hit / tt if tt else 0.0
        f = f1(p, rc)
        for name, mine, theirs in zip(("precision", "recall", "f1"), (p, rc, f), stored):
            if theirs == theirs and abs(mine - theirs) > TOL:
                problems.append("%s %s: %s stored %.6f, recomputed %.6f" % (kind, key, name, theirs, mine))
        if not (0 <= p <= 1 and 0 <= rc <= 1) or hit > min(em, tt):
            problems.append("%s %s: counts out of range (hit %g, emitted %g, truth %g)" % (kind, key, hit, em, tt))
        if key in out:
            problems.append("%s %s: more than one row" % (kind, key))
        out[key] = (p, rc, f, tt)
    return out


def cell(v, best):
    if v is None:
        return "--"
    text = "%.3f" % v
    return "\\textbf{%s}" % text if best else text


def table(ref, sections, platforms, depth, phaser):
    """sections: [(title, inputs, values)]"""
    cols = "ll" + "ccc" * len(platforms)
    lines = [
        "\\begin{table}[htbp]",
        "\\centering",
        "\\scriptsize",
        "\\setlength{\\tabcolsep}{3pt}",
        "\\caption{Het SNP calling and phasing across genomic regions, %s, %d$\\times$: precision (P), "
        "recall (R) and F1. Calling: the caller's PASS biallelic het SNPs against the truth. Phasing (%s): "
        "truth het SNPs phased correctly against the run's records (P) and against all truth het SNPs (R). "
        "Best value per region, platform and metric in bold.}" % (
            st.REFERENCE_LABEL[ref], depth, st.PHASER_LABEL[phaser]),
        "\\label{tab:strata_prf_%s}" % ref.lower(),
        "\\begin{tabular}{%s}" % cols,
        "\\toprule",
        " & & " + " & ".join("\\multicolumn{3}{c}{%s}" % PLATFORM_TEX[p] for p in platforms) + " \\\\",
        "".join("\\cmidrule(lr){%d-%d}" % (3 + 3 * i, 5 + 3 * i) for i in range(len(platforms))),
        "Region & Input & " + " & ".join(["P & R & F1"] * len(platforms)) + " \\\\",
        "\\midrule",
    ]
    for si, (title, inputs, values) in enumerate(sections):
        if si:
            lines.append("\\midrule")
        lines.append("\\multicolumn{%d}{l}{\\textit{%s}} \\\\" % (2 + 3 * len(platforms), title))
        for stratum, label in REGIONS:
            present = [c for c in inputs if any((ref, c, p, stratum) in values for p in platforms)]
            if not present:
                continue
            for n, c in enumerate(present):
                row = [label if n == 0 else "", st.CALLER_LABEL[c]]
                for p in platforms:
                    for m in range(3):
                        v = values.get((ref, c, p, stratum))
                        col = [values[(ref, x, p, stratum)][m] for x in present
                               if (ref, x, p, stratum) in values]
                        best = v is not None and len(col) > 1 and abs(v[m] - max(col)) < 5e-4
                        row.append(cell(v[m] if v else None, best))
                lines.append(" & ".join(row) + " \\\\")
            lines.append("\\addlinespace[2pt]")
    lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=st.ROOT)
    ap.add_argument("--sample", default=st.SAMPLE)
    ap.add_argument("--depth", type=int, default=30)
    ap.add_argument("--phaser", default="whatshap")
    ap.add_argument("--refs", default=" ".join(st.REFERENCES))
    ap.add_argument("--out", default=None, help="default: print")
    ap.add_argument("--force", action="store_true", help="write even if a check fails")
    args = ap.parse_args()

    base = os.path.join(sample_root(args.root, args.sample), "metrics_strata_auto")
    problems = []
    calls = collect(read(os.path.join(base, "summary_calls_strata.tsv")), "calls", args.depth,
                    args.phaser, problems)
    phase = collect(read(os.path.join(base, "summary_phase_strata.tsv")), "phase", args.depth,
                    args.phaser, problems)
    # the truth behind a region must be the same for every input of a cell
    totals = {}
    for (ref, c, p, s_), v in list(calls.items()) + list(phase.items()):
        totals.setdefault((ref, p, s_), set()).add(v[3])
    for k, t in totals.items():
        if len(t) > 1:
            problems.append("truth_total differs between inputs for %s: %s" % (k, sorted(t)))

    out = []
    for ref in args.refs.split():
        platforms = [p for p in st.PLATFORMS if any(k[0] == ref and k[2] == p for k in calls)
                     or any(k[0] == ref and k[2] == p for k in phase)]
        if not platforms:
            print("no %s rows at %dx" % (ref, args.depth), file=sys.stderr)
            continue
        callers = [c for c in st.CALLERS if c != "truth"]
        sections = [("Variant calling", callers, calls),
                    ("Phasing with %s" % st.PHASER_LABEL[args.phaser], st.CALLERS, phase)]
        out.append(table(ref, sections, platforms, args.depth, args.phaser))
        missing = [(c, p, s_) for _, cs, vals in sections for c in cs for p in platforms
                   for s_, _ in REGIONS
                   if (ref, c, p, s_) not in vals and any(k[3] == s_ and k[0] == ref for k in vals)]
        if missing:
            print("%s: %d cells without a row (written as --), e.g. %s" % (ref, len(missing), missing[:3]),
                  file=sys.stderr)

    print("checked %d calling and %d phasing rows: %d problems" % (len(calls), len(phase), len(problems)),
          file=sys.stderr)
    for p in problems[:30]:
        print("  " + p, file=sys.stderr)
    if problems and not args.force:
        sys.exit("not written: fix the problems above or pass --force")
    text = "% needs \\usepackage{booktabs}\n" + "\n".join(out)
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write(text)
        print("written %s" % args.out, file=sys.stderr)
    else:
        print(text)


if __name__ == "__main__":
    main()
