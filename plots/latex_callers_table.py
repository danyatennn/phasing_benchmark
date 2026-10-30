#!/usr/bin/env python3
"""LaTeX three-line table (booktabs: top, header and bottom rules only) of the
variant callers per genomic region at one coverage, one table per reference:
het SNP calling precision, recall and F1, and the read haplotag recall the
caller's SNPs give one phaser. Every number is recomputed from its counts and
checked before it is written; a mismatch stops the script (--force writes
anyway).

  Calling          precision = shared / emitted, recall = shared / truth
                   (biallelic het SNPs, PASS, inside the confident BED)
                   regions: metrics_strata_auto/summary_calls_strata.tsv
                   whole genome (chr1-22): metrics/summary_calls.tsv
  Haplotag recall  reads put on their own parental haplotype / reads the
                   truth labels (metrics_reads/summary_recall_strata_<truth>.tsv,
                   bin/read_recall_strata.sh READ_TRUTH=<truth>): primary
                   alignments present on both references, >= 2 truth het SNPs,
                   regions and chr1-22 from the T2T alignment - the same reads
                   in every cell of a platform. The GIAB truth row is the
                   ceiling: the phaser given the truth SNPs.

  python plots/latex_callers_table.py --phaser whatshap --read-truth matpat --out tables/callers_regions.tex
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import style as st
from data import sample_root

REGIONS = [
    ("genome", "Whole genome"),
    ("notinalldifficultregions", "Easy regions"),
    ("alldifficultregions", "All difficult"),
    ("segdups", "Segmental dup."),
    ("lowmappabilityall", "Low mappability"),
    ("AllTandemRepeatsandHomopolymers_slop5", "TR and homopol."),
    ("satellites_slop5", "Satellites"),
    ("MHC", "MHC"),
]
PLATFORM_TEX = {"ont4k": "ONT 4 kHz", "ont5k": "ONT 5 kHz", "hifi": "PacBio Sequel II",
                "revio": "PacBio Revio"}
METRICS = ["P", "R", "F1", "HR"]
TOL = 6e-4          # the TSVs carry 6 decimals
TRUTH_TEX = {"matpat": "the parental copy of the HG002 Q100 assembly they align to better",
             "yak": "parental $k$-mers", "truthtags": "the GIAB truth VCF"}


def read(path):
    if not os.path.exists(path):
        sys.exit("missing %s" % path)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def f1(p, r):
    return 2 * p * r / (p + r) if p + r > 0 else 0.0


def calls(base, depth, problems):
    """(ref, caller, platform, region) -> (P, R, F1), all checked"""
    out = {}

    def put(key, truth, emitted, shared, stored):
        if not (0 <= shared <= min(truth, emitted)) or not truth or not emitted:
            problems.append("calls %s: counts out of range (shared %g, emitted %g, truth %g)"
                            % (key, shared, emitted, truth))
            return
        p, r = shared / emitted, shared / truth
        for name, mine, theirs in zip(("precision", "recall", "f1"), (p, r, f1(p, r)), stored):
            if theirs == theirs and abs(mine - theirs) > TOL:
                problems.append("calls %s: %s stored %.6f, recomputed %.6f" % (key, name, theirs, mine))
        if key in out:
            problems.append("calls %s: more than one row" % (key,))
        out[key] = (100 * p, 100 * r, 100 * f1(p, r), truth)

    for r in read(os.path.join(base, "metrics_strata_auto", "summary_calls_strata.tsv")):
        if int(num(r.get("depth"))) != depth:
            continue
        tt, em, sh = num(r["truth_total"]), num(r["emitted"]), num(r["shared"])
        if num(r["fp"]) != em - sh or num(r["fn"]) != tt - sh:
            problems.append("calls %s %s %s %s: fp/fn do not match the counts"
                            % (r["reference"], r["caller"], r["platform"], r["stratum"]))
        put((r["reference"], r["caller"], r["platform"], r["stratum"]), tt, em, sh,
            (num(r["precision"]), num(r["recall"]), num(r["f1"])))
    for r in read(os.path.join(base, "metrics", "summary_calls.tsv")):
        if int(num(r.get("depth"))) != depth:
            continue
        tt, em, sh = num(r["truth_het_snps"]), num(r["calls"]), num(r["shared"])
        if num(r["fp"]) != em - sh or num(r["fn"]) != tt - sh:
            problems.append("calls %s %s %s genome: fp/fn do not match the counts"
                            % (r["reference"], r["caller"], r["platform"]))
        put((r["reference"], r["caller"], r["platform"], "genome"), tt, em, sh,
            (1 - num(r["fdr"]), 1 - num(r["fnr"]), math.nan))
    return out


def haplotag(base, truth, phaser, depth, problems):
    """(ref, caller, platform, region) -> (HR, reads); every input of a platform
    and region must have the same denominator, on both references for a
    parental truth"""
    out, denom = {}, {}
    for r in read(os.path.join(base, "metrics_reads", "summary_recall_strata_%s.tsv" % truth
                               if truth != "truthtags" else "summary_recall_strata.tsv")):
        if int(num(r.get("depth"))) != depth or r["phaser"] != phaser:
            continue
        key = (r["reference"], r["caller"], r["platform"], r["stratum"])
        tt, cp, cr = num(r["truth_tagged"]), num(r["compared"]), num(r["correct"])
        if not (0 <= cr <= cp <= num(r["tagged"]) <= tt) or not tt:
            problems.append("haplotag %s: counts out of order (correct %g, compared %g, tagged %s, "
                            "truth %g)" % (key, cr, cp, r["tagged"], tt))
            continue
        if abs(cr / tt - num(r["recall"])) > TOL:
            problems.append("haplotag %s: recall stored %s, recomputed %.6f" % (key, r["recall"], cr / tt))
        if key in out:
            problems.append("haplotag %s: more than one row" % (key,))
        out[key] = (100 * cr / tt, tt)
        side = (r["platform"], r["stratum"]) if truth != "truthtags" else \
            (r["reference"], r["platform"], r["stratum"])
        denom.setdefault(side, set()).add(tt)
    for side, n in denom.items():
        if len(n) > 1:
            problems.append("haplotag %s: reads the truth labels differ between inputs%s: %s"
                            % (side, "" if truth == "truthtags" else " or references", sorted(n)))
    return out


def fmt(v, best, decimals):
    if v is None:
        return "\\multicolumn{1}{c}{--}"
    text = "%.*f" % (decimals, v)
    return "\\textbf{%s}" % text if best else text


def table(ref, platforms, values, regions, callers, args):
    """values: (ref, caller, platform, region) -> {metric: value}"""
    ncol = 1 + len(METRICS) * len(platforms)
    # a wider gap between platforms than within, so the groups read without rules
    group = "@{\\hspace{6pt}}" + "r" * len(METRICS)
    lines = [
        "\\begin{table}[htbp]",
        "\\centering",
        "\\footnotesize",
        "\\setlength{\\tabcolsep}{2.4pt}",
        "\\setlength{\\belowcaptionskip}{6pt}",
        "\\caption{Het SNP calling and read haplotagging across genomic regions, %s, %d$\\times$. "
        "P, R, F1: precision, recall and F1 (\\%%) of the caller's biallelic het SNPs against the GIAB "
        "truth inside its confident regions. HR: haplotag recall (\\%%), the share of labelled reads that "
        "%s, phasing the caller's SNPs, puts on their true parental haplotype; reads labelled by %s, "
        "primary alignments present on both references with $\\geq$2 truth het SNPs, assigned to regions "
        "by their T2T-CHM13 alignment, so each platform has the same reads on both references. GIAB truth: "
        "the phaser given the truth SNPs (upper bound). Best caller per region, platform and metric in "
        "bold.}" % (
            st.REFERENCE_LABEL[ref], args.depth, st.PHASER_LABEL[args.phaser],
            TRUTH_TEX[args.read_truth]),
        "\\label{tab:callers_regions_%s}" % ref.lower(),
        "\\begin{tabular}{@{}l%s@{}}" % "".join(group for _ in platforms),
        "\\toprule",
        " & " + " & ".join("\\multicolumn{%d}{c}{%s}" % (len(METRICS), PLATFORM_TEX[p])
                             for p in platforms) + " \\\\",
        "Region / input & " + " & ".join(
            " & ".join("\\multicolumn{1}{c}{%s}" % m for m in METRICS) for _ in platforms) + " \\\\",
        "\\midrule",
    ]
    for gi, (region, label) in enumerate(regions):
        if gi:
            lines.append("\\addlinespace[4pt]")
        # the region on a line of its own, its inputs indented under it
        lines.append("\\multicolumn{%d}{@{}l}{\\textit{%s}} \\\\" % (ncol, label))
        for c in callers:
            row = ["\\quad " + st.CALLER_LABEL[c]]
            for p in platforms:
                v = values.get((ref, c, p, region), {})
                for m in METRICS:
                    x = v.get(m)
                    col = [values.get((ref, y, p, region), {}).get(m) for y in callers if y != "truth"]
                    col = [y for y in col if y is not None]
                    # best as printed: ties at the shown precision are all bold
                    best = (c != "truth" and x is not None and len(col) > 1
                            and round(x, args.decimals) == round(max(col), args.decimals))
                    row.append(fmt(x, best, args.decimals))
            lines.append(" & ".join(row) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=st.ROOT)
    ap.add_argument("--sample", default=st.SAMPLE)
    ap.add_argument("--depth", type=int, default=30)
    ap.add_argument("--phaser", default="whatshap", help="phaser behind HR")
    ap.add_argument("--read-truth", choices=["matpat", "yak", "truthtags"], default="matpat")
    ap.add_argument("--decimals", type=int, default=1, help="1 fits A4 with 2 cm margins")
    ap.add_argument("--refs", default=" ".join(st.REFERENCES))
    ap.add_argument("--out", default=None, help="default: print")
    ap.add_argument("--force", action="store_true", help="write even if a check fails")
    args = ap.parse_args()

    base = sample_root(args.root, args.sample)
    problems = []
    called = calls(base, args.depth, problems)
    tagged = haplotag(base, args.read_truth, args.phaser, args.depth, problems)
    # the truth behind a calling cell must be the same for every caller
    totals = {}
    for (ref, c, p, s_), v in called.items():
        totals.setdefault((ref, p, s_), set()).add(v[3])
    for k, t in totals.items():
        if len(t) > 1:
            problems.append("calls %s: truth het SNPs differ between callers: %s" % (k, sorted(t)))

    values = {}
    for k, (p, r, f, _) in called.items():
        values.setdefault(k, {}).update({"P": p, "R": r, "F1": f})
    for k, (hr, _) in tagged.items():
        values.setdefault(k, {})["HR"] = hr

    out = []
    for ref in args.refs.split():
        platforms = [p for p in st.PLATFORMS if any(k[0] == ref and k[2] == p for k in values)]
        if not platforms:
            print("no %s rows at %dx" % (ref, args.depth), file=sys.stderr)
            continue
        regions = [(s_, l) for s_, l in REGIONS if any(k[0] == ref and k[3] == s_ for k in values)]
        callers = [c for c in st.CALLERS if any(k[0] == ref and k[1] == c for k in values)]
        out.append(table(ref, platforms, values, regions, callers, args))
        holes = [(c, p, s_, m) for c in callers for p in platforms for s_, _ in regions for m in METRICS
                 if not (c == "truth" and m != "HR") and m not in values.get((ref, c, p, s_), {})]
        if holes:
            print("%s: %d cells without a value (written as --), e.g. %s" % (ref, len(holes), holes[:3]),
                  file=sys.stderr)

    print("checked %d calling and %d haplotag rows: %d problems" % (len(called), len(tagged), len(problems)),
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
