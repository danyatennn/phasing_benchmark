#!/usr/bin/env python3
"""The tables of make_tables.py as a Supplementary Materials document: three-line
(booktabs) tables, captions above them, group labels written once with thin rules
between groups, sections S1-S6 with a short line of text and a table of contents.

  python plots/make_latex.py --tables tables/ --out tables/supplementary.tex
  pdflatex supplementary.tex     (twice, for the contents and table numbers)

  S1 Data and truth sets          sequencing data, truth sets
  S2 Phasing with the truth       headline metrics at 30x, error types, ranking, coverage
  S3 Variant calling and phasing  caller FNR / FDR, what calling and phasing lose
  S4 Read-level metrics           reads and bases tagged, precision and recall
  S5 GRCh38 against T2T-CHM13     pooled, per platform
  S6 Software

The full per-run grid (t04_grid) is far too wide for a page and is referred to as
Supplementary Data 1 (tables.xlsx). A table whose TSV is missing is left out. The
Wilcoxon p-values of the reference comparison are recomputed here from t04_grid
when scipy is available (make_tables.py leaves them empty without it).
"""

import argparse
import csv
import math
import os
import re

try:
    from scipy.stats import wilcoxon
except ImportError:
    wilcoxon = None

PLATFORMS = ["ONT 4 kHz", "ONT 5 kHz", "PacBio Sequel II", "PacBio Revio"]
REFS = ["GRCh38", "T2T-CHM13"]
INPUTS = ["GIAB truth", "Clair3", "DeepVariant", "longcallD"]
PHASER_ORDER = ["WhatsHap", "LongPhase", "Margin", "HapCUT2"]
METRIC_LABEL = {"switch_rate_pct": "Switch error rate", "ng50_mb": "NG50",
                "truth_phased_pct": "Phased", "read_recall_pct": "Read recall"}
METRIC_KIND = {"switch_rate_pct": "sig3", "ng50_mb": "sig3", "truth_phased_pct": "pct"}

# the provider alignments on GRCh38 (BAM headers; methodology 2.1)
ALIGNER = {"ONT 4 kHz": "minimap2 2.24", "ONT 5 kHz": "minimap2 2.24",
           "PacBio Sequel II": "pbmm2 1.1", "PacBio Revio": "pbmm2 1.9"}
X = "$\\times$"


def read(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def f(v, kind):
    """pct: 2 decimals; pct1: 1 decimal; sig3: 3 significant figures, trailing zeros
    kept; ratio: 2 decimals; int: thousands separated by commas; p: p-value."""
    x = num(v)
    if x is None:
        return "--"
    if kind == "pct":
        return "%.2f" % x
    if kind == "pct1":
        return ("%.1f" % x).replace("-0.0", "0.0")
    if kind == "sig3":
        return "0" if x == 0 else (("%#.3g" % x).rstrip(".") if x < 1000 else "{:,.0f}".format(x))
    if kind == "ratio":
        return "%.2f" % x
    if kind == "int":
        return "{:,}".format(int(round(x)))
    if kind == "p":
        if x < 0.001:
            m, e = ("%.1e" % x).split("e")
            return "$%s\\times10^{%d}$" % (m, int(e))
        return "%.3f" % x
    raise ValueError(kind)


def tex(s):
    return (str(s).replace("\\", "\\textbackslash{}").replace("&", "\\&").replace("%", "\\%")
            .replace("_", "\\_").replace("#", "\\#"))


# --------------------------------------------------------------- layout

def stack(*lines):
    """A header cell over several lines, bottom-aligned."""
    return "\\begin{tabular}[b]{@{}c@{}}%s\\end{tabular}" % "\\\\".join(lines)


def layout(rows, levels):
    """rows: full cell lists sorted by their first `levels` columns. A change in
    the first column gets a \\midrule, in the second a thin rule from column 2 on,
    in the third a little space; a group label is written only on its first row."""
    out, prev = [], None
    for r in rows:
        if prev is not None:
            if levels >= 1 and r[0] != prev[0]:
                out.append(("MID",))
            elif levels >= 2 and r[1] != prev[1]:
                out.append(("CMID", 2))
            elif levels >= 3 and r[2] != prev[2]:
                out.append(("SPACE",))
        cells = list(r)
        if prev is not None:
            for i in range(levels):
                if r[:i + 1] == prev[:i + 1]:
                    cells[i] = ""
                else:
                    break
        out.append(cells)
        prev = r
    return out


def body_lines(body, ncols):
    lines = []
    for r in body:
        if isinstance(r, tuple):
            if r[0] == "MID":
                lines.append("\\midrule")
            elif r[0] == "CMID":
                lines.append("\\cmidrule(l){%d-%d}" % (r[1], ncols))
            elif r[0] == "SPACE":
                lines.append("\\addlinespace[3pt]")
        else:
            lines.append(" & ".join(r) + " \\\\")
    return lines


def table(label, caption, colspec, header, body, size="\\small", long=True, landscape=False,
          sep=None):
    """header: raw header lines (each ends with \\\\, rules allowed between them)."""
    spec = re.sub(r"\{[^}]*\}", "", colspec)        # p{9cm} is one column
    ncols = sum(spec.count(c) for c in "lcrp")
    out = []
    if landscape:
        out.append("\\begin{landscape}")
    if long:
        out += ["{" + size]
        if sep:
            out.append("\\setlength{\\tabcolsep}{%s}" % sep)
        out += ["\\begin{longtable}{@{}%s@{}}" % colspec,
                "\\caption{%s}\\label{%s}\\\\" % (caption, label), "\\toprule"] + header + \
               ["\\midrule", "\\endfirsthead",
                "\\multicolumn{%d}{@{}l}{\\tablabelfont Table \\thetable\\ \\normalfont\\small(continued)}\\\\[2pt]"
                % ncols, "\\toprule"] + header + \
               ["\\midrule", "\\endhead", "\\bottomrule", "\\endlastfoot"]
        out += body_lines(body, ncols) + ["\\end{longtable}", "}"]
    else:
        out += ["\\begin{table}[htbp]", "\\centering", size]
        if sep:
            out.append("\\setlength{\\tabcolsep}{%s}" % sep)
        out += ["\\caption{%s}\\label{%s}" % (caption, label),
                "\\begin{tabular}{@{}%s@{}}" % colspec, "\\toprule"] + header + ["\\midrule"]
        out += body_lines(body, ncols) + ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
    if landscape:
        out.append("\\end{landscape}")
    return "\n".join(out) + "\n"


def phasers_of(d):
    seen = {r["phaser"] for r in d}
    return [p for p in PHASER_ORDER if p in seen]


# ------------------------------------------------------------------ tables

def t_datasets(d):
    rows = []
    for r in d:
        if r["sample"] != "HG002":
            continue
        depth = lambda c: "%.1f" % num(r[c]) if num(r[c]) else "--"
        rows.append([tex(r["platform"]), ALIGNER.get(r["platform"], "--"),
                     depth("source_depth_grch38"), depth("source_depth_t2t"),
                     f(r["reads_30x"], "int"), f(r["mean_kb"], "pct1"), f(r["median_kb"], "pct1"),
                     f(r["read_n50_kb"], "pct1"), f(r["reads_ge10kb_pct"], "pct1"),
                     f(r.get("reads_ge50kb_pct"), "pct1")])
    cap = ("Sequencing data of HG002. Aligner: the provider's alignment to GRCh38; every data set "
           "was re-aligned to T2T-CHM13 with minimap2 2.31 (\\texttt{lr:hq} for ONT, "
           "\\texttt{map-hifi} for PacBio). Depth: mean depth of the full data set (mosdepth). "
           "Each data set was downsampled to 10, 20, 30, 40 and 50%s (\\texttt{samtools view -s}, "
           "fixed seed; achieved depth within 0.03%s of the target). Read statistics from the "
           "30%s GRCh38 data, primary alignments on chr1--22; read length counts the whole read, "
           "clipped ends included%s." % (X, X, X, " (how the reads align is in Table~\\ref{tab:alignment})"
                                            if any(num(r.get("aligned_pct")) is not None for r in d) else ""))
    header = [" & & \\multicolumn{2}{c}{Depth (%s)} & & \\multicolumn{3}{c}{Read length (kb)} & & \\\\" % X,
              "\\cmidrule(lr){3-4}\\cmidrule(lr){6-8}",
              "Platform & Aligner & GRCh38 & T2T & %s & Mean & Median & N50 & %s & %s \\\\"
              % (stack("Reads", "(30%s)" % X), stack("$\\geq$10 kb", "(\\% reads)"),
                 stack("$\\geq$50 kb", "(\\% reads)"))]
    return table("tab:datasets", cap, "llrrrrrrrr", header, rows, size="\\footnotesize", sep="4pt")


def t_alignment(d):
    rows = []
    for r in d:
        if r["sample"] != "HG002" or num(r.get("aligned_pct")) is None:
            continue
        rows.append([tex(r["platform"]), f(r["aligned_pct"], "pct1"), f(r["softclip_pct"], "pct1"),
                     f(r["read_n50_kb"], "pct1"), f(r["aligned_span_n50_kb"], "pct1"),
                     f(r["reads_ge50kb_pct"], "pct"), f(r["span_ge50kb_pct"], "pct"),
                     f(r["mapq_ge20_pct"], "pct1"), f(r["identity_pct"], "pct1")])
    if not rows:
        return None
    cap = ("How the reads of HG002 align (30%s GRCh38 data, primary alignments on chr1--22). "
           "Aligned and soft-clipped: read bases inside the primary alignment and clipped from its "
           "ends, of all read bases. Aligned span: end $-$ start of the primary alignment on the "
           "reference, the stretch over which a read can link heterozygous SNPs; a long read that "
           "aligns only in pieces has a long read length but a short aligned span. Identity = 1 $-$ "
           "NM / aligned bases (--: no NM tag in the PacBio BAMs)." % X)
    header = [" & \\multicolumn{2}{c}{Read bases (\\%)} & \\multicolumn{2}{c}{N50 (kb)} & "
              "\\multicolumn{2}{c}{$\\geq$50 kb (\\% reads)} & & \\\\",
              "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
              "Platform & Aligned & Soft-clipped & Read length & Aligned span & Read length & "
              "Aligned span & %s & %s \\\\" % (stack("MAPQ $\\geq$20", "(\\%)"), stack("Identity", "(\\%)"))]
    return table("tab:alignment", cap, "lrrrrrrrr", header, rows, size="\\footnotesize", sep="5pt")


def t_truth(d):
    rows = []
    for r in d:
        if r["sample"] != "HG002":
            continue
        rows.append([tex(r["reference"]), f(r["bed_mb"], "int"), f(r["het_snp_01"], "int"),
                     f(r["het_snp_12"], "int"), f(r["het_indel"], "int"), f(r["hom_alt_snp"], "int"),
                     f(r["hom_alt_indel"], "int"), f(r["benchmark_het_snps"], "int"),
                     f(r["benchmark_per_mb"], "pct1"), f(r["median_distance_bp"], "int")])
    cap = ("Truth sets: GIAB v5.0q small variants of HG002 (Q100-based) inside the benchmark "
           "regions, chr1--22. Het SNP 0/1: biallelic heterozygous SNPs; 1/2: heterozygous with "
           "two different ALT alleles; het indel: heterozygous records with an indel allele. "
           "Benchmark het SNPs: the phased biallelic heterozygous SNPs on which every phasing "
           "metric is computed (all het SNP 0/1 are phased in v5.0q), their number per Mb of "
           "benchmark regions and the median distance between consecutive ones.")
    header = [" & & \\multicolumn{5}{c}{Truth variants} & \\multicolumn{3}{c}{Benchmark het SNPs} \\\\",
              "\\cmidrule(lr){3-7}\\cmidrule(l){8-10}",
              "Reference & %s & %s & %s & %s & %s & %s & Number & Per Mb & %s \\\\"
              % (stack("Regions", "(Mb)"), stack("Het SNP", "0/1"), stack("Het SNP", "1/2"),
                 stack("Het", "indel"), stack("Hom-alt", "SNP"), stack("Hom-alt", "indel"),
                 stack("Median", "distance (bp)"))]
    return table("tab:truth", cap, "lrrrrrrrrr", header, rows, size="\\footnotesize", sep="3.5pt")


def t_main(d):
    cols = [("switch_rate_pct", "sig3"), ("ng50_mb", "sig3"), ("truth_phased_pct", "pct"),
            ("hamming_rate_pct", "sig3"), ("reads_tagged_pct", "pct1"), ("bases_tagged_pct", "pct1"),
            ("read_precision_snp_pct", "pct"), ("read_recall_snp_pct", "pct")]
    rows = []
    for ref in REFS:
        for p in PLATFORMS:
            for ph in phasers_of(d):
                for r in d:
                    if (r["reference"], r["platform"], r["phaser"]) == (ref, p, ph):
                        rows.append([tex(ref), tex(p), tex(ph)] + [f(r[c], k) for c, k in cols])
    cap = ("Phasing of HG002 at 30%s with the GIAB truth as input (every phaser receives the same "
           "unphased benchmark het SNPs). Switch: switch errors / assessed pairs of consecutive "
           "phased sites. NG50: phase block span at which the blocks cover half of the length of "
           "chr1--22. Phased: benchmark het SNPs in a block of $\\geq$2 SNPs. Hamming: blockwise "
           "Hamming distance / phased sites. Tagged: reads and bases assigned a haplotype in a "
           "block of $\\geq$2 SNPs (\\texttt{whatshap haplotag}). Precision and recall: reads on the "
           "truth haplotype, of reads tagged by both the run and the truth (precision) or by the "
           "truth (recall); blocks oriented by their SNPs." % X)
    header = [" & & & \\multicolumn{4}{c}{Phasing} & \\multicolumn{2}{c}{Tagged (\\%)} & "
              "\\multicolumn{2}{c}{Reads (\\%)} \\\\",
              "\\cmidrule(lr){4-7}\\cmidrule(lr){8-9}\\cmidrule(l){10-11}",
              "Reference & Platform & Phaser & Switch (\\%) & NG50 (Mb) & Phased (\\%) & Hamming (\\%) "
              "& Reads & Bases & Precision & Recall \\\\"]
    return table("tab:main", cap, "lll" + "r" * 8, header, layout(rows, 2), size="\\footnotesize",
                 long=True, sep="4pt")


def t_coverage(d):
    rows = []
    for ref in REFS:
        for r in d:
            if r["reference"] != ref:
                continue
            k = METRIC_KIND[r["metric"]]
            rows.append([tex(ref), tex(r["platform"]), tex(r["phaser"]), METRIC_LABEL[r["metric"]]] +
                        [f(r[x], k) for x in ("10x", "20x", "30x", "40x", "50x")] +
                        [f(r["change_10_30_pct"], "pct1"), f(r["change_30_50_pct"], "pct1"),
                         (r["saturation_depth"] + X) if r["saturation_depth"] else "--"])
    cap = ("Phasing metrics against coverage (GIAB truth input, chr1--22): switch error rate (\\%%), "
           "NG50 (Mb) and share of benchmark het SNPs phased (\\%%). Change: from 10%s to 30%s and "
           "from 30%s to 50%s, in \\%% of the earlier value. Saturation: the lowest coverage within "
           "5\\%% of the best value over all coverages." % (X, X, X, X))
    header = [" & & & & \\multicolumn{5}{c}{Coverage} & \\multicolumn{2}{c}{Change (\\%)} & \\\\",
              "\\cmidrule(lr){5-9}\\cmidrule(lr){10-11}",
              "Reference & Platform & Phaser & Metric & 10%s & 20%s & 30%s & 40%s & 50%s & "
              "10$\\to$30 & 30$\\to$50 & Saturation \\\\" % (X, X, X, X, X)]
    return table("tab:coverage", cap, "llll" + "r" * 7 + "c", header, layout(rows, 3),
                 size="\\footnotesize", long=True, landscape=True, sep="5pt")


def t_error_types(d):
    idx = {(r["reference"], r["platform"], r["phaser"]): r for r in d}
    rows = []
    for p in PLATFORMS:
        for ph in phasers_of(d):
            cells = [idx.get((ref, p, ph)) for ref in REFS]
            if not any(cells):
                continue
            row = [tex(p), tex(ph)]
            for c in cells:
                row += [f(c["switch_only_per100"], "sig3"), f(c["flip_per100"], "sig3"),
                        f(c["flip_share_pct"], "pct1")] if c else ["--"] * 3
            rows.append(row)
    cap = ("Switch and flip errors at 30%s (GIAB truth input, chr1--22), per 100 assessed pairs of "
           "consecutive phased sites. A flip is a single site phased the wrong way (two switches "
           "in a row, counted once); switch: the remaining switches, a stretch phased the wrong "
           "way. Flip share: flips / all errors." % X)
    header = [" & & \\multicolumn{3}{c}{GRCh38} & \\multicolumn{3}{c}{T2T-CHM13} \\\\",
              "\\cmidrule(lr){3-5}\\cmidrule(l){6-8}",
              "Platform & Phaser & Switch & Flip & Flip share (\\%) & Switch & Flip & Flip share (\\%) \\\\"]
    return table("tab:error-types", cap, "ll" + "r" * 6, header, layout(rows, 1))


def t_ranking(d):
    phasers = phasers_of(d)
    idx = {(r["reference"], r["input"], r["metric"], r["phaser"]): r for r in d}
    rows = []
    for ref in REFS:
        for i in INPUTS:
            for m in ["switch_rate_pct", "ng50_mb", "truth_phased_pct", "read_recall_pct"]:
                cells = [idx.get((ref, i, m, ph)) for ph in phasers]
                if not any(cells):
                    continue
                rows.append([tex(ref), tex(i), METRIC_LABEL[m]] +
                            ["%s/%s (%s)" % (c["wins"], c["cells"], "%g" % num(c["median_rank"]))
                             if c else "--" for c in cells])
    cap = ("Phaser ranking. Within every platform $\\times$ coverage cell (20 per input) the phasers "
           "are ranked by the metric: lowest switch error rate, longest NG50, highest share phased, "
           "highest read recall. Each cell: number of cells where the phaser ranks first / cells "
           "(median rank).")
    header = ["Reference & Input & Metric & " + " & ".join(tex(p) for p in phasers) + " \\\\"]
    return table("tab:ranking", cap, "lll" + "c" * len(phasers), header, layout(rows, 2),
                 size="\\footnotesize", long=True)


def t_callers(d):
    depths = ["10", "20", "30", "40", "50"]
    idx = {(r["reference"], r["caller"], r["platform"], r["depth"]): r for r in d}
    rows = []
    for ref in REFS:
        for c in INPUTS[1:]:
            for p in PLATFORMS:
                cells = [idx.get((ref, c, p, x)) for x in depths]
                if not any(cells):
                    continue
                rows.append([tex(ref), tex(c), tex(p)] +
                            [f(x["fnr_pct"], "pct") if x else "--" for x in cells] +
                            [f(x["fdr_pct"], "pct") if x else "--" for x in cells])
    cap = ("Variant callers: false negative rate (FNR) and false discovery rate (FDR) of biallelic "
           "heterozygous PASS SNPs against the benchmark het SNPs, inside the benchmark regions on "
           "chr1--22, positions with more than one record excluded. FNR = FN / benchmark het SNPs; "
           "FDR = FP / calls.")
    header = [" & & & \\multicolumn{5}{c}{FNR (\\%)} & \\multicolumn{5}{c}{FDR (\\%)} \\\\",
              "\\cmidrule(lr){4-8}\\cmidrule(l){9-13}",
              "Reference & Caller & Platform & " + " & ".join(["%s%s" % (x, X) for x in depths] * 2) + " \\\\"]
    return table("tab:callers", cap, "lll" + "r" * 10, header, layout(rows, 2), size="\\footnotesize",
                 long=True, sep="4pt")


def t_call_vs_phase(d):
    idx = {(r["reference"], r["platform"], r["input"], r["phaser"]): r for r in d}
    rows = []
    for p in PLATFORMS:
        for i in INPUTS:
            for ph in phasers_of(d):
                cells = [idx.get((ref, p, i, ph)) for ref in REFS]
                if not any(cells):
                    continue
                row = [tex(p), tex(i), tex(ph)]
                for c in cells:
                    row += [f(c["called_pct"], "pct"), f(c["correct_pct"], "pct"),
                            f(c["lost_calling_pct"], "pct"), f(c["lost_phasing_pct"], "pct")] if c else ["--"] * 4
                rows.append(row)
    cap = ("What calling and what phasing lose at 30%s, all as a share of the benchmark het SNPs "
           "(chr1--22). Called: found by the caller; correct: found and phased correctly; lost in "
           "calling = 100 $-$ called; lost in phasing = called $-$ correct. Caller VCFs were given "
           "to the phasers unchanged (DeepVariant and longcallD write their own phase)." % X)
    header = [" & & & \\multicolumn{4}{c}{GRCh38 (\\%)} & \\multicolumn{4}{c}{T2T-CHM13 (\\%)} \\\\",
              "\\cmidrule(lr){4-7}\\cmidrule(l){8-11}",
              "Platform & Input & Phaser & Called & Correct & Lost calling & Lost phasing & "
              "Called & Correct & Lost calling & Lost phasing \\\\"]
    return table("tab:call-vs-phase", cap, "lll" + "r" * 8, header, layout(rows, 2),
                 size="\\footnotesize", long=True, sep="4pt")


def t_reads(d):
    rows = []
    for ref in REFS:
        for p in PLATFORMS:
            for i in INPUTS:
                for ph in phasers_of(d):
                    for r in d:
                        if (r["reference"], r["platform"], r["input"], r["phaser"]) == (ref, p, i, ph):
                            rows.append([tex(ref), tex(p), tex(i), tex(ph),
                                         f(r["reads_tagged_pct"], "pct1"), f(r["bases_tagged_pct"], "pct1"),
                                         f(r["read_precision_snp_pct"], "pct"), f(r["read_recall_snp_pct"], "pct"),
                                         f(r["read_precision_anchor_pct"], "pct"), f(r["read_recall_anchor_pct"], "pct")])
    cap = ("Read-level metrics at 30%s (chr1--22). Tagged: reads and bases assigned a haplotype in "
           "a block of $\\geq$2 SNPs. Precision: reads on the truth haplotype of reads tagged by "
           "both the run and the truth; recall: of reads tagged by the truth. Blocks oriented by "
           "their SNPs against the truth genotypes (SNP), or by an anchor read, the longest read "
           "of the block with a truth label (anchor; the anchor itself is not scored)." % X)
    header = [" & & & & \\multicolumn{2}{c}{Tagged (\\%)} & \\multicolumn{2}{c}{SNP orientation (\\%)} & "
              "\\multicolumn{2}{c}{Anchor orientation (\\%)} \\\\",
              "\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(l){9-10}",
              "Reference & Platform & Input & Phaser & Reads & Bases & Precision & Recall & Precision & Recall \\\\"]
    return table("tab:reads", cap, "llll" + "r" * 6, header, layout(rows, 3), size="\\footnotesize",
                 long=True, sep="4pt")


def wilcoxon_p(grid, col, platforms, phasers):
    if wilcoxon is None or grid is None:
        return None
    v = {}
    for r in grid:
        if r["input"] == "GIAB truth" and r["platform"] in platforms and r["phaser"] in phasers:
            x = num(r[col])
            if x is not None and x > 0:
                v[(r["reference"], r["platform"], r["depth"], r["phaser"])] = x
    d = [math.log(v[("T2T-CHM13",) + k[1:]] / v[k]) for k in v
         if k[0] == "GRCh38" and ("T2T-CHM13",) + k[1:] in v]
    if len(d) < 6:
        return None
    return 1.0 if all(x == 0 for x in d) else float(wilcoxon(d).pvalue)


def t_reference_pooled(d, grid):
    phasers = [p for p in PHASER_ORDER if grid and any(r["phaser"] == p for r in grid)]
    groups = {p: ([p], PLATFORMS) for p in phasers}
    groups.update({"ONT, all phasers": (phasers, PLATFORMS[:2]),
                   "PacBio, all phasers": (phasers, PLATFORMS[2:]),
                   "All runs": (phasers, PLATFORMS)})
    rows = []
    for r in d:
        ph, pl = groups.get(r["group"], (None, None))
        p = num(r["wilcoxon_p"])
        if p is None and ph is not None:
            p = wilcoxon_p(grid, r["metric"], pl, ph)
        rows.append([tex(r["group"]), METRIC_LABEL.get(r["metric"], tex(r["metric"])), r["pairs"],
                     f(r["median_ratio"], "ratio"),
                     "%s--%s" % (f(r["q1_ratio"], "ratio"), f(r["q3_ratio"], "ratio")),
                     r["t2t_better"], f(p, "p") if p is not None else "--"])
    cap = ("GRCh38 against T2T-CHM13, pooled. Ratio T2T-CHM13 / GRCh38 of the same run (same "
           "platform, coverage and phaser; GIAB truth input, chr1--22), pooled per phaser over the "
           "four platforms and five coverages, and over phasers per technology. T2T better: pairs "
           "with a lower switch error rate, a longer NG50 or a higher share phased on T2T-CHM13. "
           "$p$: two-sided Wilcoxon signed-rank test of log ratio = 0; descriptive, since the runs "
           "of one platform share their reads across coverages. Reads were aligned differently to "
           "the two references (provider alignments on GRCh38, minimap2 2.31 on T2T-CHM13), so the "
           "ratio contains both the reference and the alignment.")
    header = ["Group & Metric & Pairs & Median ratio & IQR & T2T better & $p$ \\\\"]
    return table("tab:reference", cap, "llrrcrr", header, layout(rows, 1))


def t_reference_platform(d):
    rows = []
    for r in d:
        k = METRIC_KIND[r["metric"]]
        rows.append([tex(r["platform"]), tex(r["phaser"]), METRIC_LABEL[r["metric"]],
                     f(r["grch38_30x"], k), f(r["t2t_30x"], k), f(r["ratio_30x"], "ratio"),
                     f(r["median_ratio"], "ratio"),
                     "%s--%s" % (f(r["min_ratio"], "ratio"), f(r["max_ratio"], "ratio")),
                     "%s/%s" % (r["t2t_better"], r["depths"])])
    cap = ("GRCh38 against T2T-CHM13 per platform and phaser (GIAB truth input, chr1--22). Values at "
           "30%s on both references (switch error rate and share phased in \\%%, NG50 in Mb) and their "
           "ratio T2T-CHM13 / GRCh38; over the five coverages, the median ratio, its range and the "
           "number of coverages at which T2T-CHM13 is better." % X)
    header = [" & & & \\multicolumn{3}{c}{30%s} & \\multicolumn{3}{c}{10--50%s} \\\\" % (X, X),
              "\\cmidrule(lr){4-6}\\cmidrule(l){7-9}",
              "Platform & Phaser & Metric & GRCh38 & T2T & Ratio & Median ratio & Range & T2T better \\\\"]
    return table("tab:reference-platform", cap, "lll" + "rrrrcc", header, layout(rows, 2),
                 size="\\footnotesize", long=True)


ROLE = [("whatshap", "WhatsHap", "phasing, haplotagging, evaluation (\\texttt{compare}, \\texttt{stats})"),
        ("LongPhase", "LongPhase", "phasing"), ("Margin", "Margin", "phasing"),
        ("hapcut2", "HapCUT2", "phasing"), ("Clair3", "Clair3", "variant calling"),
        ("DeepVariant", "DeepVariant", "variant calling"), ("longcallD", "longcallD", "variant calling"),
        ("minimap2", "minimap2", "alignment to T2T-CHM13 and to the Q100 assembly"),
        ("samtools", "SAMtools", "downsampling, BAM processing"),
        ("bcftools", "BCFtools", "VCF processing, truth preparation"), ("mosdepth", "mosdepth", "coverage"),
        ("hifiasm", "hifiasm", "assembly"), ("verkko", "Verkko", "trio assembly"),
        ("yak", "yak", "parental k-mers for hifiasm trio"), ("meryl", "Meryl", "parental k-mers for Verkko"),
        ("merqury", "Merqury", "hap-mers"), ("dipcall", "dipcall", "variants from assemblies")]


def t_software(d):
    have = {r["tool"]: r["version"] for r in d}
    rows = [[name, tex(have[k]) if have[k] not in ("", "?") else "to be filled", role]
            for k, name, role in ROLE if k in have]
    cap = ("Software. Clair3 models \\texttt{r1041\\_e82\\_400bps\\_sup\\_v410} (ONT 4 kHz), "
           "\\texttt{r1041\\_e82\\_400bps\\_sup\\_v500} (ONT 5 kHz), \\texttt{hifi\\_sequel2} and "
           "\\texttt{hifi\\_revio}; DeepVariant models \\texttt{ONT\\_R104} and \\texttt{PACBIO}; "
           "longcallD presets \\texttt{--ont} and \\texttt{--hifi}; Margin with its default "
           "parameter files (\\texttt{allParams.phase\\_vcf.ont.json}, "
           "\\texttt{allParams.phase\\_vcf.pb-hifi.json}).")
    return table("tab:software", cap, "llp{9cm}", ["Tool & Version & Used for \\\\"], rows)


PREAMBLE = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=2.4cm,headheight=14pt]{geometry}
\usepackage[T1]{fontenc}
\usepackage{newtxtext,newtxmath}
\usepackage[scaled=0.92]{helvet}
\usepackage{xcolor}
\definecolor{accent}{RGB}{31,95,107}
\usepackage{booktabs,longtable,array,pdflscape,caption,titlesec,fancyhdr}
\usepackage[hidelinks]{hyperref}
\usepackage{tocloft}

\captionsetup{font=small,labelfont={sf,bf,color=accent},labelsep=period,
  justification=justified,singlelinecheck=false,skip=6pt}
\captionsetup[table]{position=top}
\setlength{\LTcapwidth}{\linewidth}
\captionsetup[longtable]{width=\linewidth,margin=0pt}
\newcommand{\tablabelfont}{\small\sffamily\bfseries\color{accent}}
\renewcommand{\arraystretch}{1.05}
\setlength{\heavyrulewidth}{0.08em}
\setlength{\lightrulewidth}{0.05em}

\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thetable}{S\arabic{table}}
\titleformat{\section}{\Large\sffamily\bfseries\color{accent}}{\thesection}{0.8em}{}
\titlespacing*{\section}{0pt}{16pt}{8pt}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\scriptsize\sffamily\color{accent}Supplementary Materials}
\fancyfoot[C]{\small\thepage}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\headrule}{{\color{accent!40}\hrule height\headrulewidth}}

\renewcommand{\cfttoctitlefont}{\large\sffamily\bfseries\color{accent}}
\renewcommand{\cftsecfont}{\bfseries}
\renewcommand{\cftsecpagefont}{\bfseries}
\renewcommand{\cftsecleader}{\cftdotfill{\cftdotsep}}
\setlength{\cftsecnumwidth}{2.6em}

\begin{document}
\begin{center}{\LARGE\sffamily\bfseries\color{accent}Supplementary Materials}\end{center}
\vspace{2pt}{\color{accent!40}\hrule}\vspace{6pt}
\tableofcontents
\clearpage
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tables", required=True, help="folder with the TSVs of make_tables.py")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    T = lambda n: read(os.path.join(args.tables, n + ".tsv"))
    grid = T("t04_grid")
    margin = bool(grid) and any(r["phaser"] == "Margin" for r in grid)

    def ref(label):
        return "Table~\\ref{%s}" % label

    S = [PREAMBLE]
    # the alignment table needs the columns added to t01 in October 2026
    alignment = t_alignment(T("t01_datasets")) if T("t01_datasets") else None
    S.append("\\section{Data and truth sets}\n%s describes the sequencing data of HG002, %sand %s the "
             "truth sets. All metrics are computed on chr1--22.%s The metrics of every run "
             "(reference $\\times$ input $\\times$ platform $\\times$ coverage $\\times$ phaser) are "
             "in Supplementary Data~1 (\\texttt{tables.xlsx}, sheet \\texttt{t04\\_grid}).\n"
             % (ref("tab:datasets"), (ref("tab:alignment") + " how its reads align ") if alignment else "",
                ref("tab:truth"),
                "" if margin else " Margin results will be added."))
    if T("t01_datasets"):
        S.append(t_datasets(T("t01_datasets")))
    if alignment:
        S.append(alignment)
    if T("t02_truth"):
        S.append(t_truth(T("t02_truth")))

    S.append("\\clearpage\n\\section{Phasing with the truth as input}\nWith the GIAB truth as input "
             "every phaser receives the same unphased benchmark het SNPs, so only the phaser "
             "differs. %s gives the metrics at 30%s, %s the error types, %s how often each phaser "
             "ranks first and %s the change with coverage.\n"
             % (ref("tab:main"), X, ref("tab:error-types"), ref("tab:ranking"), ref("tab:coverage")))
    for name, make in [("t05_main_30x", t_main), ("t10_error_types", t_error_types),
                       ("t07_ranking", t_ranking), ("t08_coverage", t_coverage)]:
        if T(name):
            S.append(make(T(name)))

    S.append("\\clearpage\n\\section{Variant calling and phasing}\n%s gives the accuracy of the "
             "het SNP calls, %s what calling and what phasing lose on the benchmark het SNPs.\n"
             % (ref("tab:callers"), ref("tab:call-vs-phase")))
    for name, make in [("t03_callers", t_callers), ("t09_call_vs_phase", t_call_vs_phase)]:
        if T(name):
            S.append(make(T(name)))

    S.append("\\clearpage\n\\section{Read-level metrics}\nReads were assigned a haplotype with "
             "\\texttt{whatshap haplotag}, by the VCF of every run and by the phased truth. %s "
             "compares two ways to orient phase blocks before scoring.\n" % ref("tab:reads"))
    if T("t11_reads_30x"):
        S.append(t_reads(T("t11_reads_30x")))

    S.append("\\clearpage\n\\section{GRCh38 against T2T-CHM13}\n%s pools the ratio T2T-CHM13 / GRCh38 "
             "per phaser and technology, %s gives it per platform.\n"
             % (ref("tab:reference"), ref("tab:reference-platform")))
    if T("t06b_reference_pooled"):
        S.append(t_reference_pooled(T("t06b_reference_pooled"), grid))
    if T("t06a_reference"):
        S.append(t_reference_platform(T("t06a_reference")))

    if T("t14_software"):
        S.append("\\clearpage\n\\section{Software}\n" + t_software(T("t14_software")))
    S.append("\\end{document}\n")
    with open(args.out, "w") as fh:
        fh.write("\n".join(S))
    print("-> %s%s" % (args.out, "" if wilcoxon else " (no scipy: p-values left empty)"))


if __name__ == "__main__":
    main()
