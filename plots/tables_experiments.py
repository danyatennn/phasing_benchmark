#!/usr/bin/env python3
"""Tables of the experiments beyond the main grid, the numbers behind
plots/figures_experiments.py (the same loader and value functions). Rows appear
as runs finish; a table with nothing yet is skipped. Each table is a TSV in
--out, all of them together in experiments.xlsx, and README_experiments.tsv
says what each holds.

  tA1_inputs             experiment A, every input, sample, platform and phaser at 30x
  tA2_inputs_change      truthindel and truthall against the truth input of the same run
  tA3_added_variants     how the added indels / other variants are phased, reads gained and lost
  tB1_methylation        each methylation-aware run against its phaser alone
  tC1_assemblies         every assembly run: phasing and read-level metrics
  tC2_assemblies_vs_reads  assemblies and read-based phasers (truth and DeepVariant input) per cell

  python plots/tables_experiments.py
  python plots/tables_experiments.py --samples "HG002 HG001" --out tables
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import style as st
from data import ANCHOR_TRUTHS, NAN, Metrics, ratio, sample_root
from figures import (n50_mb, read_precision_pct, read_recall_pct, switch_pct, tagged_pct,
                     truth_phased_pct)
from figures_experiments import (ASM_DEPTHS, ASSEMBLERS, DEPTH, INPUTS_A, hamming_pct,
                                 phasers_on, samples_of, variants_table)
from make_tables import (Table, caller_label, finite, phaser_label, plat, ref_label, write_tsv,
                         write_xlsx)

ASM_LABEL = {"hifiasm": "hifiasm (no parents)", "hifiasmtrio": "hifiasm trio",
             "verkkotrio": "Verkko trio"}


def anchor_cells(db, key):
    """precision, recall against each read truth (blocks oriented by an anchor read)"""
    cells = []
    for truth in ANCHOR_TRUTHS:
        cells += [100 * db.get(key, "read_precision_anchor_%s" % truth),
                  100 * db.get(key, "read_recall_anchor_%s" % truth)]
    return cells


ANCHOR_HEADER = ["%s_%s_pct" % (m, truth) for truth in ANCHOR_TRUTHS for m in ("precision", "recall")]


def f1(p, r):
    return 2 * p * r / (p + r) if finite(p) and finite(r) and p + r > 0 else NAN


def run_cells(db, key):
    """The values every run table shares."""
    return [switch_pct(db, key), hamming_pct(db, key), n50_mb(db, key), truth_phased_pct(db, key),
            tagged_pct(db, key), read_precision_pct(db, key), read_recall_pct(db, key)]


RUN_HEADER = ["switch_rate_pct", "hamming_rate_pct", "ng50_mb", "truth_phased_pct",
              "reads_tagged_pct", "read_precision_pct", "read_recall_pct"]


# -------------------------------------------------------------------- A

def tA1(dbs):
    t = Table("tA1_inputs", "Experiment A: phasing with extra variants in the input, %dx" % DEPTH,
              "The same phaser on the same cell with three inputs: truth = GIAB truth het SNPs; "
              "truthindel = the same SNPs + het indels; truthall = + every other het variant (1/2 "
              "sites split, complex variants atomized; bin/make_truth_variants.sh). Phasing metrics "
              "are over the SNPs only (whatshap compare --only-snvs), chr1-22, blocks of >= 2 SNPs. "
              "Read precision/recall: blocks oriented by their SNPs; *_truthtags / matpat / yak: "
              "blocks oriented by an anchor read, against each read truth.",
              ["sample", "reference", "platform", "phaser", "input"] + RUN_HEADER + ANCHOR_HEADER)
    for sample, (db, refs) in dbs.items():
        for ref in refs:
            for p in st.PLATFORMS:
                for ph in phasers_on(p):
                    if not any(db.has((ref, c, p, DEPTH, ph)) for c in INPUTS_A[1:]):
                        continue
                    for c in INPUTS_A:
                        k = (ref, c, p, DEPTH, ph)
                        if db.has(k):
                            t.add(sample, ref_label(ref), plat(p), phaser_label(ph), caller_label(c),
                                  *(run_cells(db, k) + anchor_cells(db, k)))
    return t


def tA2(dbs):
    t = Table("tA2_inputs_change", "Experiment A: change against the SNP input, %dx" % DEPTH,
              "Each truthindel / truthall run against the truth-input run of the same sample, "
              "reference, platform and phaser. Ratios are extra / SNP input; differences are in "
              "percentage points.",
              ["sample", "reference", "platform", "phaser", "input", "ng50_ratio",
               "switch_rate_diff_pp", "hamming_rate_ratio", "truth_phased_diff_pp",
               "read_recall_diff_pp"])
    for sample, (db, refs) in dbs.items():
        for ref in refs:
            for p in st.PLATFORMS:
                for ph in phasers_on(p):
                    base = (ref, "truth", p, DEPTH, ph)
                    if not db.has(base):
                        continue
                    for c in INPUTS_A[1:]:
                        k = (ref, c, p, DEPTH, ph)
                        if not db.has(k):
                            continue
                        t.add(sample, ref_label(ref), plat(p), phaser_label(ph), caller_label(c),
                              ratio(n50_mb(db, k), n50_mb(db, base)),
                              switch_pct(db, k) - switch_pct(db, base),
                              ratio(hamming_pct(db, k), hamming_pct(db, base)),
                              truth_phased_pct(db, k) - truth_phased_pct(db, base),
                              read_recall_pct(db, k) - read_recall_pct(db, base))
    return t


def tA3(root, samples):
    t = Table("tA3_added_variants", "Experiment A: the added variants and the reads they tag, %dx" % DEPTH,
              "bin/variants_eval.sh. A block is oriented by the majority of its truth SNPs; an "
              "added variant in it is right if its genotype, so oriented, is the truth's. indel = "
              "length change; other = the remaining added records (truthall: halves of 1/2 sites, "
              "parts of complex variants). Reads against the same phaser on the SNP input: gained "
              "= tagged here, not there; of those, on no truth het SNP; lost = tagged there, not here.",
              ["sample", "reference", "platform", "phaser", "input",
               "indels_in_input", "indels_phased_pct", "indels_wrong_pct_of_oriented",
               "other_in_input", "other_phased_pct", "other_wrong_pct_of_oriented",
               "reads_tagged_snp_input", "reads_tagged", "reads_gained", "reads_gained_no_snp",
               "reads_lost", "reads_net_change_pct"])
    for sample in samples:
        table = variants_table(root, sample)
        order = lambda item: (st.REFERENCES.index(item[0][0]) if item[0][0] in st.REFERENCES else 9,
                              item[0][1], st.PLATFORMS.index(item[0][2]) if item[0][2] in st.PLATFORMS else 9,
                              st.PHASERS.index(item[0][4]) if item[0][4] in st.PHASERS else 9)
        for (ref, c, p, depth, ph), r in sorted(table.items(), key=order):
            t.add(sample, ref_label(ref), plat(p), phaser_label(ph), caller_label(c),
                  r.get("indel_in_input"), 100 * ratio(r.get("indel_phased"), r.get("indel_in_input")),
                  100 * ratio(r.get("indel_wrong"), r.get("indel_oriented")),
                  r.get("other_in_input"), 100 * ratio(r.get("other_phased"), r.get("other_in_input")),
                  100 * ratio(r.get("other_wrong"), r.get("other_oriented")),
                  r.get("reads_tagged_base"), r.get("reads_tagged"), r.get("reads_rescued"),
                  r.get("reads_rescued_no_snp"), r.get("reads_lost"),
                  100 * ratio(r.get("reads_tagged", NAN) - r.get("reads_tagged_base", NAN),
                              r.get("reads_tagged_base")))
    return t


# -------------------------------------------------------------------- B

def tB1(dbs):
    t = Table("tB1_methylation", "Experiment B: methylation-aware phasing, truth input, %dx" % DEPTH,
              "lpmod = LongPhase with its 5mC mode (longphase modcall, phase --mod-file), against "
              "LongPhase; mp<phaser> = MethPhaser on top of that phaser's VCF, against the phaser "
              "alone. Only cells whose reads carry 5mC calls. MethPhaser does not write unphased "
              "records before the first and after the last blocks of a chromosome, so compare "
              "phased counts, not phased fractions. Ratios are method / phaser alone.",
              ["sample", "reference", "platform", "method", "phaser_alone"]
              + ["method_" + h for h in RUN_HEADER] + ["alone_" + h for h in RUN_HEADER]
              + ["ng50_ratio", "hamming_rate_ratio", "switch_rate_diff_pp"])
    for sample, (db, refs) in dbs.items():
        for ref in refs:
            for p in st.PLATFORMS:
                for m in st.METH_PHASERS:
                    if not st.applies(p, m):
                        continue
                    base_phaser = "longphase" if m == "lpmod" else m[2:]
                    k = (ref, "truth", p, DEPTH, m)
                    base = (ref, "truth", p, DEPTH, base_phaser)
                    if not db.has(k):
                        continue
                    t.add(sample, ref_label(ref), plat(p), phaser_label(m), phaser_label(base_phaser),
                          *(run_cells(db, k) + run_cells(db, base)),
                          ratio(n50_mb(db, k), n50_mb(db, base)),
                          ratio(hamming_pct(db, k), hamming_pct(db, base)),
                          switch_pct(db, k) - switch_pct(db, base))
    return t


# -------------------------------------------------------------------- C

def tC1(dbs):
    t = Table("tC1_assemblies", "Assemblies through dipcall, every run",
              "Phase sets: hifiasm without parents - a run of SNPs on one pair of contigs; trio - a "
              "chromosome. SNPs: the assembly's own (dipcall), biallelic het, chr1-22, inside the "
              "confident BED and dip.bed; snps_compared = SNP pairs whatshap compare assessed. Read "
              "metrics: reads tagged with the assembly's VCF, blocks oriented by their SNPs (and by "
              "an anchor read against each read truth).",
              ["sample", "reference", "assembler", "platform", "depth"] + RUN_HEADER
              + ["read_f1_pct", "snps_compared", "blocks_oriented", "blocks_flipped"] + ANCHOR_HEADER)
    for sample, (db, refs) in dbs.items():
        for ref in refs:
            for a in ASSEMBLERS:
                for p in st.PLATFORMS:
                    for d in ASM_DEPTHS:
                        k = (ref, a, p, d, "dipcall")
                        if not db.has(k):
                            continue
                        t.add(sample, ref_label(ref), ASM_LABEL[a], plat(p), d,
                              *(run_cells(db, k)
                                + [f1(read_precision_pct(db, k), read_recall_pct(db, k)),
                                   db.get(k, "assessed"), db.get(k, "blocks_oriented"),
                                   db.get(k, "blocks_flipped")]
                                + anchor_cells(db, k)))
    return t


def tC2(dbs):
    t = Table("tC2_assemblies_vs_reads", "Assemblies against read-based phasers, per cell",
              "For every platform and coverage with an assembly: each assembly (its own SNPs) and "
              "each read-based phaser on the truth input and on DeepVariant's calls. The SNP sets "
              "differ: compare snps_compared before reading small differences.",
              ["sample", "reference", "platform", "depth", "method", "input"] + RUN_HEADER
              + ["read_f1_pct", "snps_compared"])
    for sample, (db, refs) in dbs.items():
        for ref in refs:
            for p in st.PLATFORMS:
                for d in ASM_DEPTHS:
                    asm = [a for a in ASSEMBLERS if db.has((ref, a, p, d, "dipcall"))]
                    if not asm:
                        continue
                    rows = [((ref, a, p, d, "dipcall"), ASM_LABEL[a], "own SNPs (dipcall)") for a in asm]
                    rows += [((ref, c, p, d, ph), phaser_label(ph), caller_label(c))
                             for c in ("truth", "deepvariant") for ph in phasers_on(p)]
                    for k, method, inp in rows:
                        if db.has(k):
                            t.add(sample, ref_label(ref), plat(p), d, method, inp,
                                  *(run_cells(db, k)
                                    + [f1(read_precision_pct(db, k), read_recall_pct(db, k)),
                                       db.get(k, "assessed")]))
    return t


# ------------------------------------------------------------------ main

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--out", default=None, help="default <root>/tables")
    parser.add_argument("--samples", default="", help="default: HG002 and samples/*")
    parser.add_argument("--refs", default=" ".join(st.REFERENCES),
                        help="references of HG002; the other samples are GRCh38 only")
    parser.add_argument("--rules", default="auto_min2")
    ctx = parser.parse_args()
    ctx.out = ctx.out or os.path.join(ctx.root, "tables")
    os.makedirs(ctx.out, exist_ok=True)

    samples = [s for s in samples_of(ctx) if os.path.isdir(sample_root(ctx.root, s))]
    refs_of = {s: ctx.refs.split() if s == "HG002" else [r for r in ctx.refs.split() if r == "GRCh38"]
               for s in samples}
    ab = {s: (Metrics(ctx.root, s, st.REGION, readacc_rules=ctx.rules).load(
        refs_of[s], INPUTS_A, st.PLATFORMS, [DEPTH], st.PHASERS + st.METH_PHASERS), refs_of[s])
          for s in samples}
    c = {s: (Metrics(ctx.root, s, st.REGION, readacc_rules=ctx.rules).load(
        refs_of[s], ["truth", "deepvariant"] + ASSEMBLERS, st.PLATFORMS, ASM_DEPTHS,
        st.PHASERS + ["dipcall"]), refs_of[s]) for s in samples}

    tables = [tA1(ab), tA2(ab), tA3(ctx.root, samples), tB1(ab), tC1(c), tC2(c)]
    written = [t for t in tables if t.rows]
    for t in tables:
        print("  %-8s %-26s %d rows" % ("ok" if t.rows else "empty", t.name, len(t.rows)))
    if not written:
        print("nothing to write yet")
        return 0
    readme = Table("README_experiments", "What each table holds", "", ["table", "title", "definition"])
    for t in written:
        readme.add(t.name, t.title, t.note)
    for t in [readme] + written:
        write_tsv(t, ctx.out)
    write_xlsx(os.path.join(ctx.out, "experiments.xlsx"),
               [(t.name, t.header, t.rows) for t in [readme] + written])
    print("%d tables -> %s (TSV each, and experiments.xlsx)" % (len(written), ctx.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
