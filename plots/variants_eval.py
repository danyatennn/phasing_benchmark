#!/usr/bin/env python3
"""Experiment A, what eval_run.sh does not measure: how well the variants added
to the truth SNPs (indels; for truthall also the halves of 1/2 sites and the
parts of complex variants) are phased, and which reads they let haplotag tag.
chr1-22 only.

Phase of an added variant. A phase block's orientation against the truth is
set by the majority of its truth SNPs (the SNPs every input shares); an added
variant in that block is right if its genotype, so oriented, is the truth's.
This works where two records share a position, which whatshap compare does
not (it keeps the first). Per class:
  in_input   records of the class in the truth input
  phased     phased in the run (a phased het genotype with a PS)
  oriented   phased in a block whose truth SNPs decide its orientation
  correct / wrong   of the oriented ones

Reads. Against the same phaser on the SNP-only input (caller truth):
  reads_tagged_base   reads tagged H1/H2 there
  reads_tagged        reads tagged H1/H2 here
  reads_rescued       tagged here, not there
  reads_rescued_no_snp  of those, reads covering no truth het SNP (the allele
                        matrix of the SNP input has no line for them)
  reads_lost          tagged there, not here

  python plots/variants_eval.py --run run.vcf.gz --truth-snps truth.het_snps.phased.vcf.gz \\
      --truth-extra truth.indel.phased.vcf.gz --tags run.haplotag.tsv \\
      --base-tags base.haplotag.tsv --alleles reads/<cell>_truth.reads.tsv.gz \\
      --label "GRCh38 truthindel ont4k 30 whatshap" --out run.variants.tsv
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import haplotag_list, is_autosome, open_text, phased_het, vcf_records

COLUMNS = ["reference", "caller", "platform", "depth", "phaser",
           "indel_in_input", "indel_phased", "indel_oriented", "indel_correct", "indel_wrong",
           "other_in_input", "other_phased", "other_oriented", "other_correct", "other_wrong",
           "reads_tagged_base", "reads_tagged", "reads_rescued", "reads_rescued_no_snp", "reads_lost"]


def truth_gts(path):
    """(chrom, pos, ref, alt) -> GT of the phased truth, chr1-22"""
    out = {}
    for chrom, pos, ref, alt, gt, _ in vcf_records(path, autosomes_only=True):
        if phased_het(gt) is not None:
            out[(chrom, pos, ref, alt)] = gt
    return out


def flip(gt):
    a, b = gt.split("|")
    return b + "|" + a


def phase_scores(run, snps, extra):
    """per class: [in_input, phased, oriented, correct, wrong]"""
    blocks, calls = {}, []
    for chrom, pos, ref, alt, gt, ps in vcf_records(run, autosomes_only=True):
        k = (chrom, pos, ref, alt)
        if phased_het(gt) is None or ps in (".", ""):
            continue
        b = (chrom, ps)
        if k in snps:
            v = blocks.setdefault(b, [0, 0])
            if gt == snps[k]:
                v[0] += 1
            elif gt == flip(snps[k]):
                v[1] += 1
        elif k in extra:
            calls.append((k, b, gt))
    score = {c: [0, 0, 0, 0, 0] for c in ("indel", "other")}
    for k, (cls, _) in extra.items():
        score[cls][0] += 1
    for k, b, gt in calls:
        cls, truth = extra[k]
        s = score[cls]
        s[1] += 1
        same, opposite = blocks.get(b, (0, 0))
        if same == opposite:
            continue                      # no truth SNP, or a tie: orientation unknown
        s[2] += 1
        oriented = gt if same > opposite else flip(gt)
        if oriented == truth:
            s[3] += 1
        else:
            s[4] += 1
    return score


def tagged(path):
    return {r for r, (label, _, _) in haplotag_list(path, autosomes_only=True).items()
            if label in ("H1", "H2")}


def reads_with_snps(path):
    out = set()
    if not path or not os.path.exists(path):
        return None
    with open_text(path) as fh:
        fh.readline()
        for line in fh:
            f = line.split("\t", 3)
            if len(f) > 2 and int(f[2]) > 0:
                out.add(f[0])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--run", required=True)
    parser.add_argument("--truth-snps", required=True)
    parser.add_argument("--truth-extra", required=True, help="the phased truth of this input")
    parser.add_argument("--tags", required=True)
    parser.add_argument("--base-tags", required=True)
    parser.add_argument("--alleles", default="", help="reads.tsv.gz of the SNP input's allele matrix")
    parser.add_argument("--label", required=True, help="reference caller platform depth phaser")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    snps = truth_gts(args.truth_snps)
    extra = {}
    for k, gt in truth_gts(args.truth_extra).items():
        if k in snps:
            continue
        _, _, ref, alt = k
        extra[k] = ("indel" if len(ref) != len(alt) else "other", gt)
    score = phase_scores(args.run, snps, extra)

    here, base = tagged(args.tags), tagged(args.base_tags)
    rescued = here - base
    with_snps = reads_with_snps(args.alleles)
    no_snp = float("nan") if with_snps is None else len(rescued - with_snps)

    values = args.label.split() + score["indel"] + score["other"] + \
        [len(base), len(here), len(rescued), no_snp, len(base - here)]
    tmp = "%s.%d" % (args.out, os.getpid())
    with open(tmp, "w") as fh:
        fh.write("\t".join(str(v) for v in values) + "\n")
    os.replace(tmp, args.out)
    i, o = score["indel"], score["other"]
    print("%s: indels %d/%d phased, %d wrong of %d oriented; other %d/%d phased, %d wrong of %d; "
          "reads +%d -%d" % (args.label, i[1], i[0], i[4], i[2], o[1], o[0], o[4], o[2],
                             len(rescued), len(base - here)), file=sys.stderr)


if __name__ == "__main__":
    main()
