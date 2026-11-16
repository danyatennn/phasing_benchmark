#!/usr/bin/env python3
"""Pairwise phasing precision and recall of one run against the phased truth, chr1-22.

Pairs are all pairs of het SNPs, not only neighbouring ones. A SNP of the run
counts when it is a truth het SNP (same position, REF, ALT) phased by the run
(GT with '|'); in a block (PS) with a such SNPs in the truth's orientation and
b in the other, the a(a-1)/2 + b(b-1)/2 pairs inside each group are phased
correctly relative to each other and the a*b pairs across them are not.

  predicted_pairs  sum over blocks of n(n-1)/2, n = a + b: pairs the run phases
  correct_pairs    sum over blocks of a(a-1)/2 + b(b-1)/2
  truth_pairs      sum over chromosomes of N(N-1)/2, N = truth het SNPs: every
                   pair of a chromosome has a known relative phase in the truth
  pairwise_precision_pct  100 * correct_pairs / predicted_pairs
  pairwise_recall_pct     100 * correct_pairs / truth_pairs

  python plots/pairwise.py --truth truth.phased.vcf.gz --vcf run.vcf.gz --name <run>
prints one header line and one row: the fields of the run name, then the above.
"""

import argparse
import gzip

AUTOSOMES = {"chr%d" % i for i in range(1, 23)}


def snps(path, truth):
    """chrom -> {pos: (ref, alt, first allele, ps)} of the phased het SNPs"""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 10)
            if f[0] not in AUTOSOMES or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            fmt = f[8].split(":")
            val = f[9].rstrip("\n").split(":")
            gt = val[0]
            if "|" not in gt or gt not in ("0|1", "1|0"):
                continue
            ps = f[0]
            if not truth and "PS" in fmt:
                i = fmt.index("PS")
                ps = val[i] if i < len(val) and val[i] not in ("", ".") else f[0]
            out.setdefault(f[0], {})[int(f[1])] = (f[3].upper(), f[4].upper(), gt[0], ps)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--truth", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--name", required=True, help="<platform>_<ref>_<sample>_<depth>x_<region>_<caller>_<phaser>")
    a = ap.parse_args()

    truth, run = snps(a.truth, True), snps(a.vcf, False)
    predicted = correct = truth_pairs = 0
    for chrom, t in truth.items():
        n = len(t)
        truth_pairs += n * (n - 1) // 2
        blocks = {}                     # block -> [same orientation as the truth, opposite]
        for pos, x in run.get(chrom, {}).items():
            y = t.get(pos)
            if y and x[0] == y[0] and x[1] == y[1]:
                blocks.setdefault(x[3], [0, 0])[0 if x[2] == y[2] else 1] += 1
        for same, opp in blocks.values():
            predicted += (same + opp) * (same + opp - 1) // 2
            correct += same * (same - 1) // 2 + opp * (opp - 1) // 2

    p = a.name.split("_")
    head = ["platform", "reference", "sample", "depth", "region", "caller", "phaser",
            "truth_pairs", "predicted_pairs", "correct_pairs", "pairwise_precision_pct", "pairwise_recall_pct"]
    row = p[:4] + [p[4], p[5], "_".join(p[6:])] if len(p) >= 7 else [a.name] + [""] * 6
    row[3] = row[3].rstrip("x")
    row += [truth_pairs, predicted, correct,
            "%.4f" % (100.0 * correct / predicted if predicted else 0),
            "%.6f" % (100.0 * correct / truth_pairs if truth_pairs else 0)]
    print("\t".join(head))
    print("\t".join(str(v) for v in row))


if __name__ == "__main__":
    main()
