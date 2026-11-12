#!/usr/bin/env python3
"""WhatsHap haplotag against LongPhase haplotag on the same phased VCF and BAM.

Reads: the primary alignments of the BAM on the chromosome (reads.tsv: name,
mapq, length). Both tools tag a read with a haplotype (1/2) of a phase block;
on the same VCF, haplotype 1 is the first allele of the GT for both, so the
tags compare directly. Per read one class:

  none          neither tags it
  wh_only       WhatsHap only
  lp_only       LongPhase only
  same          both, same haplotype and same block
  opposite      both, same block, the other haplotype
  other_block   both, different blocks

  python plots/haplotag_compare.py --reads reads.tsv --whatshap wh.tsv \
      --longphase lp.tsv --chrom chr20 --label <run> --out out.tsv
out.tsv (appended): label class reads pct median_mapq median_length_kb
"""

import argparse
import os
import statistics
from collections import defaultdict

CLASSES = ["none", "wh_only", "lp_only", "same", "opposite", "other_block"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reads", required=True, help="name mapq length, primary alignments")
    ap.add_argument("--whatshap", required=True, help="whatshap --output-haplotag-list")
    ap.add_argument("--longphase", required=True, help="name HP PS of the LongPhase BAM")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    reads = {}
    with open(a.reads) as fh:
        for line in fh:
            n, q, l = line.split()
            reads[n] = (int(q), int(l))
    wh = {}
    with open(a.whatshap) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 4 or f[3] != a.chrom or f[1] not in ("H1", "H2"):
                continue
            wh[f[0]] = (f[1][1], f[2])
    lp = {}
    with open(a.longphase) as fh:
        for line in fh:
            n, hp, ps = line.rstrip("\n").split("\t")
            if hp in ("1", "2"):
                lp[n] = (hp, ps)

    by_class = defaultdict(list)
    for n, (q, l) in reads.items():
        w, p = wh.get(n), lp.get(n)
        if not w and not p:
            c = "none"
        elif not p:
            c = "wh_only"
        elif not w:
            c = "lp_only"
        elif w[1] != p[1]:
            c = "other_block"
        else:
            c = "same" if w[0] == p[0] else "opposite"
        by_class[c].append((q, l))

    new = not os.path.exists(a.out)
    with open(a.out, "a") as out:
        if new:
            out.write("label\tclass\treads\tpct\tmedian_mapq\tmedian_length_kb\n")
        for c in CLASSES:
            v = by_class[c]
            out.write("%s\t%s\t%d\t%.3f\t%s\t%s\n" % (
                a.label, c, len(v), 100.0 * len(v) / max(1, len(reads)),
                statistics.median(q for q, _ in v) if v else "",
                "%.1f" % (statistics.median(l for _, l in v) / 1e3) if v else ""))


if __name__ == "__main__":
    main()
