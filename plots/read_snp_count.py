#!/usr/bin/env python3
"""Het SNPs of a VCF that each read covers with REF or ALT (n_ref + n_alt of
the truth allele matrix of the grid), from SAM records on stdin.

  samtools view reads.bam | python plots/read_snp_count.py --vcf truth.vcf.gz > read_snps.tsv
-> read  snps   (reads with none are left out)
"""

import argparse
import bisect
import gzip
import re
import sys
from collections import defaultdict

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf", required=True)
    a = ap.parse_args()
    pos = defaultdict(list)          # chrom -> sorted 0-based positions
    alleles = {}
    with gzip.open(a.vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 6)
            if len(f[3]) == 1 and len(f[4]) == 1:
                p = int(f[1]) - 1
                pos[f[0]].append(p)
                alleles[(f[0], p)] = (f[3].upper(), f[4].upper())
    for v in pos.values():
        v.sort()
    out = sys.stdout
    out.write("read\tsnps\n")
    for line in sys.stdin:
        f = line.split("\t", 11)
        chrom, start, cigar, seq = f[2], int(f[3]) - 1, f[5], f[9]
        sites = pos.get(chrom)
        if not sites or cigar == "*" or seq == "*":
            continue
        span = sum(int(n) for n, op in CIGAR.findall(cigar) if op in "MDN=X")
        i = bisect.bisect_left(sites, start)
        j = bisect.bisect_left(sites, start + span)
        if i == j:
            continue
        want = sites[i:j]
        k, n, rp, qp = 0, 0, start, 0
        for length, op in CIGAR.findall(cigar):
            length = int(length)
            if op in "M=X":
                while k < len(want) and want[k] < rp + length:
                    if want[k] >= rp:
                        b = seq[qp + want[k] - rp].upper()
                        if b in alleles[(chrom, want[k])]:
                            n += 1
                    k += 1
                rp += length
                qp += length
            elif op in "IS":
                qp += length
            elif op in "DN":
                while k < len(want) and want[k] < rp + length:
                    k += 1
                rp += length
            if k >= len(want):
                break
        if n:
            out.write("%s\t%d\n" % (f[0], n))


if __name__ == "__main__":
    main()
