#!/usr/bin/env python3
"""FNR and FDR of a caller's heterozygous SNPs, counted where they can be
counted: inside the confident regions and on chr1-22.

The row of eval_run.sh divides by every site the run emitted, inside the BED or
not, and compares with the phased truth only. Here:
  truth   every biallelic het SNP of the GIAB VCF (phased or not) inside the
          BED on chr1-22 - for a sample whose truth is partly unphased this is
          more than the phasing truth
  calls   the caller's het SNPs (the phaser input) inside the BED on chr1-22
  a position with more than one call record (a 1/2 site written as two 0/1
  records, as longcallD does) is not a biallelic het SNP: it leaves both sides
  a call matches the truth on chromosome, position, REF and ALT

One line per cell and caller:
  reference sample platform depth caller truth_het_snps calls shared fn fp fnr
  fdr calls_outside_bed calls_off_autosomes multi_record_positions

  python plots/calls_fdr.py --truth-vcf giab.vcf.gz --bed confident.bed --jobs jobs.tsv
jobs.tsv: reference sample platform depth caller vcf out
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import IntervalSet, is_autosome, open_text

CALLS_COLUMNS = ["reference", "sample", "platform", "depth", "caller", "truth_het_snps", "calls",
                 "shared", "fn", "fp", "fnr", "fdr", "calls_outside_bed", "calls_off_autosomes",
                 "multi_record_positions"]


def het_biallelic_snps(path, bed, require_het):
    """{(chrom, pos): [(ref, alt), ...]} inside the BED on chr1-22, every record
    kept so that positions with several can be found; plus the counts of
    records outside the BED and off the autosomes."""
    sites = {}
    outside = off = 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 10)
            if len(f) < 10:
                continue
            chrom, pos, ref, alt = f[0], int(f[1]), f[3].upper(), f[4].upper()
            if len(ref) != 1 or len(alt) != 1 or alt in (".", "*"):
                continue
            if require_het:
                keys = f[8].split(":")
                if "GT" not in keys:
                    continue
                gt = f[9].split(":")[keys.index("GT")].replace("|", "/")
                a = gt.split("/")
                if len(a) != 2 or a[0] == a[1] or not set(a) <= {"0", "1"}:
                    continue
            if not is_autosome(chrom):
                off += 1
                continue
            if not bed.contains(chrom, pos - 1):
                outside += 1
                continue
            sites.setdefault((chrom, pos), []).append((ref, alt))
    return sites, outside, off


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--truth-vcf", required=True, help="the GIAB benchmark VCF")
    parser.add_argument("--bed", required=True, help="its confident regions")
    parser.add_argument("--jobs", required=True, help="TSV: reference sample platform depth caller vcf out")
    args = parser.parse_args()

    bed = IntervalSet.from_bed(args.bed, autosomes_only=True)
    truth, _, _ = het_biallelic_snps(args.truth_vcf, bed, require_het=True)
    # the truth itself can hold two records at a position; such a site is not
    # a biallelic het SNP either
    truth = {k: v[0] for k, v in truth.items() if len(v) == 1}
    print("%s: %d het biallelic SNPs on chr1-22 inside the BED" % (args.truth_vcf, len(truth)),
          file=sys.stderr)

    with open(args.jobs) as fh:
        jobs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for job in jobs:
        labels, vcf, out = job[:5], job[5], job[6]
        calls, outside, off = het_biallelic_snps(vcf, bed, require_het=False)
        multi = {k for k, v in calls.items() if len(v) > 1}
        calls = {k: v[0] for k, v in calls.items() if k not in multi}
        truth_n = sum(1 for k in truth if k not in multi)
        shared = sum(1 for k, ra in calls.items() if truth.get(k) == ra)
        fn = truth_n - shared
        fp = len(calls) - shared
        fields = labels + [truth_n, len(calls), shared, fn, fp,
                           "%.6f" % (fn / truth_n if truth_n else 0),
                           "%.6f" % (fp / len(calls) if calls else 0),
                           outside, off, len(multi)]
        tmp = "%s.%d.tmp" % (out, os.getpid())
        with open(tmp, "w") as fh:
            fh.write("\t".join(str(x) for x in fields) + "\n")
        os.replace(tmp, out)
        print("%s: FNR %s FDR %s (%d outside the BED, %d off chr1-22, %d positions with several records)" % (
            os.path.basename(out), fields[10], fields[11], outside, off, len(multi)), file=sys.stderr)


if __name__ == "__main__":
    main()
