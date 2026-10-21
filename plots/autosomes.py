#!/usr/bin/env python3
"""The phasing metrics of a run on chr1-22 only, from what eval_run.sh left on
disk - nothing is phased or haplotagged again:

  switches, assessed, hamming, covered   compare.tsv, summed over chr1-22
  variants, phased, blocks               stats.tsv, summed over chr1-22
  truth_total                            phased truth SNPs on chr1-22
  emitted, shared                        the run's VCF records on chr1-22, and
                                         those matching a truth SNP (chrom, pos,
                                         REF, ALT), as bcftools isec -c none
  reads_total, reads_tagged              haplotag list, reads on chr1-22
  reads_tagged_min                       of those, tagged in a block of at least
                                         --min-snps phased het SNPs

Rates as in eval_run.sh: switch_rate = switches / assessed, hamming_rate =
hamming / covered, correct = covered - hamming, recall_phased = correct /
truth_total, fnr = 1 - shared / truth_total, fdr = 1 - shared / emitted.

  python plots/autosomes.py --truth-vcf truth.phased.vcf.gz --min-snps 2 --jobs runs.tsv --threads 8
runs.tsv: reference region caller platform depth phaser compare stats vcf tags out
"""

import argparse
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import is_autosome, open_text, phased_het, read_tsv_header, vcf_records

AUTO_COLUMNS = ["reference", "region", "caller", "platform", "depth", "phaser",
                "variants", "phased", "blocks", "switches", "assessed", "hamming", "covered",
                "truth_total", "emitted", "shared", "correct", "switch_rate", "hamming_rate",
                "recall_phased", "fnr", "fdr", "reads_total", "reads_tagged", "tagging_completeness",
                "reads_tagged_min", "tagging_completeness_min"]

TRUTH = None
MIN_SNPS = 2


def init(truth, min_snps):
    global TRUTH, MIN_SNPS
    TRUTH = truth
    MIN_SNPS = min_snps


def chrom_of(row):
    return row.get("chromosome", row.get("chrom", ""))


def summed(rows, column):
    total = 0
    for r in rows:
        if is_autosome(chrom_of(r)):
            v = r.get(column, "")
            if v not in ("", "nan", "NA"):
                total += int(float(v))
    return total


def one(job):
    labels, (cmp_path, stats_path, vcf, tags, out) = job[:6], job[6:11]
    cmp_rows = read_tsv_header(cmp_path)
    stats_rows = read_tsv_header(stats_path)
    switches = summed(cmp_rows, "all_switches")
    assessed = summed(cmp_rows, "all_assessed_pairs")
    hamming = summed(cmp_rows, "blockwise_hamming")
    covered = summed(cmp_rows, "covered_variants")
    variants = summed(stats_rows, "variants")
    phased = summed(stats_rows, "phased")
    blocks = summed(stats_rows, "blocks")

    emitted = shared = 0
    size = {}
    for chrom, pos, ref, alt, gt, ps in vcf_records(vcf, autosomes_only=True):
        emitted += 1
        if (chrom, pos, ref, alt) in TRUTH:
            shared += 1
        if ps not in (".", "") and phased_het(gt) is not None:
            size[(chrom, ps)] = size.get((chrom, ps), 0) + 1

    seen = set()
    total = tagged = tagged_min = 0
    with open_text(tags) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] in seen:
                continue
            seen.add(f[0])
            if not is_autosome(f[3]):
                continue
            total += 1
            if f[1] in ("H1", "H2"):
                tagged += 1
                if size.get((f[3], f[2]), 0) >= MIN_SNPS:
                    tagged_min += 1

    truth_total = len(TRUTH)
    correct = covered - hamming
    fields = labels + [variants, phased, blocks, switches, assessed, hamming, covered,
                       truth_total, emitted, shared, correct,
                       "%.6f" % (switches / assessed if assessed else 0),
                       "%.6f" % (hamming / covered if covered else 0),
                       "%.6f" % (correct / truth_total if truth_total else 0),
                       "%.6f" % (1 - shared / truth_total if truth_total else 0),
                       "%.6f" % (1 - shared / emitted if emitted else 0),
                       total, tagged, "%.6f" % (tagged / total if total else 0),
                       tagged_min, "%.6f" % (tagged_min / total if total else 0)]
    tmp = "%s.%d.tmp" % (out, os.getpid())
    with open(tmp, "w") as fh:
        fh.write("\t".join(str(x) for x in fields) + "\n")
    os.replace(tmp, out)
    return os.path.basename(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--truth-vcf", required=True, help="the phased truth the runs were compared with")
    parser.add_argument("--min-snps", type=int, default=2)
    parser.add_argument("--jobs", required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    truth = {(c, p, r, a) for c, p, r, a, gt, ps in vcf_records(args.truth_vcf, autosomes_only=True)}
    print("%s: %d SNPs on chr1-22" % (args.truth_vcf, len(truth)), file=sys.stderr)
    with open(args.jobs) as fh:
        jobs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    with Pool(max(1, args.threads), initializer=init, initargs=(truth, args.min_snps)) as pool:
        for name in pool.imap_unordered(one, jobs):
            print("  %s" % name, file=sys.stderr)


if __name__ == "__main__":
    main()
