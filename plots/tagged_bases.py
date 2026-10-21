#!/usr/bin/env python3
"""The share of read bases, not of reads, that a run put on a haplotype.

reads_tagged / reads_total of eval_run.sh counts a 500 bp read like a 50 kb
one, and a read without a het SNP cannot be tagged: ONT, with many short
reads, looks 10-15 points worse than PacBio by reads and is ahead by bases.
Bases are what a haplotype-resolved analysis gets to use.

Population: the primary alignments on chr1-22 of the cell's read table
(bin/read_tables.sh) - the same reads for every run of the cell, whatever
contigs the run's VCF names (the haplotag list of the truth skips chrM and
decoys, those of the callers do not). A read counts as tagged if it has H1/H2
in a block of at least --min-snps phased het SNPs, as in readacc_run.sh.

One line per run (metrics/<run>.tagbases_<rules>.tsv):
  reference caller platform depth phaser reads_total reads_tagged bases_total
  bases_tagged bases_tagged_fraction

  python plots/tagged_bases.py --reads-table <cell>.reads.tsv.gz --min-snps 2 --jobs runs.tsv
runs.tsv: reference caller platform depth phaser tags vcf out
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import haplotag_list, phase_blocks, read_table

TAGBASES_COLUMNS = ["reference", "caller", "platform", "depth", "phaser", "reads_total",
                    "reads_tagged", "bases_total", "bases_tagged", "bases_tagged_fraction"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reads-table", required=True)
    parser.add_argument("--min-snps", type=int, default=2)
    parser.add_argument("--jobs", required=True)
    args = parser.parse_args()

    reads = read_table(args.reads_table, autosomes_only=True)
    lengths = {r: v[3] for r, v in reads.items()}
    del reads
    bases_total = sum(lengths.values())
    print("%s: %d reads, %d bases on chr1-22" % (args.reads_table, len(lengths), bases_total),
          file=sys.stderr)

    with open(args.jobs) as fh:
        jobs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for job in jobs:
        labels, (tags_path, vcf, out) = job[:5], job[5:8]
        sizes = {b: len(p) for b, p in phase_blocks(vcf, autosomes_only=True).items()}
        tags = haplotag_list(tags_path, autosomes_only=True)
        n = bases = 0
        for read, (label, ps, chrom) in tags.items():
            if label not in ("H1", "H2") or sizes.get((chrom, ps), 0) < args.min_snps:
                continue
            length = lengths.get(read)
            if length is None:
                continue
            n += 1
            bases += length
        fields = labels + [len(lengths), n, bases_total, bases,
                           "%.6f" % (bases / bases_total if bases_total else 0)]
        tmp = "%s.%d.tmp" % (out, os.getpid())
        with open(tmp, "w") as fh:
            fh.write("\t".join(str(x) for x in fields) + "\n")
        os.replace(tmp, out)
        print("  %s: %s of the bases tagged" % (os.path.basename(out), fields[-1]), file=sys.stderr)


if __name__ == "__main__":
    main()
