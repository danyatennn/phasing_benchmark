#!/usr/bin/env python3
"""One line per primary alignment of a BAM: where the read lies and what it is
like. The haplotag lists say only which haplotype a read got; this says how
long the read is, how well it maps and how clean it is, for the analyses of
untagged reads, of the anchor read, of tagged bases.

  python plots/read_table.py --bam <cell>/input.bam --out reads/<cell>.reads.tsv.gz

Columns (positions 0-based, half-open, as in BED):
  read chrom start end   the primary alignment
  length                 the whole read, soft and hard clips included
  mapq
  aligned                read bases in the alignment (M, I, =, X)
  softclip               soft-clipped bases
  nm                     NM tag, edit distance of the alignment; -1 if absent
  ins dels               inserted and deleted bases (CIGAR I, D)
  mean_bq                mean base quality of the read; nan without qualities

Secondary, supplementary and unmapped records are skipped: a read appears once.
Needs pysam (the phasing env).
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pysam

from benchio import READ_TABLE_COLUMNS, write_text

try:
    import numpy as np
except ImportError:
    np = None

# indices into AlignedSegment.get_cigar_stats()[0]
M, I, D, S, H, EQ, X, NM = 0, 1, 2, 4, 5, 7, 8, 10


def mean_quality(read):
    q = read.query_qualities
    if q is None or len(q) == 0:
        return float("nan")
    if np is not None:
        return float(np.frombuffer(q, dtype=np.uint8).mean())
    return sum(q) / len(q)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--bam", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--threads", type=int, default=4, help="BAM decompression threads")
    args = parser.parse_args()

    out = write_text(args.out)
    out.write("\t".join(READ_TABLE_COLUMNS) + "\n")
    n = 0
    try:
        with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                bases, _ = read.get_cigar_stats()
                aligned = bases[M] + bases[I] + bases[EQ] + bases[X]
                length = read.infer_read_length() or read.query_length or 0
                nm = read.get_tag("NM") if read.has_tag("NM") else -1
                bq = mean_quality(read)
                out.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" % (
                    read.query_name, read.reference_name, read.reference_start,
                    read.reference_end, length, read.mapping_quality, aligned,
                    bases[S], nm, bases[I], bases[D],
                    "nan" if math.isnan(bq) else "%.2f" % bq))
                n += 1
    except BaseException:
        out.abandon()
        raise
    out.commit()
    print("%s: %d primary alignments -> %s" % (args.bam, n, args.out), file=sys.stderr)


if __name__ == "__main__":
    main()
