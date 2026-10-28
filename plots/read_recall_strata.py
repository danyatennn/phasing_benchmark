#!/usr/bin/env python3
"""Read-level haplotag recall per genomic region, for every run of one cell.

  recall    = reads on their truth haplotype / reads the truth tags
  precision = reads on their truth haplotype / reads tagged by both

Reads: primary alignments on chr1-22 (the cell's read table) that whatshap
haplotag puts on H1/H2 with the phased truth (<cell>_truthtags.tsv) - the
denominator of recall. Of those, a read is right if the run tags it H1/H2 in a
block of --min-snps phased het SNPs or more, the block shares SNPs with the
truth (which orient it, majority of them), and so oriented its label is the
truth's - the rules of readacc_run.sh, whose genome-wide read_recall this
reproduces. A read is in a region if half its alignment or more lies there;
"genome" is every read.

One line per run x region (metrics_reads/<run>.recall_strata.tsv):
  reference sample caller platform depth phaser stratum truth_tagged
  tagged compared correct unoriented recall precision
    tagged      of the truth-tagged reads, those the run tags (blocks >= min)
    compared    of those, in a block the truth orients
    unoriented  tagged, in a block that shares no SNP with the truth

  python plots/read_recall_strata.py --reads-table <cell>.reads.tsv.gz \\
      --truth-tags <cell>_truthtags.tsv --truth-vcf truth.phased.vcf.gz \\
      --stratum segdups=prepared_auto/GRCh38_segdups.bed ... --min-snps 2 --jobs runs.tsv
runs.tsv: reference sample caller platform depth phaser tags vcf out
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import (IntervalSet, haplotag_list, is_autosome, open_text, phased_het, truth_labels,
                     vcf_records)
from read_classes import blocks_and_orientation

STRATUM_FRACTION = 0.5
COLUMNS = ["reference", "sample", "caller", "platform", "depth", "phaser", "stratum", "truth_tagged",
           "tagged", "compared", "correct", "unoriented", "recall", "precision"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reads-table", required=True)
    ap.add_argument("--truth-tags", required=True)
    ap.add_argument("--truth-vcf", required=True)
    ap.add_argument("--stratum", action="append", default=[], help="name=bed, repeatable")
    ap.add_argument("--min-snps", type=int, default=2)
    ap.add_argument("--jobs", required=True)
    args = ap.parse_args()

    truth = truth_labels(args.truth_tags, "truthtags", autosomes_only=True)
    strata = []
    for item in args.stratum:
        name, _, path = item.partition("=")
        if os.path.exists(path):
            strata.append((name, IntervalSet.from_bed(path, autosomes_only=True)))
    # the regions of every truth-tagged read, from its primary alignment
    region_of = {}
    with open_text(args.reads_table) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ic, i_s, ie = h.index("chrom"), h.index("start"), h.index("end")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] not in truth or not is_autosome(f[ic]):
                continue
            c, s, e = f[ic], int(f[i_s]), int(f[ie])
            span = max(1, e - s)
            region_of[f[0]] = ("genome",) + tuple(
                n for n, iv in strata if iv.overlap(c, s, e) >= STRATUM_FRACTION * span)
    # a truth-tagged read missing from the read table (not a primary alignment
    # on chr1-22 there) still counts genome-wide, as in readacc_run.sh
    for r in truth:
        region_of.setdefault(r, ("genome",))
    print("%s: %d truth-tagged reads on chr1-22" % (args.reads_table, len(truth)), file=sys.stderr)

    truth_gt = {(c, p, ref, alt): gt for c, p, ref, alt, gt, ps in
                vcf_records(args.truth_vcf, autosomes_only=True) if phased_het(gt) is not None}
    order = ["genome"] + [n for n, _ in strata]
    with open(args.jobs) as fh:
        jobs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for job in jobs:
        labels, (tags_path, vcf, out) = job[:6], job[6:9]
        size, flip = blocks_and_orientation(vcf, truth_gt)
        tags = haplotag_list(tags_path, autosomes_only=True)
        acc = {n: [0, 0, 0, 0, 0] for n in order}     # truth_tagged tagged compared correct unoriented
        for r, lab in truth.items():
            regions = region_of[r]
            t = tags.get(r)
            tagged = t is not None and t[0] in ("H1", "H2") and size.get((t[2], t[1]), 0) >= args.min_snps
            comp = corr = unor = 0
            if tagged:
                f = flip.get((t[2], t[1]))
                if f is None:
                    unor = 1
                else:
                    comp = 1
                    corr = int(((0 if t[0] == "H1" else 1) ^ f) == lab)
            for n in regions:
                a = acc[n]
                a[0] += 1
                a[1] += tagged
                a[2] += comp
                a[3] += corr
                a[4] += unor
        lines = []
        for n in order:
            tt, tg, cp, cr, un = acc[n]
            if not tt:
                continue
            lines.append("\t".join(labels + [n, str(tt), str(tg), str(cp), str(cr), str(un),
                                             "%.6f" % (cr / tt), "%.6f" % (cr / cp if cp else 0)]) + "\n")
        tmp = "%s.%d.tmp" % (out, os.getpid())
        with open(tmp, "w") as fh:
            fh.write("\t".join(COLUMNS) + "\n")
            fh.writelines(lines)
        os.replace(tmp, out)
        g = acc["genome"]
        print("  %s: genome recall %.4f (%d of %d)" % (labels[5], g[3] / g[0] if g[0] else 0, g[3], g[0]),
              file=sys.stderr)


if __name__ == "__main__":
    main()
