#!/usr/bin/env python3
"""Read-level haplotag recall per genomic region, for every run of one cell,
counted in reads and in read bases.

  recall    = reads on their truth haplotype / reads the truth tags
  precision = reads on their truth haplotype / reads tagged by both
  recall_bases = the same as recall, each read weighted by its length

Reads: primary alignments on chr1-22 (the cell's read table) that whatshap
haplotag puts on H1/H2 with the phased truth (<cell>_truthtags.tsv) and that
cover --min-read-snps truth het SNPs or more with REF or ALT (the truth allele
matrix, reads/<cell>_truth.reads.tsv.gz: n_ref + n_alt) - the denominator of
recall, the same for every run of the cell. Of those, a read is right if the
run tags it H1/H2 in a block of --min-snps phased het SNPs or more, the block
shares SNPs with the truth (which orient it, majority of them), and so oriented
its label is the truth's - the rules of readacc_run.sh. A read is in a region
if half its alignment or more lies there; "genome" is every read. Bases: the
whole read (column length of the read table), as plots/tagged_bases.py.

Truth (--truth-source):
  truthtags  whatshap haplotag against the phased truth VCF of the run's
             reference (<cell>_truthtags.tsv)
  matpat     the parental copy of the HG002 Q100 assembly a read fits better
             (bin/matpat_truth.sh, reads/<cell>.matpat.tsv.gz), PAT -> 0
  yak        parental k-mers (bin/yak_truth.sh, reads/<cell>.yak.tsv.gz)
matpat and yak do not depend on a reference or a VCF: given the same reads,
the read table, allele matrix and strata of one reference (T2T) for runs on
both, every run is scored against the same reads, labels and regions.
--require-reads <table> keeps only reads with a primary alignment in every such
read table too (the cell on each reference), in case downsampling kept a few
reads on one reference only.

Block orientation (--orient): the label H1/H2 of a block is arbitrary, so it is
matched to the truth first.
  snps   majority of the block's SNPs against the truth VCF (readacc_run.sh)
  reads  majority of the other truth-labelled reads tagged in the block (the
         read itself left out, so it cannot vote for itself); a tie, or no
         other labelled read, leaves the read unoriented. Uses no VCF.

One line per run x region (metrics_reads/<run>.recall_strata.tsv):
  reference sample caller platform depth phaser stratum truth_tagged
  tagged compared correct unoriented recall precision
  truth_bases tagged_bases compared_bases correct_bases recall_bases
  truth_source orient
    tagged      of the truth-tagged reads, those the run tags (blocks >= min)
    compared    of those, in a block the truth orients
    unoriented  tagged, in a block that shares no SNP with the truth

  python plots/read_recall_strata.py --reads-table <cell>.reads.tsv.gz \\
      --truth-alleles reads/<cell>_truth.reads.tsv.gz --truth-tags <cell>_truthtags.tsv \\
      --truth-source truthtags --truth-vcf truth.phased.vcf.gz --stratum segdups=prepared_auto/GRCh38_segdups.bed ... \\
      --min-snps 2 --min-read-snps 2 --jobs runs.tsv
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
           "tagged", "compared", "correct", "unoriented", "recall", "precision",
           "truth_bases", "tagged_bases", "compared_bases", "correct_bases", "recall_bases",
           "truth_source", "orient"]


def snps_per_read(path):
    """read -> truth het SNPs it covers with REF or ALT"""
    out = {}
    with open_text(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        i_ref, i_alt = h.index("n_ref"), h.index("n_alt")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = int(f[i_ref]) + int(f[i_alt])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reads-table", required=True)
    ap.add_argument("--truth-alleles", required=True, help="reads/<cell>_truth.reads.tsv.gz")
    ap.add_argument("--truth-tags", required=True, help="the labels, as --truth-source says")
    ap.add_argument("--truth-source", choices=["truthtags", "matpat", "yak"], default="truthtags")
    ap.add_argument("--orient", choices=["snps", "reads"], default=None,
                    help="default: snps for truthtags, reads for matpat and yak")
    ap.add_argument("--truth-vcf", required=True, help="read with --orient snps only")
    ap.add_argument("--stratum", action="append", default=[], help="name=bed, repeatable")
    ap.add_argument("--min-snps", type=int, default=2, help="phased het SNPs of the run's block")
    ap.add_argument("--min-read-snps", type=int, default=2, help="truth het SNPs on the read")
    ap.add_argument("--require-reads", action="append", default=[],
                    help="read table; keep only reads in it, repeatable")
    ap.add_argument("--jobs", required=True)
    args = ap.parse_args()

    orient = args.orient or ("snps" if args.truth_source == "truthtags" else "reads")
    if args.truth_source == "truthtags":
        truth = truth_labels(args.truth_tags, "truthtags", autosomes_only=True)
    else:
        # chr1-22 by the read table's alignment, below
        truth = truth_labels(args.truth_tags, args.truth_source, autosomes_only=False)
    snps = snps_per_read(args.truth_alleles)
    keep = None
    for path in args.require_reads:
        with open_text(path) as fh:
            fh.readline()
            names = {line.split("\t", 1)[0] for line in fh}
        keep = names if keep is None else keep & names
    strata = []
    for item in args.stratum:
        name, _, path = item.partition("=")
        if os.path.exists(path):
            strata.append((name, IntervalSet.from_bed(path, autosomes_only=True)))
    # the regions and length of every truth-tagged read with enough truth SNPs,
    # from its primary alignment; a read missing from the read table or the
    # allele matrix is left out
    region_of, length_of, autosomal = {}, {}, set()
    with open_text(args.reads_table) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ic, i_s, ie, il = h.index("chrom"), h.index("start"), h.index("end"), h.index("length")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            r = f[0]
            if r not in truth or not is_autosome(f[ic]) or (keep is not None and r not in keep):
                continue
            autosomal.add(r)
            if snps.get(r, 0) < args.min_read_snps:
                continue
            c, s, e = f[ic], int(f[i_s]), int(f[ie])
            span = max(1, e - s)
            region_of[r] = ("genome",) + tuple(
                n for n, iv in strata if iv.overlap(c, s, e) >= STRATUM_FRACTION * span)
            length_of[r] = int(f[il])
    print("%s: %s labels, %d reads labelled on chr1-22, %d with >= %d truth het SNPs; blocks oriented by %s"
          % (args.reads_table, args.truth_source, len(autosomal), len(region_of), args.min_read_snps, orient),
          file=sys.stderr)
    del snps, keep

    truth_gt = {}
    if orient == "snps":
        truth_gt = {(c, p, ref, alt): gt for c, p, ref, alt, gt, ps in
                    vcf_records(args.truth_vcf, autosomes_only=True) if phased_het(gt) is not None}
    order = ["genome"] + [n for n, _ in strata]
    with open(args.jobs) as fh:
        jobs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for job in jobs:
        labels, (tags_path, vcf, out) = job[:6], job[6:9]
        size, flip = blocks_and_orientation(vcf, truth_gt)
        tags = haplotag_list(tags_path, autosomes_only=True)
        if orient == "reads":
            # per block: labelled reads whose tag index equals the label, and the rest
            votes = {}
            for r in autosomal:
                t = tags.get(r)
                if t is None or t[0] not in ("H1", "H2") or size.get((t[2], t[1]), 0) < args.min_snps:
                    continue
                v = votes.setdefault((t[2], t[1]), [0, 0])
                v[0 if (0 if t[0] == "H1" else 1) == truth[r] else 1] += 1
        # reads: truth_tagged tagged compared correct unoriented, then the same in bases
        acc = {n: [0] * 10 for n in order}
        for r, regions in region_of.items():
            lab, ln = truth[r], length_of[r]
            t = tags.get(r)
            tagged = int(t is not None and t[0] in ("H1", "H2") and size.get((t[2], t[1]), 0) >= args.min_snps)
            comp = corr = unor = 0
            if tagged:
                if orient == "snps":
                    f = flip.get((t[2], t[1]))
                else:
                    # the block's orientation from the other labelled reads
                    same, opp = votes[(t[2], t[1])]
                    if (0 if t[0] == "H1" else 1) == lab:
                        same -= 1
                    else:
                        opp -= 1
                    f = None if same == opp else opp > same
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
                a[5] += ln
                a[6] += tagged * ln
                a[7] += comp * ln
                a[8] += corr * ln
        lines = []
        for n in order:
            tt, tg, cp, cr, un, tb, gb, cb, rb, _ = acc[n]
            if not tt:
                continue
            lines.append("\t".join(labels + [n, str(tt), str(tg), str(cp), str(cr), str(un),
                                             "%.6f" % (cr / tt), "%.6f" % (cr / cp if cp else 0),
                                             str(tb), str(gb), str(cb), str(rb),
                                             "%.6f" % (rb / tb if tb else 0),
                                             args.truth_source, orient]) + "\n")
        tmp = "%s.%d.tmp" % (out, os.getpid())
        with open(tmp, "w") as fh:
            fh.write("\t".join(COLUMNS) + "\n")
            fh.writelines(lines)
        os.replace(tmp, out)
        g = acc["genome"]
        print("  %s: genome recall %.4f (%d of %d reads), by bases %.4f"
              % (labels[5], g[3] / g[0] if g[0] else 0, g[3], g[0], g[8] / g[5] if g[5] else 0),
              file=sys.stderr)


if __name__ == "__main__":
    main()
