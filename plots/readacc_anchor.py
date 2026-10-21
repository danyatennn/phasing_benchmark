#!/usr/bin/env python3
"""Read-level phasing accuracy with every phase block oriented by an anchor
read, for all runs of one cell (one BAM, so one read table and one truth).

Haplotype labels are arbitrary within a block, so a block is oriented before
its reads are compared with the truth. readacc_run.sh orients it by its SNPs
against the truth genotypes; here it is oriented by one read:

  1. a block is (chromosome, PS) of the run's VCF; blocks of fewer than
     --min-snps phased het SNPs orient nothing, and their reads count as
     untagged (Margin writes a phase set even for a single SNP)
  2. candidates: reads tagged H1/H2 in the block that also have a truth label
  3. the anchor is the longest candidate; ties go to the one whose alignment
     covers more of the block's SNPs, then to the smaller read name
  4. the block is flipped if the anchor's label disagrees with its truth label
  5. every other candidate is scored with that orientation; the anchor is not
     scored (it would always be right), so a block with one candidate scores
     nothing

Truth labels: whatshap haplotag against the phased truth (truthtags), or the
parental-haplotype tables of bin/matpat_truth.sh (matpat) and bin/yak_truth.sh
(yak). With --autosomes, reads and blocks off chr1-22 are left out.

One output line per run:
  reference region caller platform depth phaser truth_source
  reads_truth_labelled   reads with a truth label (the recall denominator)
  reads_tagged           reads tagged by the run in a block that counts
  reads_scored           candidates other than anchors
  reads_correct          of those, on the truth haplotype after orientation
  anchors                blocks with an anchor, i.e. reads not scored as anchors
  blocks_flipped         of those blocks, the ones flipped
  read_precision         reads_correct / reads_scored
  read_recall            reads_correct / reads_truth_labelled

  python plots/readacc_anchor.py --reads-table <cell>.reads.tsv.gz \\
      --truth <cell>_truthtags.tsv --truth-source truthtags \\
      --min-snps 2 --autosomes 1 --jobs runs.tsv
runs.tsv: reference region caller platform depth phaser tags vcf out
"""

import argparse
import bisect
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import haplotag_list, is_autosome, phase_blocks, read_table, truth_labels

ANCHOR_COLUMNS = ["reference", "region", "caller", "platform", "depth", "phaser", "truth_source",
                  "reads_truth_labelled", "reads_tagged", "reads_scored", "reads_correct",
                  "anchors", "blocks_flipped", "read_precision", "read_recall"]


def score_run(reads, truth, tags_path, vcf_path, min_snps, autosomes):
    blocks = phase_blocks(vcf_path, autosomes_only=autosomes)
    tags = haplotag_list(tags_path, autosomes_only=autosomes)

    members = {}          # block -> [(read, label 0/1)]
    tagged = 0
    for read, (label, ps, chrom) in tags.items():
        if label not in ("H1", "H2"):
            continue
        positions = blocks.get((chrom, ps))
        if positions is None or len(positions) < min_snps:
            continue
        tagged += 1
        if read in truth:
            members.setdefault((chrom, ps), []).append((read, 0 if label == "H1" else 1))

    scored = correct = anchors = flipped = 0
    for block, cands in members.items():
        positions = blocks[block]

        def rank(item):
            info = reads.get(item[0])
            if info is None:
                return (0, 0, item[0])
            start, end, length = info[1], info[2], info[3]
            covered = bisect.bisect_left(positions, end + 1) - bisect.bisect_left(positions, start + 1)
            # longest first, then most SNPs, then the smaller name
            return (-length, -covered, item[0])

        cands.sort(key=rank)
        anchor_read, anchor_label = cands[0]
        flip = anchor_label != truth[anchor_read]
        anchors += 1
        flipped += flip
        for read, label in cands[1:]:
            scored += 1
            if (label != flip) == truth[read]:     # label xor flip == truth
                correct += 1
    return len(truth), tagged, scored, correct, anchors, flipped


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reads-table", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--truth-source", required=True, choices=["truthtags", "matpat", "yak"])
    parser.add_argument("--min-snps", type=int, default=2)
    parser.add_argument("--autosomes", type=int, default=1, choices=[0, 1])
    parser.add_argument("--jobs", required=True, help="TSV: reference region caller platform depth phaser tags vcf out")
    args = parser.parse_args()

    reads = read_table(args.reads_table, autosomes_only=bool(args.autosomes))
    read_chrom = {r: v[0] for r, v in reads.items()}
    truth = truth_labels(args.truth, args.truth_source, bool(args.autosomes), read_chrom)
    print("%s: %d reads, %d with a %s label" % (args.reads_table, len(reads), len(truth),
                                               args.truth_source), file=sys.stderr)

    with open(args.jobs) as fh:
        jobs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for job in jobs:
        labels, (tags, vcf, out) = job[:6], job[6:9]
        n_truth, tagged, scored, correct, anchors, flipped = score_run(
            reads, truth, tags, vcf, args.min_snps, bool(args.autosomes))
        fields = labels + [args.truth_source, n_truth, tagged, scored, correct, anchors, flipped,
                           "%.6f" % (correct / scored if scored else 0),
                           "%.6f" % (correct / n_truth if n_truth else 0)]
        tmp = "%s.%d.tmp" % (out, os.getpid())
        with open(tmp, "w") as fh:
            fh.write("\t".join(str(f) for f in fields) + "\n")
        os.replace(tmp, out)
        print("%s: precision %s, recall %s" % (os.path.basename(out), fields[-2], fields[-1]),
              file=sys.stderr)


if __name__ == "__main__":
    main()
