#!/usr/bin/env python3
"""How many het SNPs a read hits, and how often its haplotag is right, for
every run of one cell and one phaser input - the numbers behind
figR3 of plots/figures_reads.py.

Reads: those of the run's haplotag list on chr1-22 (whatshap lists every read
it saw, tagged or not). For each read, two counts:
  snps    het SNPs of the input it covers with REF or ALT (allele matrix,
          reads/<cell>_<input>.reads.tsv.gz: n_ref + n_alt); a read missing
          from the matrix covers none. The same for every phaser.
  votes   phased SNPs of the phaser's block that carry an allele of it (votes
          matrix, h1 + h2): of the block the read was tagged in, or for an
          untagged read the counting block with most votes. What haplotag
          actually decided on.
A read is tagged if it has H1/H2 in a block of --min-snps phased het SNPs or
more. It is compared if it also has a truth label (whatshap haplotag against
the phased truth) and its block shares SNPs with the truth, which orient it
(majority of the shared SNPs, as readacc_run.sh); correct if, so oriented, its
label is the truth's.

Output, one line per run x count x value (counts above --cap go to --cap):
  reference sample caller platform depth phaser count n reads tagged compared correct

  python plots/snps_accuracy.py --alleles reads/<cell>_truth --truth-tags <cell>_truthtags.tsv \\
      --truth-vcf truth.phased.vcf.gz --min-snps 2 --jobs runs.tsv --out table.tsv
runs.tsv: reference sample caller platform depth phaser tags vcf
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import haplotag_list, open_text, phased_het, truth_labels, vcf_records
from read_classes import blocks_and_orientation, votes


def snps_per_read(prefix):
    """read -> het SNPs of the input it covers with REF or ALT"""
    out = {}
    path = prefix + ".reads.tsv.gz"
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_ref, i_alt = header.index("n_ref"), header.index("n_alt")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = int(f[i_ref]) + int(f[i_alt])
    return out


def count_run(tags_path, vcf, vts, snps, truth, truth_gt, min_snps, cap):
    size, flip = blocks_and_orientation(vcf, truth_gt)
    tags = haplotag_list(tags_path, autosomes_only=True)
    # (count, n) -> [reads, tagged, compared, correct]
    acc = {}
    for r, (label, ps, chrom) in tags.items():
        tagged = label in ("H1", "H2") and size.get((chrom, ps), 0) >= min_snps
        rv = vts.get(r, {}) if vts is not None else {}
        if tagged:
            h = rv.get(ps, (0, 0))
        else:
            counting = [v for p, v in rv.items() if size.get((chrom, p), 0) >= min_snps]
            h = max(counting, key=lambda v: v[0] + v[1]) if counting else (0, 0)
        comp = corr = 0
        if tagged and r in truth:
            f = flip.get((chrom, ps))
            if f is not None:
                comp = 1
                corr = 1 if ((0 if label == "H1" else 1) != f) == truth[r] else 0
        for count, n in (("snps", snps.get(r, 0)), ("votes", h[0] + h[1])):
            x = acc.setdefault((count, min(n, cap)), [0, 0, 0, 0])
            x[0] += 1
            x[1] += tagged
            x[2] += comp
            x[3] += corr
    return acc


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--alleles", required=True, help="allele matrix prefix reads/<cell>_<input>")
    parser.add_argument("--truth-tags", required=True)
    parser.add_argument("--truth-vcf", required=True)
    parser.add_argument("--min-snps", type=int, default=2)
    parser.add_argument("--cap", type=int, default=1000)
    parser.add_argument("--jobs", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    snps = snps_per_read(args.alleles)
    truth = truth_labels(args.truth_tags, "truthtags", autosomes_only=True)
    truth_gt = {(c, p, ref, alt): gt for c, p, ref, alt, gt, ps in
                vcf_records(args.truth_vcf, autosomes_only=True) if phased_het(gt) is not None}
    print("%s: %d reads with SNPs, %d with a truth label" % (args.alleles, len(snps), len(truth)),
          file=sys.stderr)

    lines = []
    with open(args.jobs) as fh:
        jobs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for job in jobs:
        labels, (tags, vcf) = job[:6], job[6:8]
        vts = votes(args.alleles, labels[5])
        if vts is None:
            print("WARNING: no votes for %s in %s, 'votes' will read 0" % (labels[5], args.alleles),
                  file=sys.stderr)
        acc = count_run(tags, vcf, vts, snps, truth, truth_gt, args.min_snps, args.cap)
        for (count, n), x in sorted(acc.items()):
            lines.append("\t".join(labels + [count, str(n)] + [str(v) for v in x]) + "\n")
        g = [sum(x[i] for (c, _), x in acc.items() if c == "snps") for i in range(4)]
        print("  %s: %d reads, %d tagged, %d of %d correct" % (labels[5], g[0], g[1], g[3], g[2]),
              file=sys.stderr)

    tmp = "%s.%d.tmp" % (args.out, os.getpid())
    with open(tmp, "w") as fh:
        fh.write("reference\tsample\tcaller\tplatform\tdepth\tphaser\tcount\tn"
                 "\treads\ttagged\tcompared\tcorrect\n")
        fh.writelines(lines)
    os.replace(tmp, args.out)


if __name__ == "__main__":
    main()
