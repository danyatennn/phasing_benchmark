#!/usr/bin/env python3
"""Where along the chromosomes a run's phase switches, for every run of one
cell - the numbers behind plots/figures_switch.py. Four files per run, in
--out-dir, all over chr1-22:

  <run>.switch_windows.tsv   chrom start end assessed switches
      assessed: pairs of consecutive phased het SNPs of one block that are both
      truth het SNPs (same alleles); switches: those pairs phased the other way
      round relative to each other than in the truth. A pair goes to the window
      of its midpoint. This is whatshap compare's switch count, per window:
      their sums match all_switches / all_assessed_pairs of <run>.compare.tsv.
  <run>.switch_positions.tsv chrom pos1 pos2 - the two SNPs of every switch
  <run>.blocks.tsv           chrom ps start end snps - every phase block of
                             --min-snps phased het SNPs or more
  <run>.read_depth.tsv.gz    chrom start end concordant discordant - mean
      depth of the tagged reads (H1/H2 in a block of --min-snps or more, with a
      truth label) whose label is the truth's (concordant) or not (discordant),
      once each block is oriented against the truth by the majority of its
      SNPs. So a block boundary does not flip the colour; a switch inside a
      block does: past it the reads are on the other haplotype. Needs the
      cell's read table (positions); left out without one.

  python plots/switch_tracks.py --truth-vcf truth.phased.vcf.gz \\
      --truth-tags <cell>_truthtags.tsv --reads-table <cell>.reads.tsv.gz \\
      --jobs runs.tsv --out-dir metrics_switchtracks
runs.tsv: run name, haplotag list, phased VCF
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from benchio import (AUTOSOMES, haplotag_list, is_autosome, open_text, phased_het,
                     truth_labels, vcf_records, write_text)


def load_truth(path):
    """(chrom, pos, ref, alt) -> GT of the phased truth het SNPs"""
    return {(c, p, ref, alt): gt for c, p, ref, alt, gt, ps in vcf_records(path, autosomes_only=True)
            if phased_het(gt) is not None}


def switches(vcf, truth_gt, min_snps):
    """Blocks, the pairs assessed and the switches of one phased VCF.
    Returns blocks {(chrom, ps): [start, end, snps]}, pairs [(chrom, mid, switch)],
    positions [(chrom, pos1, pos2)] and the block orientation {(chrom, ps): flip}."""
    blocks, pairs, positions = {}, [], []
    last = {}                      # (chrom, ps) -> (pos, same) of the last assessed SNP
    same_n, opp_n = {}, {}
    for chrom, pos, ref, alt, gt, ps in vcf_records(vcf, autosomes_only=True):
        if ps in (".", "") or phased_het(gt) is None:
            continue
        b = (chrom, ps)
        x = blocks.setdefault(b, [pos, pos, 0])
        x[1] = pos
        x[2] += 1
        t = truth_gt.get((chrom, pos, ref, alt))
        if t is None:
            continue
        a, c = gt.split("|")
        if gt == t:
            same = True
        elif "%s|%s" % (c, a) == t:
            same = False
        else:
            continue
        (same_n if same else opp_n)[b] = (same_n if same else opp_n).get(b, 0) + 1
        prev = last.get(b)
        if prev is not None:
            switched = prev[1] != same
            pairs.append((chrom, (prev[0] + pos) // 2, switched))
            if switched:
                positions.append((chrom, prev[0], pos))
        last[b] = (pos, same)
    # a pair needs two SNPs of a block, so every pair is in a block of two or
    # more; blocks under --min-snps are dropped from the block list and the
    # orientation only
    small = {b for b, x in blocks.items() if x[2] < min_snps}
    blocks = {b: x for b, x in blocks.items() if b not in small}
    flip = {b: opp_n.get(b, 0) > same_n.get(b, 0) for b in set(same_n) | set(opp_n) if b not in small}
    return blocks, pairs, positions, flip


def windows(pairs, lengths, size):
    """chrom -> (assessed, switches) arrays over windows of `size`"""
    out = {c: (np.zeros(lengths[c] // size + 1, dtype=np.int64),
               np.zeros(lengths[c] // size + 1, dtype=np.int64)) for c in lengths}
    for chrom, mid, switched in pairs:
        if chrom not in out:
            continue
        w = mid // size
        out[chrom][0][w] += 1
        out[chrom][1][w] += switched
    return out


def read_depth(reads_table, tags, truth, size_of, flip, lengths, bin_size, min_snps):
    """chrom -> (concordant, discordant) depth arrays over bins of bin_size,
    from the read positions of the read table"""
    conc = {c: np.zeros(lengths[c] // bin_size + 1) for c in lengths}
    disc = {c: np.zeros(lengths[c] // bin_size + 1) for c in lengths}
    with open_text(reads_table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_chrom, i_start, i_end = header.index("chrom"), header.index("start"), header.index("end")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom = f[i_chrom]
            if chrom not in conc:
                continue
            t = tags.get(f[0])
            if t is None or t[0] not in ("H1", "H2") or f[0] not in truth:
                continue
            b = (chrom, t[1])
            if size_of.get(b, 0) < min_snps or b not in flip:
                continue
            label = (0 if t[0] == "H1" else 1) ^ flip[b]
            arr = conc[chrom] if label == truth[f[0]] else disc[chrom]
            s, e = int(f[i_start]), int(f[i_end])
            if e <= s:
                continue
            b0, b1 = s // bin_size, (e - 1) // bin_size
            if b0 == b1:
                arr[b0] += e - s
            else:
                arr[b0] += (b0 + 1) * bin_size - s
                arr[b0 + 1:b1] += bin_size
                arr[b1] += e - b1 * bin_size
    for c in conc:
        conc[c] /= bin_size
        disc[c] /= bin_size
    return conc, disc


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--truth-vcf", required=True)
    parser.add_argument("--truth-tags", required=True)
    parser.add_argument("--reads-table", default=None, help="for the read depth; left out without it")
    parser.add_argument("--chrom-lengths", required=True, help="chrom <tab> length")
    parser.add_argument("--jobs", required=True, help="run name, haplotag list, phased VCF")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--window", type=int, default=1_000_000)
    parser.add_argument("--depth-bin", type=int, default=10_000)
    parser.add_argument("--min-snps", type=int, default=2)
    args = parser.parse_args()

    lengths = {}
    with open(args.chrom_lengths) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2 and is_autosome(f[0]):
                lengths[f[0]] = int(f[1])
    order = [c for c in AUTOSOMES if c in lengths]
    truth_gt = load_truth(args.truth_vcf)
    truth = truth_labels(args.truth_tags, "truthtags", autosomes_only=True)
    print("%d truth het SNPs, %d reads with a truth label" % (len(truth_gt), len(truth)), file=sys.stderr)
    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.jobs) as fh:
        jobs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for name, tags_path, vcf in jobs:
        base = os.path.join(args.out_dir, name)
        blocks, pairs, positions, flip = switches(vcf, truth_gt, args.min_snps)

        win = windows(pairs, lengths, args.window)
        out = write_text(base + ".switch_windows.tsv")
        out.write("chrom\tstart\tend\tassessed\tswitches\n")
        for c in order:
            assessed, sw = win[c]
            for w in range(len(assessed)):
                out.write("%s\t%d\t%d\t%d\t%d\n" % (c, w * args.window, min(lengths[c], (w + 1) * args.window),
                                                     assessed[w], sw[w]))
        out.commit()

        out = write_text(base + ".switch_positions.tsv")
        out.write("chrom\tpos1\tpos2\n")
        out.write("".join("%s\t%d\t%d\n" % p for p in positions))
        out.commit()

        out = write_text(base + ".blocks.tsv")
        out.write("chrom\tps\tstart\tend\tsnps\n")
        for (c, ps), (s, e, n) in sorted(blocks.items(), key=lambda kv: (kv[0][0], kv[1][0])):
            out.write("%s\t%s\t%d\t%d\t%d\n" % (c, ps, s, e, n))
        out.commit()

        msg = "%s: %d blocks, %d switches in %d pairs" % (
            name, len(blocks), len(positions), len(pairs))

        if args.reads_table and os.path.exists(args.reads_table):
            tags = haplotag_list(tags_path, autosomes_only=True)
            size_of = {b: x[2] for b, x in blocks.items()}
            conc, disc = read_depth(args.reads_table, tags, truth, size_of, flip, lengths,
                                    args.depth_bin, args.min_snps)
            out = write_text(base + ".read_depth.tsv.gz")
            out.write("chrom\tstart\tend\tconcordant\tdiscordant\n")
            for c in order:
                for i in range(len(conc[c])):
                    out.write("%s\t%d\t%d\t%.2f\t%.2f\n" % (c, i * args.depth_bin,
                              min(lengths[c], (i + 1) * args.depth_bin), conc[c][i], disc[c][i]))
            out.commit()
            msg += ", read depth written"
        print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
