#!/usr/bin/env python3
"""Why reads are not assigned a haplotype, how sure an assignment is, and what
the reads in each group look like - for every run of one cell and one phaser
input, genome-wide and in each stratum.

Reads: primary alignments on chr1-22 (the read table of the cell).

Class - exactly one per read, the first rule that fits:
  tagged              H1/H2 from whatshap haplotag, in a block of at least
                      --min-snps phased het SNPs
  lowmapq_paralog     MAPQ < 20, half the alignment or more in segmental
                      duplications or low-mappability regions
  lowmapq_other       MAPQ < 20 elsewhere
  no_snp_outside_bed  no het SNP of the input on the read, and the read lies
                      wholly outside the confident regions
  no_snp              no het SNP of the input on the read
  snps_unphased       het SNPs on it, but none in a block that counts: left
                      unphased, or phased only as single-SNP blocks
  phased_unassigned   SNPs of a block that counts on it, yet no haplotype

Confidence - from the read's votes (allele matrix): at each phased het SNP the
read carries the allele of haplotype 1 or of haplotype 2 of the block. The
block is the one the read was tagged in, or, untagged, the counting block with
most votes.
  high    two or more votes, all for one haplotype
  medium  a majority, with votes against
  single  one vote
  tie     as many votes for each
  none    no vote

Accuracy: a tagged read with a truth label is compared after its block is
oriented by its SNPs against the truth (as readacc_run.sh does).

Features, per group: reads, bases, mean length and MAPQ, soft-clipped share of
the bases, edit distance (NM), mismatches and indels per aligned base, mean
base quality, mean |VAF - 0.5| of the SNPs on the reads, share of the
neighbouring-SNP links on them showing more than two allele combinations.

A read belongs to a stratum if half its alignment or more lies in it.

  python plots/read_classes.py --reads-table <cell>.reads.tsv.gz \\
      --alleles reads/<cell>_<caller> --truth-tags <cell>_truthtags.tsv \\
      --truth-vcf truth.phased.vcf.gz --confident-bed bed \\
      --stratum segdups=prepared/GRCh38_segdups.bed ... \\
      --paralog segdups --paralog lowmappabilityall \\
      --min-snps 2 --jobs runs.tsv
runs.tsv: reference sample caller platform depth phaser tags vcf out
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import (IntervalSet, haplotag_list, open_text, phased_het, read_table,
                     truth_labels, vcf_records)

CLASSES = ["tagged", "lowmapq_paralog", "lowmapq_other", "no_snp_outside_bed", "no_snp",
           "snps_unphased", "phased_unassigned"]
CONFIDENCE = ["high", "medium", "single", "tie", "none"]
MIN_MAPQ = 20
STRATUM_FRACTION = 0.5

OUT_COLUMNS = ["reference", "sample", "caller", "platform", "depth", "phaser", "stratum",
               "group_type", "group", "reads", "bases", "mean_length", "mean_mapq",
               "softclip_frac", "nm_rate", "mismatch_rate", "indel_rate", "mean_bq",
               "vaf_dev", "paralog_link_frac", "compared", "correct"]

# accumulator slots
(N, BASES, LEN, MAPQ, SOFT, ALN, NM, NM_ALN, MM, INDEL, BQ, NBQ, VAF, NVAF, LINKS, PLINKS,
 COMP, CORR) = range(18)


def new_acc():
    return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0, 0.0, 0, 0, 0, 0, 0]


def allele_reads(prefix):
    """read -> (n_sites, vaf_dev or None, links, paralog_links)"""
    out = {}
    path = prefix + ".reads.tsv.gz"
    if not os.path.exists(path):
        return out
    with open_text(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = (int(f[2]), None if f[7] == "nan" else float(f[7]), int(f[9]), int(f[10]))
    return out


def votes(prefix, phaser):
    """read -> {PS: (h1, h2)}"""
    out = {}
    path = "%s.votes.%s.tsv.gz" % (prefix, phaser)
    if not os.path.exists(path):
        return None
    with open_text(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.setdefault(f[0], {})[f[2]] = (int(f[3]), int(f[4]))
    return out


def blocks_and_orientation(vcf, truth_gt):
    """(chrom, PS) -> number of phased het SNPs, and -> flip (True/False) for
    blocks that share a SNP with the truth, majority of their SNPs deciding."""
    size, same, opp = {}, {}, {}
    for chrom, pos, ref, alt, gt, ps in vcf_records(vcf, autosomes_only=True):
        if ps in (".", "") or phased_het(gt) is None:
            continue
        b = (chrom, ps)
        size[b] = size.get(b, 0) + 1
        t = truth_gt.get((chrom, pos, ref, alt))
        if t is None:
            continue
        a, c = gt.split("|")
        if gt == t:
            same[b] = same.get(b, 0) + 1
        elif "%s|%s" % (c, a) == t:
            opp[b] = opp.get(b, 0) + 1
    flip = {b: opp.get(b, 0) > same.get(b, 0) for b in set(same) | set(opp)}
    return size, flip


def confidence(h1, h2):
    n = h1 + h2
    if n == 0:
        return "none"
    if n == 1:
        return "single"
    if h1 == h2:
        return "tie"
    if min(h1, h2) == 0:
        return "high"
    return "medium"


def classify_run(reads, read_strata, in_bed, paralog, alle, truth, tags_path, vcf, vts,
                 truth_gt, min_snps):
    size, flip = blocks_and_orientation(vcf, truth_gt)
    tags = haplotag_list(tags_path, autosomes_only=True)
    acc = {}

    def add(stratum, gtype, group, r, info, a, comp, corr):
        k = (stratum, gtype, group)
        x = acc.get(k)
        if x is None:
            x = acc[k] = new_acc()
        chrom, start, end, length, mapq, aligned, soft, nm, ins, dels, bq = info
        x[N] += 1
        x[BASES] += length
        x[LEN] += length
        x[MAPQ] += mapq
        x[SOFT] += soft
        x[ALN] += aligned
        if nm >= 0 and aligned > 0:
            x[NM] += nm
            x[NM_ALN] += aligned
            x[MM] += max(0, nm - ins - dels)
            x[INDEL] += ins + dels
        if bq == bq:
            x[BQ] += bq
            x[NBQ] += 1
        if a is not None:
            if a[1] is not None:
                x[VAF] += a[1]
                x[NVAF] += 1
            x[LINKS] += a[2]
            x[PLINKS] += a[3]
        x[COMP] += comp
        x[CORR] += corr

    for r, info in reads.items():
        chrom, start, end, length, mapq = info[0], info[1], info[2], info[3], info[4]
        a = alle.get(r)
        rv = vts.get(r, {}) if vts is not None else {}
        label, ps, _ = tags.get(r, ("none", "none", chrom))
        tagged = label in ("H1", "H2") and size.get((chrom, ps), 0) >= min_snps

        counting = {p: v for p, v in rv.items() if size.get((chrom, p), 0) >= min_snps}
        if tagged:
            cls = "tagged"
            h = rv.get(ps, (0, 0))
        else:
            if mapq < MIN_MAPQ:
                cls = "lowmapq_paralog" if r in paralog else "lowmapq_other"
            elif a is None or a[0] == 0:
                cls = "no_snp" if r in in_bed else "no_snp_outside_bed"
            elif not any(h1 + h2 for h1, h2 in counting.values()):
                cls = "snps_unphased"
            else:
                cls = "phased_unassigned"
            h = max(counting.values(), key=lambda v: v[0] + v[1]) if counting else (0, 0)
        conf = confidence(h[0], h[1])

        comp = corr = 0
        if tagged and r in truth:
            f = flip.get((chrom, ps))
            if f is not None:
                comp = 1
                lab = 0 if label == "H1" else 1
                corr = 1 if (lab != f) == truth[r] else 0

        for stratum in ("genome",) + read_strata.get(r, ()):
            add(stratum, "all", "all", r, info, a, comp, corr)
            add(stratum, "class", cls, r, info, a, comp, corr)
            add(stratum, "confidence", conf, r, info, a, comp, corr)
            if tagged:
                add(stratum, "tagged_confidence", conf, r, info, a, comp, corr)
    return acc


def fmt(x, digits=6):
    return ("%%.%df" % digits) % x


def rows(labels, acc, strata_order):
    order = {"all": ["all"], "class": CLASSES, "confidence": CONFIDENCE,
             "tagged_confidence": CONFIDENCE}
    out = []
    for stratum in ["genome"] + strata_order:
        for gtype in ("all", "class", "confidence", "tagged_confidence"):
            for group in order[gtype]:
                x = acc.get((stratum, gtype, group))
                if x is None:
                    continue
                n = x[N]
                out.append("\t".join(labels + [stratum, gtype, group, str(n), str(x[BASES]),
                    fmt(x[LEN] / n, 1), fmt(x[MAPQ] / n, 2),
                    fmt(x[SOFT] / x[LEN]) if x[LEN] else "nan",
                    fmt(x[NM] / x[NM_ALN]) if x[NM_ALN] else "nan",
                    fmt(x[MM] / x[NM_ALN]) if x[NM_ALN] else "nan",
                    fmt(x[INDEL] / x[NM_ALN]) if x[NM_ALN] else "nan",
                    fmt(x[BQ] / x[NBQ], 2) if x[NBQ] else "nan",
                    fmt(x[VAF] / x[NVAF], 4) if x[NVAF] else "nan",
                    fmt(x[PLINKS] / x[LINKS]) if x[LINKS] else "nan",
                    str(x[COMP]), str(x[CORR])]) + "\n")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reads-table", required=True)
    parser.add_argument("--alleles", required=True, help="allele matrix prefix reads/<cell>_<caller>")
    parser.add_argument("--truth-tags", required=True)
    parser.add_argument("--truth-vcf", required=True)
    parser.add_argument("--confident-bed", required=True)
    parser.add_argument("--stratum", action="append", default=[], help="name=bed, repeatable")
    parser.add_argument("--paralog", action="append", default=[], help="stratum names counted as paralogous")
    parser.add_argument("--min-snps", type=int, default=2)
    parser.add_argument("--jobs", required=True, help="TSV: reference sample caller platform depth phaser tags vcf out")
    args = parser.parse_args()

    reads = read_table(args.reads_table, autosomes_only=True)
    print("%s: %d reads on chr1-22" % (args.reads_table, len(reads)), file=sys.stderr)

    strata = []
    for item in args.stratum:
        name, _, path = item.partition("=")
        if os.path.exists(path):
            strata.append((name, IntervalSet.from_bed(path, autosomes_only=True)))
    bed = IntervalSet.from_bed(args.confident_bed, autosomes_only=True)
    par_names = set(args.paralog)
    par_rows = []
    for name, s in strata:
        if name in par_names:
            for c in s.starts:
                par_rows += [(c, a, b) for a, b in zip(s.starts[c], s.ends[c])]
    paralog_set = IntervalSet(par_rows)

    # once per cell: which strata each read belongs to, which reads are
    # paralogous, which touch the confident regions
    read_strata, paralog, in_bed = {}, set(), set()
    for r, info in reads.items():
        chrom, start, end = info[0], info[1], info[2]
        span = max(1, end - start)
        names = tuple(n for n, s in strata if s.overlap(chrom, start, end) >= STRATUM_FRACTION * span)
        if names:
            read_strata[r] = names
        if paralog_set.overlap(chrom, start, end) >= STRATUM_FRACTION * span:
            paralog.add(r)
        if bed.overlap(chrom, start, end) > 0:
            in_bed.add(r)

    truth = truth_labels(args.truth_tags, "truthtags", autosomes_only=True)
    truth_gt = {(c, p, ref, alt): gt for c, p, ref, alt, gt, ps in
                vcf_records(args.truth_vcf, autosomes_only=True) if phased_het(gt) is not None}
    alle = allele_reads(args.alleles)
    strata_order = [n for n, _ in strata]

    with open(args.jobs) as fh:
        jobs = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for job in jobs:
        labels, (tags, vcf, out) = job[:6], job[6:9]
        vts = votes(args.alleles, labels[5])
        if vts is None:
            print("WARNING: no votes for %s in %s, confidence will read 'none'" % (labels[5], args.alleles),
                  file=sys.stderr)
        acc = classify_run(reads, read_strata, in_bed, paralog, alle, truth, tags, vcf, vts,
                           truth_gt, args.min_snps)
        tmp = "%s.%d.tmp" % (out, os.getpid())
        with open(tmp, "w") as fh:
            fh.write("\t".join(OUT_COLUMNS) + "\n")
            fh.writelines(rows(labels, acc, strata_order))
        os.replace(tmp, out)
        g = acc.get(("genome", "class", "tagged"))
        allr = acc.get(("genome", "all", "all"))
        print("%s: %d reads, %d tagged" % (os.path.basename(out), allr[N] if allr else 0, g[N] if g else 0),
              file=sys.stderr)


if __name__ == "__main__":
    main()
