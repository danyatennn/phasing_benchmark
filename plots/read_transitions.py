#!/usr/bin/env python3
"""The state of each read on GRCh38 and on T2T-CHM13, the same reads, for one
platform and coverage: how many reads go from each state on GRCh38 to each
state on T2T-CHM13, per phaser. The numbers behind figR20 of
plots/figures_summary.py.

Reads: those with a parental label (bin/matpat_truth.sh, PAT or MAT; the label
does not depend on the reference) whose primary alignment is on chr1-22 on at
least one reference. A read missing from the read table of a reference counts
as unmapped there only if that cell's BAM has it unmapped (--unmapped-*);
otherwise the downsampling of that cell did not draw it and the read is left
out.

State on a reference, the first rule that fits (as plots/read_classes.py):
  unmapped   no alignment
  correct    tagged H1/H2 in a block of --min-snps phased het SNPs or more,
             on its parental haplotype once the block is oriented by its SNPs
             against the truth VCF (PAT/MAT put on the truth's haplotypes per
             chromosome through the truth haplotag list, as read_recall_strata.py)
  wrong      tagged, on the other haplotype
  lowmapq    MAPQ < 20
  no_snp     no truth het SNP on the read, and it touches the benchmark BED
  outside    no truth het SNP on the read, wholly outside the benchmark BED
  other      every other untagged read: het SNPs on it but not tagged, or its
             primary alignment off chr1-22
A tagged read in a block the truth cannot orient is left out (counted on stderr).

Output, --out <file>: phaser grch38 t2t reads

  python plots/read_transitions.py --run GRCh38 whatshap tags.tsv phased.vcf.gz \\
      --run T2T whatshap ... --labels reads/<cell>.matpat.tsv.gz \\
      --reads-grch38 ... --unmapped-grch38 ... --alleles-grch38 ... --truth-vcf-grch38 ... \\
      --truth-tags-grch38 ... --bench-grch38 ... (and the same for t2t) --out out.tsv
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import IntervalSet, haplotag_list, is_autosome, open_text, phased_het, vcf_records
from read_classes import blocks_and_orientation

MIN_MAPQ = 20
REFS = ("GRCh38", "T2T")
STATES = ["unmapped", "other", "no_snp", "outside", "lowmapq", "wrong", "correct"]


def parental_labels(paths):
    """read -> 0 (PAT) or 1 (MAT), from every matpat table given"""
    labels = {}
    for path in paths:
        with open_text(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[1] in ("PAT", "MAT"):
                    labels.setdefault(f[0], 0 if f[1] == "PAT" else 1)
    return labels


def alignments(path, wanted):
    """read -> (chrom, start, end, mapq) of its primary alignment, wanted reads only"""
    out = {}
    with open_text(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ic, i_s, ie, im = h.index("chrom"), h.index("start"), h.index("end"), h.index("mapq")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in wanted:
                out[f[0]] = (f[ic], int(f[i_s]), int(f[ie]), int(f[im]))
    return out


def names(path, wanted):
    """the wanted reads listed in a file of read names, one per line"""
    if not path or not os.path.exists(path):
        return set()
    with open_text(path) as fh:
        return {line.strip() for line in fh if line.strip() in wanted}


def snps_per_read(path, wanted):
    """read -> truth het SNPs it covers with REF or ALT"""
    out = {}
    with open_text(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        i_ref, i_alt = h.index("n_ref"), h.index("n_alt")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in wanted:
                out[f[0]] = int(f[i_ref]) + int(f[i_alt])
    return out


def parent_swap(truth_tags, labels):
    """chrom -> True if PAT is the truth's H2 there, as most reads labelled by both say"""
    both = {}
    for r, (hap, _, c) in haplotag_list(truth_tags, autosomes_only=True).items():
        if hap in ("H1", "H2") and r in labels:
            n = both.setdefault(c, [0, 0])
            n[0 if (0 if hap == "H1" else 1) == labels[r] else 1] += 1
    return {c: opp > same for c, (same, opp) in both.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", nargs=4, action="append", required=True,
                    metavar=("REF", "PHASER", "TAGS", "VCF"), help="GRCh38 or T2T, phaser, haplotag list, phased VCF")
    ap.add_argument("--labels", action="append", required=True, help="reads/<cell>.matpat.tsv.gz, repeatable")
    for ref in REFS:
        k = ref.lower()
        ap.add_argument("--reads-%s" % k, required=True, help="read table of the %s cell" % ref)
        ap.add_argument("--unmapped-%s" % k, default=None, help="names of the reads unmapped in the %s BAM" % ref)
        ap.add_argument("--alleles-%s" % k, required=True, help="<cell>_truth.reads.tsv.gz on %s" % ref)
        ap.add_argument("--truth-vcf-%s" % k, required=True)
        ap.add_argument("--truth-tags-%s" % k, required=True, help="<cell>_truthtags.tsv on %s" % ref)
        ap.add_argument("--bench-%s" % k, required=True, help="benchmark BED of %s" % ref)
    ap.add_argument("--min-snps", type=int, default=2)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    opt = lambda name, ref: getattr(args, "%s_%s" % (name, ref.lower()))

    phasers = []
    for _, phaser, _, _ in args.run:
        if phaser not in phasers:
            phasers.append(phaser)
    runs = {(ref, phaser): (tags, vcf) for ref, phaser, tags, vcf in args.run}
    for ref in REFS:
        for p in phasers:
            if (ref, p) not in runs:
                sys.exit("no %s run for %s" % (ref, p))

    labels = parental_labels(args.labels)
    pos = {ref: alignments(opt("reads", ref), labels) for ref in REFS}
    unmapped = {ref: names(opt("unmapped", ref), labels) for ref in REFS}
    reads = [r for r in labels
             if any(r in pos[ref] and is_autosome(pos[ref][r][0]) for ref in REFS)
             and all(r in pos[ref] or r in unmapped[ref] for ref in REFS)]
    print("%d labelled reads, %d on chr1-22 on a reference and drawn on both" % (len(labels), len(reads)),
          file=sys.stderr)
    wanted = set(reads)
    labels = {r: labels[r] for r in reads}
    snps = {ref: snps_per_read(opt("alleles", ref), wanted) for ref in REFS}

    # what a read is on a reference, whatever the phaser: None when it can be tagged
    base = {}
    for ref in REFS:
        bench = IntervalSet.from_bed(opt("bench", ref), autosomes_only=True)
        for r in reads:
            q = pos[ref].get(r)
            if q is None:
                s = "unmapped"
            elif not is_autosome(q[0]):
                s = "other"
            elif q[3] < MIN_MAPQ:
                s = "lowmapq"
            elif not snps[ref].get(r, 0):
                s = "no_snp" if bench.overlap(q[0], q[1], q[2]) > 0 else "outside"
            else:
                s = None
            base[(ref, r)] = s

    def states(ref, phaser, truth_gt, swap):
        """read -> state on ref for one run; None for a tagged read the truth cannot orient"""
        tags_path, vcf = runs[(ref, phaser)]
        size, flip = blocks_and_orientation(vcf, truth_gt)
        tags = haplotag_list(tags_path, autosomes_only=True)
        out = {}
        for r in reads:
            q = pos[ref].get(r)
            t = tags.get(r) if q is not None and is_autosome(q[0]) else None
            if t is not None and t[0] in ("H1", "H2") and size.get((t[2], t[1]), 0) >= args.min_snps:
                f, s_ = flip.get((t[2], t[1])), swap.get(t[2])
                if f is None or s_ is None:
                    out[r] = None
                else:
                    right = ((0 if t[0] == "H1" else 1) ^ f) == (labels[r] ^ s_)
                    out[r] = "correct" if right else "wrong"
            else:
                # untagged: the reason is the run's only for a read that could be tagged
                out[r] = base[(ref, r)] or "other"
        return out

    per_ref = {}
    for ref in REFS:
        truth_gt = {(c, p, a, b): gt for c, p, a, b, gt, ps in
                    vcf_records(opt("truth_vcf", ref), autosomes_only=True) if phased_het(gt) is not None}
        swap = parent_swap(opt("truth_tags", ref), labels)
        for p in phasers:
            per_ref[(ref, p)] = states(ref, p, truth_gt, swap)
        del truth_gt

    lines = []
    for p in phasers:
        g, t = per_ref[("GRCh38", p)], per_ref[("T2T", p)]
        counts, skipped = {}, 0
        for r in reads:
            a, b = g[r], t[r]
            if a is None or b is None:
                skipped += 1
                continue
            counts[(a, b)] = counts.get((a, b), 0) + 1
        for a in STATES:
            for b in STATES:
                lines.append("%s\t%s\t%s\t%d\n" % (p, a, b, counts.get((a, b), 0)))
        n = sum(counts.values())
        same = sum(v for (a, b), v in counts.items() if a == b)
        print("%s: %d reads, %.2f%% in the same state on both, %d left out (block not oriented)"
              % (p, n, 100 * same / n if n else 0, skipped), file=sys.stderr)
    tmp = "%s.%d.tmp" % (args.out, os.getpid())
    with open(tmp, "w") as fh:
        fh.write("phaser\tgrch38\tt2t\treads\n")
        fh.writelines(lines)
    os.replace(tmp, args.out)
    print("written %s" % args.out, file=sys.stderr)


if __name__ == "__main__":
    main()
