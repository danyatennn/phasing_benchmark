#!/usr/bin/env python3
"""VAF of the het SNPs of an input as the reads of each cell see them, binned,
genome-wide and in each stratum - the numbers behind the VAF panels of
plots/figures_regions.py, small enough to copy anywhere.

The VAF of a site is ALT / (REF + ALT) from the allele matrix of the cell
(bin/allele_matrix.sh: samtools mpileup, every primary read, MAPQ 0 included,
so the reads of a paralogous copy show up as a VAF away from 0.5). Sites with
fewer than --min-reads REF + ALT reads are left out: at a handful of reads the
VAF is mostly noise. Only chr1-22.

"genome" is every site of the input on chr1-22 (for the truth: inside the
confident regions). A site is in a stratum if its position is in the stratum's
BED (regions/prepared_auto, cut to the confident regions).

A VAF of k ALT reads out of n can only be one of n + 1 values, and at odd n
never exactly 0.5; binned as points, that lattice shows up as spikes and dips
around 0.5 that are not in the data. So every site is spread evenly over
[(k - 1/2) / n, (k + 1/2) / n], cut to [0, 1] - a continuity correction: the
shape is the same, the lattice is gone, and each site still weighs 1.

Output, long format, one line per cell x input x stratum x bin:
  reference sample platform depth input stratum vaf_lo vaf_hi sites
sites is a weight, so not a whole number; it sums to the number of sites.

  python plots/vaf_strata.py --reference GRCh38 --sample HG002 \\
      --stratum segdups=prepared_auto/GRCh38_segdups.bed ... \\
      --jobs jobs.tsv --out metrics_strata_auto/vaf_GRCh38.tsv
jobs.tsv lines:
  platform depth input sites.tsv.gz
"""

import argparse
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import IntervalSet, is_autosome, open_text

# set in every worker by init(), so it works under fork and spawn alike
STRATA = []        # [(name, IntervalSet)]
BINS = 100
MIN_READS = 5


def load_strata(items):
    out = []
    for name, path in items:
        out.append((name, IntervalSet.from_bed(path, autosomes_only=True)))
    return out


def init(items, bins, min_reads):
    global STRATA, BINS, MIN_READS
    STRATA, BINS, MIN_READS = load_strata(items), bins, min_reads


def spread(k, n):
    """[(bin, weight)]: one site, k ALT reads of n, spread evenly over
    [(k - 1/2) / n, (k + 1/2) / n] cut to [0, 1]."""
    lo, hi = max(0.0, (k - 0.5) / n), min(1.0, (k + 0.5) / n)
    out = []
    b = min(BINS - 1, int(lo * BINS))
    while b < BINS and b / BINS < hi:
        overlap = min(hi, (b + 1) / BINS) - max(lo, b / BINS)
        if overlap > 0:
            out.append((b, overlap / (hi - lo)))
        b += 1
    return out


def add(counts, share):
    for b, w in share:
        counts[b] += w


def histogram(job):
    platform, depth, source, sites = job
    counts = {"genome": [0.0] * BINS}
    for name, _ in STRATA:
        counts[name] = [0.0] * BINS
    with open_text(sites) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_ref, i_alt = header.index("n_ref"), header.index("n_alt")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not is_autosome(f[0]):
                continue
            n_ref, n_alt = int(f[i_ref]), int(f[i_alt])
            if n_ref + n_alt < MIN_READS:
                continue
            share = spread(n_alt, n_ref + n_alt)
            add(counts["genome"], share)
            pos0 = int(f[1]) - 1
            for name, s in STRATA:
                if s.contains(f[0], pos0):
                    add(counts[name], share)
    print("%s %sx %s: %d sites" % (platform, depth, source, round(sum(counts["genome"]))), file=sys.stderr)
    return platform, depth, source, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reference", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--stratum", action="append", default=[], help="name=bed")
    parser.add_argument("--jobs", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--bins", type=int, default=BINS, help="bins over [0, 1]")
    parser.add_argument("--min-reads", type=int, default=MIN_READS, help="REF + ALT reads at a site")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    items = []
    for item in args.stratum:
        name, _, path = item.partition("=")
        if os.path.exists(path):
            items.append((name, path))
        else:
            print("WARNING: no BED for %s: %s" % (name, path), file=sys.stderr)

    jobs = []
    with open(args.jobs) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and os.path.exists(f[3]):
                jobs.append(tuple(f[:4]))
    if not jobs:
        sys.exit("no allele matrix found for any job in %s" % args.jobs)

    with Pool(min(args.threads, len(jobs)), initializer=init,
              initargs=(items, args.bins, args.min_reads)) as pool:
        results = pool.map(histogram, jobs)

    order = ["genome"] + [n for n, _ in items]
    lines = ["reference\tsample\tplatform\tdepth\tinput\tstratum\tvaf_lo\tvaf_hi\tsites\n"]
    for platform, depth, source, counts in results:
        for name in order:
            for b, n in enumerate(counts[name]):
                lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%.4g\t%.4g\t%.3f\n" % (
                    args.reference, args.sample, platform, depth, source, name,
                    b / args.bins, (b + 1) / args.bins, n))
    tmp = "%s.%d.tmp" % (args.out, os.getpid())
    with open(tmp, "w") as fh:
        fh.writelines(lines)
    os.replace(tmp, args.out)
    print("%s: %d cells, %d strata" % (args.out, len(results), len(order)), file=sys.stderr)


if __name__ == "__main__":
    main()
