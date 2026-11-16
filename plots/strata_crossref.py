#!/usr/bin/env python3
"""Where the reads of a stratum on one reference land on the other: for the
upset plot of plots/figures_upset.py (bin/strata_crossref.sh runs it).

A read of a stratum on the source reference: a consensus label (PAT/MAT), its
primary alignment on chr1-22 half or more inside the stratum cut to the
source's confident regions (the strata of the tables). On the other reference
three yes/no:

  aligned     a primary alignment, MAPQ >= --min-mapq, on the same chromosome
  confident   half or more of it inside the other reference's confident BED
  same        half or more of it inside the same stratum there (its whole BED,
              not cut to the confident regions, so the two questions stay apart)

A read with no primary alignment on chr1-22 there is no to all three. Counted
for every labelled read and for the reads of the recall denominator (>= 2
truth het SNPs on the source).

  python plots/strata_crossref.py --labels L --src-table S --dst-table D --src-alleles A \
      --dst-conf C --stratum name=src.bed=dst.bed ... --direction T2T_to_GRCh38 --platform ont5k --out OUT
appends to OUT: direction platform stratum reads aligned confident same count
"""

import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchio import IntervalSet, is_autosome, open_text


def table(path):
    """read -> (chrom, start, end, mapq), primary alignments on chr1-22"""
    out = {}
    with open_text(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ic, i_s, ie, iq = h.index("chrom"), h.index("start"), h.index("end"), h.index("mapq")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if is_autosome(f[ic]):
                out[f[0]] = (f[ic], int(f[i_s]), int(f[ie]), int(f[iq]))
    return out


def snps(path):
    out = {}
    with open_text(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        a, b = h.index("n_ref"), h.index("n_alt")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = int(f[a]) + int(f[b])
    return out


def inside(iv, c, s, e):
    return iv.overlap(c, s, e) >= 0.5 * max(1, e - s)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for a in ("labels", "src-table", "dst-table", "src-alleles", "dst-conf", "direction", "platform", "out"):
        ap.add_argument("--" + a, required=True)
    ap.add_argument("--stratum", action="append", required=True, help="name=source BED=other reference's BED")
    ap.add_argument("--min-mapq", type=int, default=20)
    a = ap.parse_args()

    labels = set()
    with open_text(a.labels) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 1 and f[1] in ("PAT", "MAT"):
                labels.add(f[0])
    src, dst = table(a.src_table), table(a.dst_table)
    nsnp = snps(a.src_alleles)
    conf = IntervalSet.from_bed(a.dst_conf, autosomes_only=True)
    rows = []
    for item in a.stratum:
        name, src_bed, dst_bed = item.split("=")
        if not (os.path.exists(src_bed) and os.path.exists(dst_bed)):
            print("skip %s: no BED" % name, file=sys.stderr)
            continue
        here = IntervalSet.from_bed(src_bed, autosomes_only=True)
        there = IntervalSet.from_bed(dst_bed, autosomes_only=True)
        n = {"labelled": Counter(), "min2snps": Counter()}
        for r, (c, s, e, _) in src.items():
            if r not in labels or not inside(here, c, s, e):
                continue
            d = dst.get(r)
            key = (int(d is not None and d[3] >= a.min_mapq and d[0] == c),
                   int(d is not None and inside(conf, d[0], d[1], d[2])),
                   int(d is not None and inside(there, d[0], d[1], d[2])))
            n["labelled"][key] += 1
            if nsnp.get(r, 0) >= 2:
                n["min2snps"][key] += 1
        for reads, counts in n.items():
            for (al, co, sa), k in sorted(counts.items()):
                rows.append([a.direction, a.platform, name, reads, str(al), str(co), str(sa), str(k)])
    new = not os.path.exists(a.out)
    with open(a.out, "a") as fh:
        if new:
            fh.write("direction\tplatform\tstratum\treads\taligned\tconfident\tsame\tcount\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print("%s %s: %d strata" % (a.direction, a.platform, len({r[2] for r in rows})), file=sys.stderr)


if __name__ == "__main__":
    main()
