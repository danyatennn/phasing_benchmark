#!/usr/bin/env python3
"""Parental-haplotype truth for reads: which copy of the diploid assembly a
read aligns to better. Reads minimap2 PAF (with -c, so every record has its
alignment score AS) on stdin, records of one read together, as minimap2 writes
them; the assembly holds both copies of every chromosome, named
chrN_MATERNAL / chrN_PATERNAL as in HG002 Q100.

Per read: the best AS on a maternal contig and on a paternal one. The read is
MAT or PAT if that copy wins by at least --min-diff (about one mismatch on
either preset), and none if the copies are as good as each other - the read
sits where the two haplotypes do not differ. A read with alignments to one copy
only takes that copy.

Output (gzip), one line per read: read label chrom as_mat as_pat
  chrom  the chromosome of the best alignment, the copy suffix removed

  minimap2 -c --secondary=yes -N 10 -p 0.5 -x map-hifi q100.mmi reads.fq \\
    | python plots/matpat_assign.py --out reads/<cell>.matpat.tsv.gz
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import write_text

COPY = re.compile(r"^(.*?)[_.#-](MATERNAL|PATERNAL|MAT|PAT)$", re.IGNORECASE)


def copy_of(contig):
    m = COPY.match(contig)
    if m is None:
        return None, contig
    return ("MAT" if m.group(2).upper().startswith("MAT") else "PAT"), m.group(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-diff", type=int, default=10,
                        help="AS margin to call a copy (default %(default)s)")
    args = parser.parse_args()

    out = write_text(args.out)
    counts = {"MAT": 0, "PAT": 0, "none": 0}
    unknown = set()
    current = None
    best = {}

    def flush():
        if current is None:
            return
        mat = best.get("MAT")
        pat = best.get("PAT")
        if mat is None and pat is None:
            return
        if pat is None or (mat is not None and mat[0] - pat[0] >= args.min_diff):
            label = "MAT"
        elif mat is None or pat[0] - mat[0] >= args.min_diff:
            label = "PAT"
        else:
            label = "none"
        top = max((v for v in (mat, pat) if v is not None), key=lambda v: v[0])
        counts[label] += 1
        out.write("%s\t%s\t%s\t%s\t%s\n" % (current, label, top[1],
                                            mat[0] if mat else "NA", pat[0] if pat else "NA"))

    try:
        out.write("#read\tlabel\tchrom\tas_mat\tas_pat\n")
        for line in sys.stdin:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            name = f[0]
            if name != current:
                flush()
                current = name
                best = {}
            copy, chrom = copy_of(f[5])
            if copy is None:
                unknown.add(f[5])
                continue
            score = None
            for tag in f[12:]:
                if tag.startswith("AS:i:"):
                    score = int(tag[5:])
                    break
            if score is None:
                sys.exit("PAF record without AS: run minimap2 with -c")
            if copy not in best or score > best[copy][0]:
                best[copy] = (score, chrom)
        flush()
    except BaseException:
        out.abandon()
        raise
    out.commit()
    total = sum(counts.values())
    print("%s: %d reads, MAT %d, PAT %d, none %d" % (args.out, total, counts["MAT"], counts["PAT"],
                                                       counts["none"]), file=sys.stderr)
    if unknown:
        print("WARNING: %d contigs not named as maternal or paternal, e.g. %s" % (
            len(unknown), sorted(unknown)[:3]), file=sys.stderr)


if __name__ == "__main__":
    main()
