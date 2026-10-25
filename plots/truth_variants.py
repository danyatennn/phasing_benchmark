#!/usr/bin/env python3
"""The inputs of experiment A: the truth het SNPs of the main grid, unchanged,
plus extra heterozygous variants (bin/make_truth_variants.sh normalises them
first). Writes an unsorted VCF; the caller sorts it once.

Kept from --extra: phased heterozygous records (0|1 or 1|0) whose ALT is a
base sequence. Dropped:
  - a record equal to one of the truth SNPs (chrom, pos, ref, alt): the SNP
    set must be exactly that of the main grid
  - a record that occurs twice: after a 1/2 site is atomized, a part shared by
    both haplotypes comes out once per haplotype - it is homozygous, not het
  - '*' or '.' alleles, and genotypes that are not het after atomizing

  python plots/truth_variants.py --snps truth.het_snps.phased.vcf.gz \\
      --extra extra.norm.vcf.gz --sample HG002 --out out.vcf

--sample names the one sample column as in the truth SNPs (the source VCF may
call it NA12878 or HG002 differently).
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import open_text

BASES = re.compile(r"^[ACGTN]+$")


def key_of(f):
    return f[0], f[1], f[3].upper(), f[4].upper()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--snps", required=True)
    parser.add_argument("--extra", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    snps, snp_lines = set(), []
    with open_text(args.snps) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            snps.add(key_of(f))
            snp_lines.append(line)

    header, extra, seen = [], [], {}
    dropped = {"equal to a truth SNP": 0, "not a phased het": 0, "no base allele": 0}
    with open_text(args.extra) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    line = "\t".join(line.rstrip("\n").split("\t")[:9] + [args.sample]) + "\n"
                header.append(line)
                continue
            f = line.rstrip("\n").split("\t")
            gt = f[9].split(":")[0] if len(f) > 9 else "."
            if gt not in ("0|1", "1|0"):
                dropped["not a phased het"] += 1
                continue
            if not BASES.match(f[3].upper()) or not BASES.match(f[4].upper()):
                dropped["no base allele"] += 1
                continue
            k = key_of(f)
            if k in snps:
                dropped["equal to a truth SNP"] += 1
                continue
            seen[k] = seen.get(k, 0) + 1
            # bcftools norm takes the bases it shifts from the reference, and
            # CHM13 is soft-masked; whatshap matches alleles case by case
            f[3], f[4] = f[3].upper(), f[4].upper()
            extra.append((k, "\t".join(f) + "\n"))

    twice = {k for k, n in seen.items() if n > 1}
    kept = [line for k, line in extra if k not in twice]
    dropped["twice (homozygous part)"] = len(extra) - len(kept)

    with open(args.out, "w") as out:
        out.writelines(header)
        out.writelines(snp_lines)
        out.writelines(kept)
    print("%s: %d truth SNPs + %d extra variants; dropped %s" % (
        args.out, len(snp_lines), len(kept),
        ", ".join("%s %d" % (k, v) for k, v in dropped.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
