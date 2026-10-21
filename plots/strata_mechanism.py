#!/usr/bin/env python3
"""Why a stratum is hard, in numbers - per cell, and per run for the errors.

Per cell (caller and phaser "-"), from the truth's allele matrix and mosdepth:
  truth_snps           phased truth het SNPs in the stratum
  snp_share            their share of all of them on chr1-22
  vaf_dev              mean |VAF - 0.5| over those SNPs: 0 for a clean het,
                       larger where reads of another copy pile up
  paralog_link_frac    neighbouring-SNP pairs (both in the stratum) showing
                       more than two allele combinations
  depth_ratio          mean depth in the stratum / in the whole confident BED
                       on chr1-22; above 1 for collapsed copies
  lowmapq_frac         1 - depth with MAPQ >= 20 / depth with every read

Per run (caller x phaser):
  switches, switch_share     switch errors whose midpoint lies in the stratum,
                             and their share of all the run's switch errors
  switch_enrichment          switch_share / snp_share: above 1, the stratum
                             takes more errors than its SNPs would explain
  breaks, break_share, break_enrichment
                             the same for phase block ends: the midpoint
                             between the last SNP of a block and the first of
                             the next (blocks of --min-snps or more)

Every count is over chr1-22. Output, long format:
  reference sample platform depth caller phaser stratum metric value

  python plots/strata_mechanism.py --reference GRCh38 --sample HG002 \\
      --truth-vcf truth.phased.vcf.gz --stratum segdups=bed ... --jobs jobs.tsv --out mechanism.tsv
jobs.tsv lines:
  cell  platform depth sites.tsv.gz links.tsv.gz coverage.tsv
  run   platform depth caller phaser switches.bed phased.vcf.gz
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import IntervalSet, is_autosome, open_text, phase_blocks, read_tsv_header, vcf_records


def cell_metrics(strata, truth_pos, total_snps, sites_path, links_path, coverage_path):
    out = []
    vaf = {n: [0.0, 0] for n, _ in strata}
    if sites_path and os.path.exists(sites_path):
        with open_text(sites_path) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[9] == "nan" or not is_autosome(f[0]):
                    continue
                pos0 = int(f[1]) - 1
                d = abs(float(f[9]) - 0.5)
                for n, s in strata:
                    if s.contains(f[0], pos0):
                        vaf[n][0] += d
                        vaf[n][1] += 1
    par = {n: [0, 0] for n, _ in strata}
    if links_path and os.path.exists(links_path):
        with open_text(links_path) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if not is_autosome(f[0]):
                    continue
                a, b, p = int(f[1]) - 1, int(f[2]) - 1, int(f[9])
                for n, s in strata:
                    if s.contains(f[0], a) and s.contains(f[0], b):
                        par[n][0] += p
                        par[n][1] += 1
    cov = {}
    if coverage_path and os.path.exists(coverage_path):
        for r in read_tsv_header(coverage_path):
            cov[r["stratum"]] = (float(r["depth_mapq0"]), float(r["depth_mapq20"]))
    genome_depth = cov.get("confident", (float("nan"), float("nan")))[0]

    for n, s in strata:
        snps = truth_pos[n]
        out.append((n, "truth_snps", snps))
        out.append((n, "snp_share", snps / total_snps if total_snps else float("nan")))
        out.append((n, "vaf_dev", vaf[n][0] / vaf[n][1] if vaf[n][1] else float("nan")))
        out.append((n, "paralog_link_frac", par[n][0] / par[n][1] if par[n][1] else float("nan")))
        if n in cov:
            d0, d20 = cov[n]
            out.append((n, "depth_ratio", d0 / genome_depth if genome_depth else float("nan")))
            out.append((n, "lowmapq_frac", 1 - d20 / d0 if d0 else float("nan")))
    return out


def run_metrics(strata, truth_share, switches_path, vcf, min_snps):
    out = []
    # switch errors, by the midpoint of the interval between their two SNPs
    sw = {n: 0 for n, _ in strata}
    n_sw = 0
    if switches_path and os.path.exists(switches_path):
        with open_text(switches_path) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) < 3 or not is_autosome(f[0]):
                    continue
                n_sw += 1
                mid = (int(f[1]) + int(f[2])) // 2
                for n, s in strata:
                    if s.contains(f[0], mid):
                        sw[n] += 1
    # block ends
    br = {n: 0 for n, _ in strata}
    n_br = 0
    if vcf and os.path.exists(vcf):
        per = {}
        for (chrom, ps), positions in phase_blocks(vcf, autosomes_only=True).items():
            if len(positions) >= min_snps:
                per.setdefault(chrom, []).append((positions[0], positions[-1]))
        for chrom, blocks in per.items():
            blocks.sort()
            reach = blocks[0][1]
            for first, last in blocks[1:]:
                if first > reach:          # a real gap, not a block nested in the previous one
                    n_br += 1
                    mid = (reach + first) // 2 - 1
                    for n, s in strata:
                        if s.contains(chrom, mid):
                            br[n] += 1
                reach = max(reach, last)
    out.append(("genome", "switches", n_sw))
    out.append(("genome", "breaks", n_br))
    for n, _ in strata:
        share = truth_share[n]
        s_share = sw[n] / n_sw if n_sw else float("nan")
        b_share = br[n] / n_br if n_br else float("nan")
        out.append((n, "switches", sw[n]))
        out.append((n, "switch_share", s_share))
        out.append((n, "switch_enrichment", s_share / share if share else float("nan")))
        out.append((n, "breaks", br[n]))
        out.append((n, "break_share", b_share))
        out.append((n, "break_enrichment", b_share / share if share else float("nan")))
    return out


def fmt(v):
    if isinstance(v, int):
        return str(v)
    return "nan" if v != v else "%.6g" % v


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reference", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--truth-vcf", required=True)
    parser.add_argument("--stratum", action="append", default=[], help="name=bed")
    parser.add_argument("--min-snps", type=int, default=2)
    parser.add_argument("--jobs", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    strata = []
    for item in args.stratum:
        name, _, path = item.partition("=")
        if os.path.exists(path):
            strata.append((name, IntervalSet.from_bed(path, autosomes_only=True)))
    truth_pos = {n: 0 for n, _ in strata}
    total = 0
    for c, p, ref, alt, gt, ps in vcf_records(args.truth_vcf, autosomes_only=True):
        total += 1
        for n, s in strata:
            if s.contains(c, p - 1):
                truth_pos[n] += 1
    truth_share = {n: truth_pos[n] / total if total else float("nan") for n, _ in strata}

    lines = ["reference\tsample\tplatform\tdepth\tcaller\tphaser\tstratum\tmetric\tvalue\n"]
    with open(args.jobs) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "cell":
                platform, depth, sites, links, coverage = f[1:6]
                for n, metric, value in cell_metrics(strata, truth_pos, total, sites, links, coverage):
                    lines.append("\t".join([args.reference, args.sample, platform, depth, "-", "-",
                                            n, metric, fmt(value)]) + "\n")
            elif f[0] == "run":
                platform, depth, caller, phaser, switches, vcf = f[1:7]
                for n, metric, value in run_metrics(strata, truth_share, switches, vcf, args.min_snps):
                    lines.append("\t".join([args.reference, args.sample, platform, depth, caller, phaser,
                                            n, metric, fmt(value)]) + "\n")
    tmp = "%s.%d.tmp" % (args.out, os.getpid())
    with open(tmp, "w") as fh:
        fh.writelines(lines)
    os.replace(tmp, args.out)
    print("%s: %d lines" % (args.out, len(lines) - 1), file=sys.stderr)


if __name__ == "__main__":
    main()
