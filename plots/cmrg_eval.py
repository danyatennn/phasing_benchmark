#!/usr/bin/env python3
"""CMRG phasing per gene, against the CMRG small-variant VCF as released (no
confident BED anywhere).

Truth: the het biallelic SNPs of the CMRG VCF inside the CMRG BED, with their
phase. Per run (one phased VCF) and gene (gene BED: chrom start end name):

  truth_snps     truth het SNPs of the gene (inside the CMRG BED)
  phased         of them, phased in the run (GT with '|', same REF/ALT)
  blocks         phase blocks (PS) among the phased ones
  switches       switch errors inside blocks (a change of orientation between
                 two neighbouring phased SNPs of one block)
  hamming        per block the smaller of agreeing / opposite SNPs, summed
  all_phased     every truth SNP of the gene phased (genes with >= 2 truth SNPs)
  one_block      ... and all in one block
  correct        ... and no switch and no Hamming error: the gene phased
                 correctly end to end

  python plots/cmrg_eval.py --truth CMRG.vcf.gz --bed CMRG.bed --genes genes.bed \
      --run ont5k_whatshap=phased.vcf.gz ... --out cmrg_<ref>
-> <out>.genes.tsv (run x gene) and <out>.tsv (one line per run)
"""

import argparse
import gzip
import sys
from collections import defaultdict

AUTOSOMES = {"chr%d" % i for i in range(1, 23)}


def read_bed(path):
    out = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out[f[0]].append((int(f[1]), int(f[2]), f[3] if len(f) > 3 else ""))
    for v in out.values():
        v.sort()
    return out


def inside(bed, chrom, pos0):
    for s, e, _ in bed.get(chrom, ()):
        if s <= pos0 < e:
            return True
        if s > pos0:
            break
    return False


def het_snps(path, keep=None):
    """(chrom, pos) -> (ref, alt, gt as (a, b) or None if unphased, ps)"""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] not in AUTOSOMES or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            pos = int(f[1])
            if keep is not None and not keep(f[0], pos - 1):
                continue
            fmt = f[8].split(":")
            val = dict(zip(fmt, f[9].split(":")))
            gt = val.get("GT", ".")
            alleles = gt.replace("|", "/").split("/")
            if sorted(alleles) != ["0", "1"]:
                continue
            phased = "|" in gt
            ps = val.get("PS", "") if phased else ""
            if phased and ps in ("", "."):
                ps = f[0]          # phased without PS: one block per chromosome (dipcall-style)
            out[(f[0], pos)] = (f[3].upper(), f[4].upper(), tuple(alleles) if phased else None, ps)
    return out


def gene_stats(truth, run, chrom, start, end):
    snps = sorted(p for (c, p) in truth if c == chrom and start < p <= end)
    phased = []
    for p in snps:
        r = run.get((chrom, p))
        t = truth[(chrom, p)]
        if r and r[2] and r[0] == t[0] and r[1] == t[1]:
            phased.append((p, r[3], r[2] == t[2]))     # same orientation as the truth?
    by_block = defaultdict(list)
    for p, ps, same in phased:
        by_block[ps].append(same)
    switches = sum(sum(a != b for a, b in zip(v, v[1:])) for v in by_block.values())
    hamming = sum(min(sum(v), len(v) - sum(v)) for v in by_block.values())
    n = len(snps)
    all_phased = n >= 2 and len(phased) == n
    one_block = all_phased and len(by_block) == 1
    return {"truth_snps": n, "phased": len(phased), "blocks": len(by_block),
            "switches": switches, "hamming": hamming, "all_phased": int(all_phased),
            "one_block": int(one_block), "correct": int(one_block and switches == 0 and hamming == 0)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--run", action="append", default=[], help="label=phased.vcf.gz, repeated")
    ap.add_argument("--out", required=True)
    ap.add_argument("--conf", help="regions/genes/CMRG/genes.tsv: the share of each gene in the confident BEDs")
    a = ap.parse_args()

    conf = {}
    if a.conf:
        with open(a.conf) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(h, line.rstrip("\n").split("\t")))
                conf[r["name"]] = (r.get("GRCh38_confident_pct", ""), r.get("T2T_confident_pct", ""))

    bed = read_bed(a.bed)
    truth = {k: v for k, v in het_snps(a.truth, lambda c, p: inside(bed, c, p)).items() if v[2]}
    genes = [(c, s, e, n) for c, v in read_bed(a.genes).items() if c in AUTOSOMES for s, e, n in v]
    print("%d truth het SNPs in the CMRG BED, %d genes" % (len(truth), len(genes)), file=sys.stderr)

    cols = ["truth_snps", "phased", "blocks", "switches", "hamming", "all_phased", "one_block", "correct"]
    with open(a.out + ".genes.tsv", "w") as g, open(a.out + ".tsv", "w") as s:
        g.write("run\tgene\tchrom\tstart\tend\t" + "\t".join(cols)
                + "\tGRCh38_confident_pct\tT2T_confident_pct\n")
        s.write("run\tgenes\tgenes_2snps\tall_phased\tone_block\tcorrect\tcorrect_pct"
                "\ttruth_snps\tphased_snps\tphased_pct\tswitches\tswitch_pct\thamming\thamming_pct\n")
        for item in a.run:
            label, path = item.split("=", 1)
            run = het_snps(path, lambda c, p: inside(bed, c, p))
            tot = defaultdict(int)
            for c, st, en, name in genes:
                r = gene_stats(truth, run, c, st, en)
                g.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\n" % (label, name, c, st, en, "\t".join(str(r[k]) for k in cols),
                                                       "\t".join(conf.get(name, ("", "")))))
                for k in cols:
                    tot[k] += r[k]
                tot["genes_2snps"] += r["truth_snps"] >= 2
            pairs = max(1, tot["phased"] - tot["blocks"])
            s.write("%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%d\t%.4f\t%d\t%.4f\n" % (
                label, len(genes), tot["genes_2snps"], tot["all_phased"], tot["one_block"], tot["correct"],
                100.0 * tot["correct"] / max(1, tot["genes_2snps"]), tot["truth_snps"], tot["phased"],
                100.0 * tot["phased"] / max(1, tot["truth_snps"]), tot["switches"], 100.0 * tot["switches"] / pairs,
                tot["hamming"], 100.0 * tot["hamming"] / max(1, tot["phased"])))
            print("%s done" % label, file=sys.stderr)


if __name__ == "__main__":
    main()
