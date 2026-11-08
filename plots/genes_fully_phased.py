#!/usr/bin/env python3
"""Genes fully phased: per run, how many genes of a gene BED one phase block covers.

A phase block is what whatshap stats counts: the phased het genotypes (GT with
'|') of the run's VCF that share a PS, two or more of them; its span runs from
the first to the last of them. Per gene (BED: chrom start end name, 0-based):

  full_span        the gene interval lies inside the span of one block - the
                   "genes fully phased" of HiPhase (Holt et al. 2024)
  all_hets_1block  the gene has >= 2 het genotypes in the run's VCF and all of
                   them are phased, in one block; genes with < 2 hets are
                   left out of its denominator (genes_2het)
  blocks           phase blocks with a het of the gene

No truth is used: the metric depends on the run's VCF alone and the genes are
not cut here. But the callers call inside the confident BED only, so a run's VCF
has no het outside it - a gene partly outside the confident BED
is fully phased only if one block spans the gap.

  python plots/genes_fully_phased.py --genes regions/genes/CMRG/GRCh38_genes.bed \\
      --set CMRG_GRCh38 --refs GRCh38 --depths 30 --out-dir metrics_genes/CMRG/body

<out-dir>/genes_fully_phased_<set>.tsv          one line per run (out-dir: metrics/)
  reference caller platform depth phaser gene_set genes full_span full_span_pct
  genes_2het all_hets_1block all_hets_1block_pct
<out-dir>/genes_fully_phased_<set>.genes.tsv.gz one line per run x gene
  reference caller platform depth phaser gene chrom start end hets phased
  blocks full_span all_hets_1block
"""

import argparse
import bisect
import gzip
import os
import sys
from multiprocessing import Pool

ROOT = "/autofs/bal36/dten/phasing_benchmark"
REGION = "wgs"
CALLERS = ["truth", "clair3", "deepvariant", "longcalld"]
PLATFORMS = ["ont4k", "ont5k", "hifi", "revio"]
PHASERS = ["whatshap", "longphase", "margin", "hapcut2", "hiphase"]
PACBIO_ONLY = {"hiphase", "mphiphase"}


def read_genes(path):
    genes = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            name = f[3] if len(f) > 3 else "%s:%s-%s" % (f[0], f[1], f[2])
            genes.append((f[0], int(f[1]), int(f[2]), name))
    return genes


def het_sites(vcf, chroms):
    """per chromosome: sorted positions and, for each, its PS if phased (else None)"""
    sites = {}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 10)
            if f[0] not in chroms:
                continue
            keys = f[8].split(":")
            if "GT" not in keys:
                continue
            vals = f[9].rstrip("\n").split(":")
            gt = vals[keys.index("GT")]
            sep = "|" if "|" in gt else "/"
            a = gt.split(sep)
            if len(a) != 2 or a[0] == a[1] or "." in a:
                continue
            ps = None
            if sep == "|" and "PS" in keys:
                i = keys.index("PS")
                if i < len(vals) and vals[i] not in (".", ""):
                    ps = vals[i]
            sites.setdefault(f[0], []).append((int(f[1]), ps))
    for c in sites:
        sites[c].sort()
    return sites


def one(job):
    root, sample, ref, caller, platform, depth, phaser, genes = job
    name = "%s_%s_%s_%dx_%s_%s_%s" % (platform, ref, sample, depth, REGION, caller, phaser)
    vcf = os.path.join(root, "phased", name + ".vcf.gz")
    if not os.path.exists(vcf):
        return None, "no VCF   " + name
    sites = het_sites(vcf, {g[0] for g in genes})
    # block spans: (chrom, PS) -> [first, last, n]
    span = {}
    for c, lst in sites.items():
        for pos, ps in lst:
            if ps is None:
                continue
            b = span.get((c, ps))
            if b is None:
                span[(c, ps)] = [pos, pos, 1]
            else:
                b[0] = min(b[0], pos)
                b[1] = max(b[1], pos)
                b[2] += 1
    by_chrom = {}
    for (c, ps), (s, e, n) in span.items():
        if n >= 2:
            by_chrom.setdefault(c, []).append((s, e, ps))
    rows = []
    for c, s0, e0, gname in genes:
        lst = sites.get(c, [])
        pos = [p for p, _ in lst]
        i, j = bisect.bisect_left(pos, s0 + 1), bisect.bisect_right(pos, e0)
        hets = lst[i:j]
        pss = {ps for _, ps in hets if ps is not None}
        # whole gene (1-based s0+1 .. e0) inside one block's span
        full = any(bs <= s0 + 1 and be >= e0 for bs, be, _ in by_chrom.get(c, []))
        one_block = (len(hets) >= 2 and all(ps is not None for _, ps in hets) and len(pss) == 1
                     and span[(c, next(iter(pss)))][2] >= 2)
        rows.append((gname, c, s0, e0, len(hets), sum(ps is not None for _, ps in hets),
                     len(pss), int(full), int(one_block)))
    return (ref, caller, platform, depth, phaser, rows), "done     " + name


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=ROOT)
    ap.add_argument("--sample", default="HG002")
    ap.add_argument("--genes", required=True, help="BED: chrom start end name")
    ap.add_argument("--set", required=True, help="name of the gene set, in the output names")
    ap.add_argument("--refs", default="GRCh38")
    ap.add_argument("--callers", default=" ".join(CALLERS))
    ap.add_argument("--platforms", default=" ".join(PLATFORMS))
    ap.add_argument("--depths", default="30")
    ap.add_argument("--phasers", default=" ".join(PHASERS))
    ap.add_argument("--out-dir", default=None, help="default: <root>/metrics")
    ap.add_argument("-j", type=int, default=8)
    args = ap.parse_args()

    root = args.root if args.sample == "HG002" else os.path.join(args.root, "samples", args.sample)
    genes = read_genes(args.genes)
    if not genes:
        sys.exit("no gene in %s" % args.genes)
    jobs = [(root, args.sample, ref, c, p, int(d), ph, genes)
            for ref in args.refs.split() for c in args.callers.split()
            for p in args.platforms.split() for d in args.depths.split()
            for ph in args.phasers.split()
            if not (ph in PACBIO_ONLY and p not in ("hifi", "revio"))]
    print("%d genes of %s, %d runs" % (len(genes), args.genes, len(jobs)), file=sys.stderr)

    out_dir = args.out_dir or os.path.join(root, "metrics")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "genes_fully_phased_%s.tsv" % args.set)
    per_gene = os.path.join(out_dir, "genes_fully_phased_%s.genes.tsv.gz" % args.set)
    results = []
    with Pool(args.j) as pool:
        for res, msg in pool.imap_unordered(one, jobs):
            print(msg, file=sys.stderr)
            if res is not None:
                results.append(res)
    results.sort(key=lambda r: r[:5])

    tmp = "%s.%d" % (out, os.getpid())
    with open(tmp, "w") as fh:
        fh.write("reference\tcaller\tplatform\tdepth\tphaser\tgene_set\tgenes\tfull_span\tfull_span_pct"
                 "\tgenes_2het\tall_hets_1block\tall_hets_1block_pct\n")
        for ref, c, p, d, ph, rows in results:
            n = len(rows)
            full = sum(r[7] for r in rows)
            n2 = sum(r[4] >= 2 for r in rows)
            ob = sum(r[8] for r in rows)
            fh.write("%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n" % (
                ref, c, p, d, ph, args.set, n, full, 100.0 * full / n, n2, ob,
                100.0 * ob / n2 if n2 else float("nan")))
    os.replace(tmp, out)

    tmp = "%s.%d" % (per_gene, os.getpid())
    with gzip.open(tmp, "wt") as fh:
        fh.write("reference\tcaller\tplatform\tdepth\tphaser\tgene\tchrom\tstart\tend\thets\tphased"
                 "\tblocks\tfull_span\tall_hets_1block\n")
        for ref, c, p, d, ph, rows in results:
            for g, ch, s, e, h, ph_n, b, full, ob in rows:
                fh.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (
                    ref, c, p, d, ph, g, ch, s, e, h, ph_n, b, full, ob))
    os.replace(tmp, per_gene)
    print("genes fully phased: %s, %s" % (out, per_gene), file=sys.stderr)


if __name__ == "__main__":
    main()
