#!/usr/bin/env python3
"""CMRG read haplotagging recall against the consensus parental labels, the
reads anchored on T2T-CHM13 v1.0, as the common mode of the grid (bin/cmrg_dev.sh reads).

The reads of a gene G (per platform), on the anchor (T2T-CHM13 v1.0):
  - a consensus label (PAT or MAT);
  - half of the primary alignment or more inside G;
  - tagged by the benchmark VCF (whatshap haplotag, at least one benchmark SNP);
  - a primary alignment in the 30x data set of the other reference too
    (--require: its read table), so both references score the same reads.
A read is correct on a reference when its primary alignment there also lies
half or more inside the same gene (by name) - a read placed elsewhere, e.g.
on a false duplicate of GRCh38, is a miss there - the run tags it, its block
can be oriented and its oriented haplotype is that of its label.

PAT/MAT are put on the benchmark haplotypes per gene and reference, by
majority vote of the benchmark haplotags of the labelled reads in the gene.
A run's block is oriented by its SNPs against the benchmark (more opposite
than agreeing SNPs: flipped; no benchmark SNP: not oriented, a miss).

  python plots/cmrg_reads.py --labels ont5k=reads/ont5k_T2T_HG002_30x_wgs.consensus.tsv.gz \
      --require ont5k=reads/ont5k_GRCh38_HG002_30x_wgs.reads.tsv.gz \
      --anchor CHM13v1=cmrg_dev/CHM13v1/cmrg --other GRCh38=cmrg_dev/GRCh38/cmrg --out cmrg_dev/cmrg_reads
-> <out>.tsv        run reference reads correct recall_pct
   <out>.genes.tsv  run reference gene reads correct recall_pct
"""

import argparse
import glob
import gzip
import os
import sys
from collections import defaultdict

AUTOSOMES = {"chr%d" % i for i in range(1, 23)}
FRACTION = 0.5


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def labels(path):
    out = {}
    with opener(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 1 and f[1] in ("PAT", "MAT"):
                out[f[0]] = f[1]
    return out


def read_names(path):
    """the first column of a read table (primary alignments), header skipped"""
    with opener(path) as fh:
        fh.readline()
        return {line.split("\t", 1)[0] for line in fh}


def genes_bed(path):
    g = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0] in AUTOSOMES:
                g[f[0]].append((int(f[1]), int(f[2]), f[3]))
    return g


def read_genes(aln, genes):
    """read -> genes holding half of its primary alignment or more"""
    out = {}
    with open(aln) as fh:
        for line in fh:
            name, chrom, start, end = line.rstrip("\n").split("\t")
            s, e = int(start), int(end)
            hit = {n for gs, ge, n in genes.get(chrom, ())
                   if min(e, ge) - max(s, gs) >= FRACTION * max(1, e - s)}
            if hit:
                out[name] = hit
    return out


def haplotags(path):
    """whatshap --output-haplotag-list -> read -> (haplotype 1/2, block)"""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[1] in ("H1", "H2"):
                out[f[0]] = (f[1][1], f[2])
    return out


def phased_snps(path, with_ps):
    """(chrom, pos) -> (ref, alt, first allele, block)"""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 10)
            if f[0] not in AUTOSOMES or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            fmt, val = f[8].split(":"), f[9].rstrip("\n").split(":")
            if val[0] not in ("0|1", "1|0"):
                continue
            ps = f[0]
            if with_ps and "PS" in fmt and fmt.index("PS") < len(val):
                ps = val[fmt.index("PS")]
            out[(f[0], int(f[1]))] = (f[3].upper(), f[4].upper(), val[0][0], ps)
    return out


def flips(run_vcf, truth):
    """block -> True when flipped against the benchmark (absent: no benchmark SNP)"""
    n = defaultdict(lambda: [0, 0])
    for k, (r, a, h, ps) in phased_snps(run_vcf, True).items():
        t = truth.get(k)
        if t and t[0] == r and t[1] == a:
            n[ps][0 if h == t[2] else 1] += 1
    return {ps: c[1] > c[0] for ps, c in n.items()}


def reference(d, platform, lab):
    """genes of each read, benchmark tags, PAT -> haplotype 1 per gene, benchmark SNPs"""
    rg = read_genes(os.path.join(d, "reads", platform + ".aln.tsv"), genes_bed(os.path.join(d, "eval", "genes.bed")))
    tt = haplotags(os.path.join(d, "reads", platform + ".truth.tsv"))
    vote = defaultdict(lambda: [0, 0])
    for r, (h, _) in tt.items():
        if r in lab and r in rg:
            for g in rg[r]:
                vote[g][0 if (lab[r] == "PAT") == (h == "1") else 1] += 1
    pat_is_1 = {g: v[0] >= v[1] for g, v in vote.items()}
    truth = phased_snps(os.path.join(d, "truth", "truth.het_snps.phased.vcf.gz"), False)
    return rg, tt, pat_is_1, truth


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--labels", action="append", required=True, help="platform=consensus.tsv.gz of the anchor's cell")
    ap.add_argument("--require", action="append", required=True, help="platform=read table of the other reference's 30x cell")
    ap.add_argument("--anchor", required=True, help="name=cmrg_dev/CHM13v1/<input>")
    ap.add_argument("--other", required=True, help="name=cmrg_dev/GRCh38/<input>")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    refs = [a.anchor.split("=", 1), a.other.split("=", 1)]
    lab_paths = dict(x.split("=", 1) for x in a.labels)
    req_paths = dict(x.split("=", 1) for x in a.require)

    rows, gene_rows = [], []
    for platform in sorted(lab_paths):
        need = [os.path.join(d, "reads", platform + x) for _, d in refs for x in (".aln.tsv", ".truth.tsv")]
        if platform not in req_paths or not all(os.path.exists(f) for f in need):
            print("%s: missing inputs, skipped" % platform, file=sys.stderr)
            continue
        lab = labels(lab_paths[platform])
        on_other = read_names(req_paths[platform])
        per_ref = {name: reference(d, platform, lab) for name, d in refs}
        rg_a, tt_a, _, _ = per_ref[refs[0][0]]
        # the reads of each gene, from the anchor
        gene_of = {r: g for r, g in rg_a.items() if r in lab and r in tt_a and r in on_other}
        print("%s: %d labelled reads in CMRG genes on %s, also on %s" % (platform, len(gene_of), refs[0][0], refs[1][0]),
              file=sys.stderr)
        for name, d in refs:
            rg, _, pat_is_1, truth = per_ref[name]
            for tags_path in sorted(glob.glob(os.path.join(d, "reads", platform + "_*.tags.tsv"))):
                run = os.path.basename(tags_path)[:-len(".tags.tsv")]
                fl = flips(os.path.join(d, "phased", run + ".vcf.gz"), truth)
                tags = haplotags(tags_path)
                n, ok = defaultdict(int), defaultdict(int)
                for r, gs in gene_of.items():
                    here = rg.get(r, set())
                    t = tags.get(r)
                    any_good = 0
                    for g in gs:
                        good = 0
                        if g in here and t and t[1] in fl and g in pat_is_1:
                            h = t[0] if not fl[t[1]] else ("2" if t[0] == "1" else "1")
                            want = "1" if (lab[r] == "PAT") == pat_is_1[g] else "2"
                            good = int(h == want)
                        n[g] += 1
                        ok[g] += good
                        any_good = max(any_good, good)
                    n[""] += 1
                    ok[""] += any_good
                label = platform + "_" + run.split("_")[-1]
                rows.append((label, name, n[""], ok[""]))
                for g in sorted(k for k in n if k):
                    gene_rows.append((label, name, g, n[g], ok[g]))
    with open(a.out + ".tsv", "w") as fh:
        fh.write("run\treference\treads\tcorrect\trecall_pct\n")
        for r in rows:
            fh.write("%s\t%s\t%d\t%d\t%.2f\n" % (r + (100.0 * r[3] / r[2] if r[2] else 0,)))
    with open(a.out + ".genes.tsv", "w") as fh:
        fh.write("run\treference\tgene\treads\tcorrect\trecall_pct\n")
        for r in gene_rows:
            fh.write("%s\t%s\t%s\t%d\t%d\t%.2f\n" % (r + (100.0 * r[4] / r[3] if r[3] else 0,)))
    print("-> %s.tsv, %s.genes.tsv" % (a.out, a.out), file=sys.stderr)


if __name__ == "__main__":
    main()
