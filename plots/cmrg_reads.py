#!/usr/bin/env python3
"""CMRG read haplotagging recall against the consensus parental labels, the
rules of the genomic regions of the grid (plots/read_recall_strata.py, common
mode anchored on T2T-CHM13), with T2T-CHM13 v1.0 as the anchor and the CMRG
benchmark BED where the grid has the confident regions (bin/cmrg_dev.sh reads).

The reads of a region, on the anchor (T2T-CHM13 v1.0), the same for both references:
  - a consensus label (PAT or MAT);
  - a primary alignment on the other reference too (its 30x read table);
  - 2 or more benchmark het SNPs covered with REF or ALT (read_snp_count.py);
  - half of the primary alignment or more inside the region: a gene cut to
    the CMRG BED, or all of them ("CMRG").
A read is correct on a reference when the run tags it there (whatshap
haplotag, anywhere in the CMRG regions of that reference) in a phase block of
2 or more phased het SNPs, the block is oriented by the majority of the
benchmark SNPs it shares, and its oriented haplotype is that of its label.
PAT/MAT are put on the benchmark haplotypes per chromosome, by majority vote of
the benchmark haplotags of the labelled reads on that reference (the vote and
its agreement are printed).

  python plots/cmrg_reads.py --labels ont5k=reads/ont5k_T2T_HG002_30x_wgs.consensus.tsv.gz \
      --require ont5k=reads/ont5k_GRCh38_HG002_30x_wgs.reads.tsv.gz \
      --anchor CHM13v1=cmrg_dev/CHM13v1/cmrg --anchor-bed CMRG_v1.bed \
      --other GRCh38=cmrg_dev/GRCh38/cmrg --out cmrg_dev/cmrg_reads_cmrg
-> <out>.tsv   run reference region reads correct recall_pct  (region: CMRG, then each gene
              with --min-gene-reads reads or more)

Every gene is counted. --exclude can leave genes out for a check, e.g.
A4GALT (the benchmark's H1 is the maternal haplotype there while it is the
paternal one on the rest of chr22: none of the labelled reads agrees) and DAZL
(T2T-CHM13 v1.0 has no chrY, so the reads of the DAZ genes of chrY align to it;
on GRCh38 they go to chrY).
"""

import argparse
import glob
import gzip
import os
import sys
from collections import defaultdict

AUTOSOMES = {"chr%d" % i for i in range(1, 23)}
FRACTION = 0.5
MIN_READ_SNPS = 2
MIN_BLOCK_SNPS = 2
EXCLUDE = ""


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
    with opener(path) as fh:
        fh.readline()
        return {line.split("\t", 1)[0] for line in fh}


def bed(path, named=False):
    out = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] in AUTOSOMES:
                out[f[0]].append((int(f[1]), int(f[2]), f[3] if named and len(f) > 3 else ""))
    return out


def cut(genes, regions):
    """each gene cut to the regions: gene -> chrom -> [(start, end)]"""
    out = defaultdict(lambda: defaultdict(list))
    for c, gs in genes.items():
        for gs_, ge, name in gs:
            for rs, re_, _ in regions.get(c, ()):
                s, e = max(gs_, rs), min(ge, re_)
                if s < e:
                    out[name][c].append((s, e))
    return out


def merged(regions):
    """chrom -> merged intervals"""
    out = {}
    for c, iv in regions.items():
        iv = sorted(iv)
        m = []
        for s, e in iv:
            if m and s <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        out[c] = m
    return out


def overlap(iv, s, e):
    return sum(max(0, min(e, ie) - max(s, is_)) for is_, ie in iv)


def alignments(path):
    """read -> (chrom, start, end)"""
    out = {}
    with open(path) as fh:
        for line in fh:
            name, chrom, start, end = line.rstrip("\n").split("\t")
            out[name] = (chrom, int(start), int(end))
    return out


def two_col(path, cast=str):
    out = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = cast(f[1])
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


def orientation(run_vcf, truth):
    """block -> flipped against the benchmark? (blocks of >= 2 phased SNPs that
    share a benchmark SNP; any other block cannot be oriented)"""
    size, n = defaultdict(int), defaultdict(lambda: [0, 0])
    for k, (r, a, h, ps) in phased_snps(run_vcf, True).items():
        size[ps] += 1
        t = truth.get(k)
        if t and t[0] == r and t[1] == a:
            n[ps][0 if h == t[2] else 1] += 1
    return {ps: c[1] > c[0] for ps, c in n.items() if size[ps] >= MIN_BLOCK_SNPS}


def pat_haplotype(aln, ttags, lab):
    """chromosome -> True when PAT is haplotype 1 of the benchmark there; prints the vote"""
    vote = defaultdict(lambda: [0, 0])
    for r, (h, _) in ttags.items():
        if r in lab and r in aln:
            vote[aln[r][0]][0 if (lab[r] == "PAT") == (h == "1") else 1] += 1
    out = {}
    for c, (same, opp) in sorted(vote.items()):
        out[c] = same >= opp
        if same + opp:
            print("    %s: PAT = H%s, %.1f%% of %d labelled reads agree"
                  % (c, "1" if out[c] else "2", 100.0 * max(same, opp) / (same + opp), same + opp), file=sys.stderr)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--labels", action="append", required=True, help="platform=consensus labels of the anchor's cell")
    ap.add_argument("--require", action="append", required=True, help="platform=read table of the other reference")
    ap.add_argument("--anchor", required=True, help="name=cmrg_dev/CHM13v1/<input>")
    ap.add_argument("--anchor-bed", required=True, help="the CMRG benchmark BED of the anchor")
    ap.add_argument("--other", required=True, help="name=cmrg_dev/GRCh38/<input>")
    ap.add_argument("--out", required=True)
    ap.add_argument("--exclude", default=EXCLUDE, help="genes left out, comma-separated (default none)")
    ap.add_argument("--min-gene-reads", type=int, default=10,
                    help="a gene's own row only with this many reads or more (CMRG counts every read)")
    a = ap.parse_args()
    refs = [a.anchor.split("=", 1), a.other.split("=", 1)]
    anchor_dir = refs[0][1]
    regions = cut(bed(os.path.join(anchor_dir, "eval", "genes.bed"), named=True), bed(a.anchor_bed))
    left_out = {g for g in a.exclude.split(",") if g}
    print("left out: %s" % ", ".join(sorted(g for g in left_out if g in regions)), file=sys.stderr)
    regions = {g: iv for g, iv in regions.items() if g not in left_out}
    union = merged({c: [iv for g in regions.values() for iv in g.get(c, [])] for c in AUTOSOMES})
    lab_paths = dict(x.split("=", 1) for x in a.labels)
    req_paths = dict(x.split("=", 1) for x in a.require)

    rows = []
    for platform in sorted(lab_paths):
        if platform not in req_paths:
            continue
        need = [os.path.join(d, "reads", platform + x) for _, d in refs for x in (".aln.tsv", ".truth.tsv")]
        need.append(os.path.join(anchor_dir, "reads", platform + ".snps.tsv"))
        if not all(os.path.exists(f) for f in need):
            print("%s: missing inputs of the tags step, skipped" % platform, file=sys.stderr)
            continue
        lab = labels(lab_paths[platform])
        on_other = read_names(req_paths[platform])
        aln_a = alignments(os.path.join(anchor_dir, "reads", platform + ".aln.tsv"))
        nsnp = two_col(os.path.join(anchor_dir, "reads", platform + ".snps.tsv"), int)
        # the reads of each region, on the anchor
        members = defaultdict(set)
        for r, (c, s, e) in aln_a.items():
            if r not in lab or r not in on_other or nsnp.get(r, 0) < MIN_READ_SNPS:
                continue
            need_bp = FRACTION * max(1, e - s)
            if overlap(union.get(c, []), s, e) >= need_bp:
                members["CMRG"].add(r)
            for g, iv in regions.items():
                if overlap(iv.get(c, []), s, e) >= need_bp:
                    members[g].add(r)
        shown = sorted(g for g in members if g != "CMRG" and len(members[g]) >= a.min_gene_reads)
        print("%s: %d labelled reads with >= %d benchmark SNPs in the CMRG genes on %s, also on %s; "
              "%d genes with >= %d reads"
              % (platform, len(members["CMRG"]), MIN_READ_SNPS, refs[0][0], refs[1][0],
                 len(shown), a.min_gene_reads), file=sys.stderr)
        for name, d in refs:
            aln = alignments(os.path.join(d, "reads", platform + ".aln.tsv"))
            print("  %s %s, PAT/MAT on the benchmark haplotypes:" % (platform, name), file=sys.stderr)
            pat1 = pat_haplotype(aln, haplotags(os.path.join(d, "reads", platform + ".truth.tsv")), lab)
            truth = phased_snps(os.path.join(d, "truth", "truth.het_snps.phased.vcf.gz"), False)
            for tags_path in sorted(glob.glob(os.path.join(d, "reads", platform + "_*.tags.tsv"))):
                run = os.path.basename(tags_path)[:-len(".tags.tsv")]
                flip = orientation(os.path.join(d, "phased", run + ".vcf.gz"), truth)
                tags = haplotags(tags_path)
                correct = set()
                for r, (h, ps) in tags.items():
                    if r not in lab or ps not in flip or r not in aln or aln[r][0] not in pat1:
                        continue
                    hh = h if not flip[ps] else ("2" if h == "1" else "1")
                    want = "1" if (lab[r] == "PAT") == pat1[aln[r][0]] else "2"
                    if hh == want:
                        correct.add(r)
                label = platform + "_" + run.split("_")[-1]
                for region in ["CMRG"] + shown:
                    m = members[region]
                    rows.append((label, name, region, len(m), len(m & correct)))
    with open(a.out + ".tsv", "w") as fh:
        fh.write("run\treference\tregion\treads\tcorrect\trecall_pct\n")
        for r in rows:
            fh.write("%s\t%s\t%s\t%d\t%d\t%.2f\n" % (r + (100.0 * r[4] / r[3] if r[3] else 0,)))
    print("-> %s.tsv" % a.out, file=sys.stderr)


if __name__ == "__main__":
    main()
