#!/usr/bin/env python3
"""Gene regions on GRCh38 and T2T from one gene list, matched by NCBI GeneID.

The gene names (4th column of --names, e.g. the CMRG gene BED) are looked up in
the RefSeq annotation of GRCh38 (gene=/Name=, then gene_synonym=) to get their
NCBI GeneID; the gene is then taken by that GeneID from the RefSeq annotation
of each reference (gene and pseudogene features, chr1-22). Use the same
annotation release for both (NCBI RS_2025_08: GRCh38.p14 and T2T-CHM13v2.0).
A name with no GeneID or with several is reported and left out; so is a gene
missing on either reference.

Writes to --out:
  genes.tsv                one line per name: status, GeneID, biotype, and per
                           reference the placements, gene bp, bp and % inside
                           the confident BED
  <ref>_genes.bed          chrom start end name, the genes found on both
                           references (a gene placed twice has two lines)
  strata.tsv               tier category stratum, as regions/strata.tsv:
                           tier 1 the whole set (--set), tier 2 one per gene
  body/prepared/<ref>_<stratum>.bed   the whole gene(s)
  conf/prepared/<ref>_<stratum>.bed   the gene(s) cut to the confident BED;
                                      no file if nothing is left

  python plots/gene_regions.py --set CMRG --names GRCh38_CMRG_benchmark_gene_coordinates.bed \\
      --gff GRCh38=GRCh38.gff.gz --gff T2T=T2T.gff.gz \\
      --assembly-report GRCh38=GRCh38_assembly_report.txt --assembly-report T2T=T2T_assembly_report.txt \\
      --confident GRCh38=GRCh38_v5q.bed --confident T2T=CHM13_v5q.bed --out regions/genes/CMRG
"""

import argparse
import gzip
import os
import sys
from urllib.parse import unquote

REFS = ["GRCh38", "T2T"]
AUTOSOMES = {"chr%d" % i for i in range(1, 23)}
GENE_TYPES = {"gene", "pseudogene"}


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def pairs(items, what):
    """["GRCh38=a", "T2T=b"] -> {"GRCh38": "a", "T2T": "b"}, both references required"""
    out = dict(item.split("=", 1) for item in items)
    if sorted(out) != sorted(REFS):
        sys.exit("--%s needs GRCh38=... and T2T=..." % what)
    return out


def chrom_names(report):
    """RefSeq accession -> chr1..chr22, from the NCBI assembly report"""
    names = {}
    with open(report) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            # Sequence-Name Sequence-Role Assigned-Molecule ... RefSeq-Accn (7th)
            if f[1] == "assembled-molecule" and "chr" + f[2] in AUTOSOMES:
                names[f[6]] = "chr" + f[2]
    return names


def read_genes(gff, chroms):
    """the gene features of a RefSeq gff: list of (geneid, symbol, synonyms, biotype, chrom, start, end);
    chrom is None off chr1-22 (alt loci, patches, chrX/Y, ...)"""
    genes = []
    with open_text(gff) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in GENE_TYPES:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            geneid = None
            for x in a.get("Dbxref", "").split(","):
                if x.startswith("GeneID:"):
                    geneid = x[len("GeneID:"):]
            if geneid is None:
                continue
            symbol = unquote(a.get("gene", a.get("Name", "")))
            synonyms = [unquote(s) for s in a.get("gene_synonym", "").split(",") if s]
            # gff is 1-based, closed; BED 0-based, half-open
            genes.append((geneid, symbol, synonyms, a.get("gene_biotype", ""),
                          chroms.get(f[0]), int(f[3]) - 1, int(f[4])))
    return genes


def read_bed(path):
    """chrom -> sorted, merged [start, end) intervals"""
    per = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            per.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return {c: merge(iv) for c, iv in per.items()}


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def cut(intervals, chrom, bed):
    """the parts of [(s, e)] on chrom that lie inside bed (read_bed)"""
    out = []
    for s, e in intervals:
        for bs, be in bed.get(chrom, []):
            if be <= s:
                continue
            if bs >= e:
                break
            out.append((max(s, bs), min(e, be)))
    return out


def write_bed(path, rows):
    with open(path, "w") as fh:
        for c, s, e in sorted(rows, key=lambda r: (int(r[0][3:]), r[1])):
            fh.write("%s\t%d\t%d\n" % (c, s, e))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--set", required=True, help="name of the gene set: the tier-1 stratum")
    ap.add_argument("--names", required=True, help="BED with the gene names in column 4")
    ap.add_argument("--gff", action="append", default=[], help="REF=RefSeq gff(.gz), for GRCh38 and T2T")
    ap.add_argument("--assembly-report", action="append", default=[], help="REF=NCBI assembly report")
    ap.add_argument("--confident", action="append", default=[], help="REF=confident BED")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    gffs = pairs(args.gff, "gff")
    reports = pairs(args.assembly_report, "assembly-report")
    confident = pairs(args.confident, "confident")

    names = []
    with open_text(args.names) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 3 and not line.startswith(("#", "track", "browser")) and f[3] not in names:
                names.append(f[3])

    genes = {ref: read_genes(gffs[ref], chrom_names(reports[ref])) for ref in REFS}
    conf = {ref: read_bed(confident[ref]) for ref in REFS}

    # name -> GeneID, from GRCh38: the official symbol first, then the synonyms
    by_symbol, by_synonym, biotype = {}, {}, {}
    for geneid, symbol, synonyms, bt, _, _, _ in genes["GRCh38"]:
        by_symbol.setdefault(symbol, set()).add(geneid)
        for s in synonyms:
            by_synonym.setdefault(s, set()).add(geneid)
        biotype[geneid] = bt
    # GeneID -> {ref: [(chrom, start, end)]} on chr1-22
    places = {}
    for ref in REFS:
        for geneid, _, _, _, chrom, s, e in genes[ref]:
            if chrom is not None:
                lst = places.setdefault(geneid, {r: [] for r in REFS})[ref]
                if (chrom, s, e) not in lst:
                    lst.append((chrom, s, e))

    os.makedirs(args.out, exist_ok=True)
    report = []
    found = []          # (name, {ref: [(chrom, s, e)]})
    for name in names:
        ids = by_symbol.get(name) or by_synonym.get(name) or set()
        how = "symbol" if name in by_symbol else ("synonym" if ids else "")
        row = {"name": name, "geneid": ",".join(sorted(ids)), "matched_by": how, "biotype": ""}
        if len(ids) != 1:
            row["status"] = "no_geneid" if not ids else "several_geneids"
        else:
            geneid = next(iter(ids))
            row["biotype"] = biotype.get(geneid, "")
            p = places.get(geneid, {r: [] for r in REFS})
            missing = [r for r in REFS if not p[r]]
            row["status"] = "ok" if not missing else "missing_" + "_".join(missing)
            for ref in REFS:
                iv = p[ref]
                bp = sum(e - s for _, s, e in iv)
                in_conf = sum(e - s for c, s0, e0 in iv for s, e in cut([(s0, e0)], c, conf[ref]))
                row[ref + "_placements"] = len(iv)
                row[ref + "_loci"] = ",".join("%s:%d-%d" % (c, s + 1, e) for c, s, e in iv)
                row[ref + "_bp"] = bp
                row[ref + "_confident_bp"] = in_conf
                row[ref + "_confident_pct"] = "%.2f" % (100.0 * in_conf / bp) if bp else ""
            if not missing:
                found.append((name, p))
        report.append(row)

    cols = ["name", "status", "geneid", "matched_by", "biotype"]
    for ref in REFS:
        cols += [ref + x for x in ("_placements", "_loci", "_bp", "_confident_bp", "_confident_pct")]
    with open(os.path.join(args.out, "genes.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in report:
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")

    with open(os.path.join(args.out, "strata.tsv"), "w") as fh:
        fh.write("tier\tcategory\tstratum\n1\tgenes\t%s\n" % args.set)
        for name, _ in found:
            fh.write("2\tgene\t%s\n" % name)

    for ref in REFS:
        with open(os.path.join(args.out, "%s_genes.bed" % ref), "w") as fh:
            rows = sorted((c, s, e, n) for n, p in found for c, s, e in p[ref])
            for c, s, e, n in sorted(rows, key=lambda r: (int(r[0][3:]), r[1])):
                fh.write("%s\t%d\t%d\t%s\n" % (c, s, e, n))
        for mode in ("body", "conf"):
            d = os.path.join(args.out, mode, "prepared")
            os.makedirs(d, exist_ok=True)
            union = []
            for name, p in found:
                iv = []
                for c in sorted({c for c, _, _ in p[ref]}):
                    piece = merge([(s, e) for cc, s, e in p[ref] if cc == c])
                    if mode == "conf":
                        piece = cut(piece, c, conf[ref])
                    iv += [(c, s, e) for s, e in piece]
                if iv:
                    write_bed(os.path.join(d, "%s_%s.bed" % (ref, name)), iv)
                union += iv
            merged = []
            for c in sorted({c for c, _, _ in union}):
                merged += [(c, s, e) for s, e in merge([(s, e) for cc, s, e in union if cc == c])]
            write_bed(os.path.join(d, "%s_%s.bed" % (ref, args.set)), merged)

    n = {}
    for row in report:
        n[row["status"]] = n.get(row["status"], 0) + 1
    print("%d names: %s" % (len(names), ", ".join("%s %d" % kv for kv in sorted(n.items()))), file=sys.stderr)
    for row in report:
        if row["status"] != "ok" or row["matched_by"] == "synonym" \
                or any(row.get(r + "_placements", 1) > 1 for r in REFS):
            print("  %-12s %-16s %s %s" % (row["name"], row["status"], row["matched_by"],
                                           " ".join("%s:%s" % (r, row.get(r + "_loci", "")) for r in REFS)),
                  file=sys.stderr)
    for ref in REFS:
        bp = sum(row.get(ref + "_bp", 0) for row in report if row["status"] == "ok")
        cb = sum(row.get(ref + "_confident_bp", 0) for row in report if row["status"] == "ok")
        print("%s: %d genes, %d bp, %.1f%% inside the confident BED"
              % (ref, len(found), bp, 100.0 * cb / bp if bp else 0), file=sys.stderr)


if __name__ == "__main__":
    main()
