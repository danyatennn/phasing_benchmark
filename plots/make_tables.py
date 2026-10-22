#!/usr/bin/env python3
"""The tables of the paper, from what the pipeline has written: the same runs
and the same numbers as the figures (data.py loads them, the value functions of
figures.py turn them into percent and Mb), so a table never disagrees with a
figure.

  python plots/make_tables.py --inputs tables/inputs.tsv --software tables/software.tsv
  python plots/make_tables.py ... --phasers "whatshap longphase hapcut2"

Run through bin/make_tables.sh, which writes both inputs from config.sh. Every
table is a TSV in --out and a sheet of tables.xlsx; README.tsv says what each
holds and how its numbers are defined. A table with nothing to show yet (strata
before regions.sh, the other samples before their runs) is skipped and the rest
are still written.

  t00_completeness     which runs each table rests on
  t01_datasets         the sequencing data: source depth, read lengths, accuracy
  t02_truth            the truth sets: variant classes inside the benchmark regions
  t03_callers          FNR and FDR of each caller, every cell
  t04_grid             every run of HG002: phasing and read-level metrics
  t05_main_30x         the headline numbers: truth input, 30x
  t06a_reference       GRCh38 against T2T per platform and phaser
  t06b_reference_pooled  the same pooled per phaser and per technology
  t07_ranking          how often each phaser is the best, per metric and input
  t08_coverage         how the metrics change with coverage
  t09_call_vs_phase    what calling and what phasing cost, on the truth denominator
  t10_error_types      switch against flip errors
  t11_reads_30x        read-level metrics at 30x, every input
  t12a_strata          what the strata are
  t12b_strata_switch   switch rate per stratum at 30x
  t13_samples          the other GIAB samples at 30x
  t14_software         tool versions
"""

import argparse
import math
import os
import statistics
import sys
import zipfile
from array import array
from multiprocessing import Pool
from xml.sax.saxutils import escape

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import style as st
from benchio import IntervalSet, is_autosome, open_text, vcf_records
from data import NAN, Metrics, ratio, sample_root
from figures import (n50_mb, read_precision_pct, read_recall_pct, recall_pct, switch_pct,
                     tagged_bases_pct, tagged_pct, truth_phased_pct)

try:
    from scipy.stats import wilcoxon
except ImportError:            # the p-values are left empty without scipy
    wilcoxon = None


def finite(v):
    return isinstance(v, (int, float)) and math.isfinite(v)


def rnd(v, sig=4):
    """Rounded to sig significant figures; NaN and None stay empty."""
    if not finite(v):
        return None
    if isinstance(v, int) or (float(v).is_integer() and abs(v) >= 10 ** sig):
        return int(round(v))
    return float("%.*g" % (sig, v))


def pct(a, b):
    return 100 * ratio(a, b)


class Table:
    def __init__(self, name, title, note, header):
        self.name, self.title, self.note, self.header = name, title, note, header
        self.rows = []

    def add(self, *values):
        assert len(values) == len(self.header), (self.name, len(values), len(self.header))
        self.rows.append([rnd(v) if isinstance(v, float) else v for v in values])


# ------------------------------------------------------------------ writing

def cell_text(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return ("%.6g" % v)
    return str(v)


def write_tsv(table, out):
    path = os.path.join(out, table.name + ".tsv")
    with open(path + ".tmp", "w") as fh:
        fh.write("\t".join(table.header) + "\n")
        for row in table.rows:
            fh.write("\t".join(cell_text(v) for v in row) + "\n")
    os.replace(path + ".tmp", path)


def col_name(i):
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def sheet_xml(header, rows):
    out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
           '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
           '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
           'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews><sheetData>']
    for r, row in enumerate([header] + rows, start=1):
        cells = []
        for c, v in enumerate(row):
            ref = "%s%d" % (col_name(c), r)
            if v is None or v == "":
                continue
            if isinstance(v, (int, float)) and not isinstance(v, bool) and finite(v):
                cells.append('<c r="%s"><v>%s</v></c>' % (ref, repr(v) if isinstance(v, float) else v))
            else:
                style = ' s="1"' if r == 1 else ""
                cells.append('<c r="%s" t="inlineStr"%s><is><t>%s</t></is></c>'
                             % (ref, style, escape(str(v))))
        out.append('<row r="%d">%s</row>' % (r, "".join(cells)))
    out.append("</sheetData></worksheet>")
    return "".join(out)


def write_xlsx(path, sheets):
    """A plain workbook, one sheet per table, header in bold and frozen; no
    library needed (openpyxl is not in every env)."""
    names = []
    for name, _, _ in sheets:
        n = name[:31]
        while n in names:
            n = n[:28] + "_%d" % len(names)
        names.append(n)
    with zipfile.ZipFile(path + ".tmp", "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                   '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                   '<Default Extension="xml" ContentType="application/xml"/>'
                   '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
                   '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
                   + "".join('<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' % (i + 1)
                             for i in range(len(sheets))) + "</Types>")
        z.writestr("_rels/.rels",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
                   "</Relationships>")
        z.writestr("xl/workbook.xml",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                   'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
                   + "".join('<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (escape(n), i + 1, i + 1)
                             for i, n in enumerate(names)) + "</sheets></workbook>")
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   + "".join('<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet%d.xml"/>' % (i + 1, i + 1)
                             for i in range(len(sheets)))
                   + '<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>' % (len(sheets) + 1)
                   + "</Relationships>")
        z.writestr("xl/styles.xml",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                   '<fonts count="2"><font><sz val="11"/><name val="Arial"/></font>'
                   '<font><b/><sz val="11"/><name val="Arial"/></font></fonts>'
                   '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
                   '<borders count="1"><border/></borders><cellStyleXfs count="1"><xf/></cellStyleXfs>'
                   '<cellXfs count="2"><xf fontId="0"/><xf fontId="1" applyFont="1"/></cellXfs></styleSheet>')
        for i, (_, header, rows) in enumerate(sheets):
            z.writestr("xl/worksheets/sheet%d.xml" % (i + 1), sheet_xml(header, rows))
    os.replace(path + ".tmp", path)


# ------------------------------------------------------------------- labels

def plat(p):
    return st.PLATFORM_LABEL.get(p, p)


def phaser_label(p):
    return st.PHASER_LABEL.get(p, p)


def caller_label(c):
    return st.CALLER_LABEL.get(c, c)


def ref_label(r):
    return st.REFERENCE_LABEL.get(r, r)


def anchor_precision_pct(db, key):
    return 100 * db.get(key, "read_precision_anchor_truthtags")


def anchor_recall_pct(db, key):
    return 100 * db.get(key, "read_recall_anchor_truthtags")


def flips(db, key):
    """(switch errors, flip errors) per 100 assessed pairs, chr1-22; a flip is
    two switches in a row around one site and is counted once."""
    rates = db.switch_flip(key)
    return (100 * rates[0], 100 * rates[1]) if rates else (NAN, NAN)


# metric name, value function, True when more is better
METRICS = [
    ("switch_rate_pct", switch_pct, False),
    ("ng50_mb", n50_mb, True),
    ("truth_phased_pct", truth_phased_pct, True),
    ("read_recall_pct", read_recall_pct, True),
]


# ------------------------------------------------------------------ inputs

def read_inputs(path):
    """sample reference item value, from bin/make_tables.sh"""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == 4:
                out.setdefault((f[0], f[1]), {})[f[2]] = f[3]
    return out


# ------------------------------------------------------------- the tables

def t00_completeness(dbs):
    t = Table("t00_completeness", "Which runs the tables rest on",
              "Runs per sample and reference by state: computed (row and chr1-22 metrics there), "
              "no_auto (bin/autosomes.sh not run yet), locked (being written), missing. The other "
              "columns count computed runs that also have NG50 over chr1-22, tagged bases, read "
              "accuracy (SNP-oriented) and read accuracy (anchor-oriented).",
              ["sample", "reference", "computed", "no_auto", "locked", "missing",
               "with_ng50", "with_tagged_bases", "with_read_accuracy", "with_anchor_accuracy"])
    for sample, db in dbs.items():
        for ref in sorted({k[0] for k in db.status}):
            t.add(sample, ref_label(ref), db.count(ref, "done"), db.count(ref, "no_auto"),
                  db.count(ref, "locked"), db.count(ref, "missing"),
                  db.coverage(ref, "block_ng50"), db.coverage(ref, "bases_tagged_fraction"),
                  db.coverage(ref, "read_precision"), db.coverage(ref, "read_precision_anchor_truthtags"))
    return t


def n50(values):
    """The value at which the sorted values (largest first) cover half their sum."""
    s = sorted(values)
    half, acc = sum(s) / 2, 0
    for v in reversed(s):
        acc += v
        if acc >= half:
            return v
    return 0


def read_table_stats(path):
    """Primary alignments on chr1-22 of one read table (plots/read_table.py).
    Read length counts the whole read, clips included; the aligned span (end -
    start on the reference) is what can link het SNPs, and the two part when
    long reads align only in pieces (soft-clipped, the rest as supplementary
    alignments)."""
    lengths = array("l")
    spans = array("l")
    mapq20 = nm = aligned = aligned_all = softclip = 0
    bq = []
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not is_autosome(f[col["chrom"]]):
                continue
            lengths.append(int(f[col["length"]]))
            spans.append(int(f[col["end"]]) - int(f[col["start"]]))
            aligned_all += int(f[col["aligned"]])
            softclip += int(f[col["softclip"]])
            if int(f[col["mapq"]]) >= 20:
                mapq20 += 1
            n, a = int(f[col["nm"]]), int(f[col["aligned"]])
            if n >= 0:
                nm += n
                aligned += a
            b = f[col["mean_bq"]]
            if b not in ("", "nan"):
                bq.append(float(b))
    if not lengths:
        return None
    s = sorted(lengths)
    total = sum(s)
    return {
        "reads": len(s), "gb": total / 1e9, "mean_kb": total / len(s) / 1e3,
        "median_kb": statistics.median(s) / 1e3, "n50_kb": n50(s) / 1e3,
        "ge10kb_pct": 100 * sum(1 for v in s if v >= 10000) / len(s),
        "ge50kb_pct": 100 * sum(1 for v in s if v >= 50000) / len(s),
        "aligned_pct": 100 * aligned_all / total,
        "softclip_pct": 100 * softclip / total,
        "span_n50_kb": n50(spans) / 1e3,
        "span_ge50kb_pct": 100 * sum(1 for v in spans if v >= 50000) / len(spans),
        "mapq20_pct": 100 * mapq20 / len(s),
        "identity_pct": 100 * (1 - nm / aligned) if aligned else NAN,
        "median_bq": statistics.median(bq) if bq else NAN,
    }


def depth_table(path):
    """platform -> value of the last column of source_depth / achieved_depth"""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == len(header):
                out.setdefault(f[0], []).append(f)
    return out


def t01_datasets(root, inputs, samples, threads, with_reads=True):
    t = Table("t01_datasets", "Sequencing data",
              "Per sample and platform. Source depth: mean depth of the provider's BAM on GRCh38 "
              "(mosdepth, data/source_depth.tsv) and of our alignment to CHM13. Read statistics "
              "from the 30x GRCh38 cell, primary alignments on chr1-22 (reads/<cell>.reads.tsv.gz): "
              "read N50 over read length (clips included); aligned: read bases in the "
              "alignment, of all read bases; soft-clipped: likewise; aligned span: end - start of "
              "the primary alignment on the reference, the stretch a read can link het SNPs over; "
              "identity = 1 - NM / aligned bases. Achieved depths: "
              "mosdepth of the downsampled BAMs. Model and aligner as written in the BAM header.",
              ["sample", "platform", "basecaller_model", "aligner_grch38", "source_depth_grch38",
               "source_depth_t2t", "achieved_depths_grch38", "reads_30x", "gb_30x", "mean_kb",
               "median_kb", "read_n50_kb", "reads_ge10kb_pct", "reads_ge50kb_pct", "aligned_pct",
               "softclip_pct", "aligned_span_n50_kb", "span_ge50kb_pct", "mapq_ge20_pct",
               "identity_pct", "median_base_quality"])
    jobs = []
    for sample in samples:
        sroot = sample_root(root, sample)
        for p in st.PLATFORMS:
            path = os.path.join(sroot, "reads", "%s_GRCh38_%s_30x_wgs.reads.tsv.gz" % (p, sample))
            jobs.append((sample, p, path))
    present = [j for j in jobs if os.path.exists(j[2])] if with_reads else []
    stats = {}
    if present:
        with Pool(max(1, min(threads, len(present)))) as pool:
            stats = dict(zip([(s, p) for s, p, _ in present],
                             pool.map(read_table_stats, [path for _, _, path in present])))
    for sample in samples:
        sroot = sample_root(root, sample)
        src38 = depth_table(os.path.join(sroot, "data", "source_depth.tsv"))
        srct2t = depth_table(os.path.join(sroot, "data", "source_depth_T2T.tsv"))
        ach = depth_table(os.path.join(sroot, "data", "achieved_depth.tsv"))
        info = inputs.get((sample, "GRCh38"), {})
        for p in st.PLATFORMS:
            if "bam_" + p not in info and (sample, p) not in stats and p not in src38:
                continue
            s = stats.get((sample, p)) or {}
            achieved = ", ".join("%sx:%s" % (r[1], r[-1]) for r in ach.get(p, []))
            t.add(sample, plat(p), info.get("model_" + p, ""), info.get("aligner_" + p, ""),
                  float(src38[p][0][-1]) if p in src38 else None,
                  float(srct2t[p][0][-1]) if p in srct2t else None, achieved,
                  s.get("reads"), s.get("gb"), s.get("mean_kb"), s.get("median_kb"), s.get("n50_kb"),
                  s.get("ge10kb_pct"), s.get("ge50kb_pct"), s.get("aligned_pct"), s.get("softclip_pct"),
                  s.get("span_n50_kb"), s.get("span_ge50kb_pct"), s.get("mapq20_pct"),
                  s.get("identity_pct"), s.get("median_bq"))
    return t


def truth_stats(job):
    """Variant classes of one truth VCF inside its BED, chr1-22."""
    sample, ref, vcf, bed, phased = job
    region = IntervalSet.from_bed(bed, autosomes_only=True)
    n = {k: 0 for k in ("records", "het_snp_01", "het_snp_12", "het_indel", "hom_snp",
                        "hom_indel", "hemizygous", "other", "het_snp_phased")}
    for chrom, pos, r, alt, gt, ps in vcf_records(vcf, autosomes_only=True):
        if not region.contains(chrom, pos - 1):
            continue
        n["records"] += 1
        alleles = gt.replace("|", "/").split("/")
        alts = alt.split(",")
        if "." in alleles or not all(a.isdigit() for a in alleles):
            n["other"] += 1
            continue
        idx = [int(a) for a in alleles]
        used = [alts[i - 1] for i in idx if 0 < i <= len(alts)]
        snp = len(r) == 1 and all(len(a) == 1 for a in used)
        if len(idx) == 1:
            n["hemizygous"] += 1
        elif idx[0] == idx[1]:
            if idx[0] == 0:
                n["other"] += 1
            else:
                n["hom_snp" if snp else "hom_indel"] += 1
        elif 0 in idx:
            if snp:
                n["het_snp_01"] += 1
                if "|" in gt:
                    n["het_snp_phased"] += 1
            else:
                n["het_indel"] += 1
        else:
            if snp:
                n["het_snp_12"] += 1
            else:
                n["het_indel"] += 1
    bench, gaps, last = 0, [], {}
    if phased and os.path.exists(phased):
        for chrom, pos, r, alt, gt, ps in vcf_records(phased, autosomes_only=True):
            bench += 1
            if chrom in last:
                gaps.append(pos - last[chrom])
            last[chrom] = pos
    # no phased truth yet (make_truth.sh not run): empty, not zero
    n["benchmark"] = bench if phased and os.path.exists(phased) else NAN
    n["bed_mb"] = region.bp() / 1e6
    n["median_gap"] = statistics.median(gaps) if gaps else NAN
    return (sample, ref, os.path.basename(vcf), n)


def t02_truth(inputs, threads):
    t = Table("t02_truth", "Truth sets",
              "Per sample and reference, chr1-22, records of the GIAB truth VCF inside its benchmark "
              "BED by genotype class. het SNP 0/1: biallelic heterozygous SNP; het SNP 1/2: two "
              "different ALT SNP alleles; het indel: any heterozygous record with an indel allele. "
              "Benchmark het SNPs: phased biallelic heterozygous SNPs, the sites every phasing "
              "metric is computed on (bin/make_truth.sh); density per Mb of the BED; median "
              "distance between consecutive ones.",
              ["sample", "reference", "truth_file", "bed_mb", "records_in_bed", "het_snp_01",
               "het_snp_12", "het_indel", "hom_alt_snp", "hom_alt_indel", "hemizygous", "other",
               "het_snp_01_phased_pct", "benchmark_het_snps", "benchmark_per_mb",
               "median_distance_bp"])
    jobs = [(s, r, v["truth_vcf"], v["bed"], v.get("truth_phased", ""))
            for (s, r), v in sorted(inputs.items())
            if v.get("truth_vcf") and v.get("bed") and os.path.exists(v["truth_vcf"])
            and os.path.exists(v["bed"])]
    if not jobs:
        return t
    with Pool(max(1, min(threads, len(jobs)))) as pool:
        results = pool.map(truth_stats, jobs)
    for sample, ref, name, n in results:
        t.add(sample, ref_label(ref), name, n["bed_mb"], n["records"], n["het_snp_01"],
              n["het_snp_12"], n["het_indel"], n["hom_snp"], n["hom_indel"], n["hemizygous"],
              n["other"], pct(n["het_snp_phased"], n["het_snp_01"]), n["benchmark"],
              ratio(n["benchmark"], n["bed_mb"]), n["median_gap"])
    return t


def t03_callers(db):
    t = Table("t03_callers", "Variant callers: het SNP calls against the truth",
              "Per cell, from bin/calls_fdr.sh: biallelic heterozygous PASS SNPs inside the "
              "benchmark BED on chr1-22, positions with more than one record left out. "
              "FNR = FN / truth het SNPs; FDR = FP / calls; F1 = 2 TP / (2 TP + FP + FN). The "
              "calls are the same for every phaser of the cell.",
              ["reference", "caller", "platform", "depth", "truth_het_snps", "calls", "true_positives",
               "false_negatives", "false_positives", "fnr_pct", "fdr_pct", "f1_pct",
               "calls_outside_bed", "multi_record_positions"])
    keys = sorted(db.calls, key=lambda k: (st.REFERENCES.index(k[0]) if k[0] in st.REFERENCES else 9,
                                           st.CALLERS.index(k[1]) if k[1] in st.CALLERS else 9,
                                           st.PLATFORMS.index(k[2]) if k[2] in st.PLATFORMS else 9, k[3]))
    for k in keys:
        c = db.calls[k]
        tp, fp, fn = c.get("shared", NAN), c.get("fp", NAN), c.get("fn", NAN)
        f1 = pct(2 * tp, 2 * tp + fp + fn)
        t.add(ref_label(k[0]), caller_label(k[1]), plat(k[2]), k[3], c.get("truth_het_snps"),
              c.get("calls"), tp, fn, fp, 100 * c.get("fnr", NAN), 100 * c.get("fdr", NAN), f1,
              c.get("calls_outside_bed"), c.get("multi_record_positions"))
    return t


def sorted_keys(db, refs=None):
    order = lambda k: (st.REFERENCES.index(k[0]), st.CALLERS.index(k[1]), st.PLATFORMS.index(k[2]),
                       k[3], st.PHASERS.index(k[4]))
    return sorted((k for k in db.rows if (refs is None or k[0] in refs)
                   and k[1] in st.CALLERS and k[4] in st.PHASERS), key=order)


def t04_grid(db):
    t = Table("t04_grid", "Every run of HG002",
              "One line per reference x input x platform x depth x phaser, chr1-22. truth_phased: "
              "benchmark het SNPs in a phase block of two or more, of all of them. switch_rate: "
              "switch errors / assessed pairs of consecutive phased sites (whatshap compare). "
              "hamming_rate: blockwise Hamming distance / phased benchmark sites. NG50: phase "
              "block span at which the blocks cover half of the autosome length. correct: "
              "benchmark het SNPs called and phased correctly, of all of them. fnr / fdr: of the "
              "run's VCF against the truth. reads_tagged: reads tagged in a block of >= 2 SNPs; "
              "bases_tagged: their bases. read precision / recall: reads on the truth haplotype "
              "of reads tagged by both / by the truth, blocks oriented by their SNPs (snp) or by "
              "an anchor read (anchor).",
              ["reference", "input", "platform", "depth", "phaser", "benchmark_het_snps",
               "phased_in_blocks", "truth_phased_pct", "blocks", "ng50_mb", "switches",
               "assessed_pairs", "switch_rate_pct", "switch_only_per100", "flip_per100",
               "hamming_rate_pct", "correct_pct", "fnr_pct", "fdr_pct", "reads_tagged_pct",
               "bases_tagged_pct", "read_precision_snp_pct", "read_recall_snp_pct",
               "read_precision_anchor_pct", "read_recall_anchor_pct"])
    for k in sorted_keys(db):
        sw, fl = flips(db, k)
        t.add(ref_label(k[0]), caller_label(k[1]), plat(k[2]), k[3], phaser_label(k[4]),
              db.get(k, "truth_total"), db.get(k, "covered"), truth_phased_pct(db, k),
              db.get(k, "blocks"), n50_mb(db, k), db.get(k, "switches"), db.get(k, "assessed"),
              switch_pct(db, k), sw, fl, 100 * db.get(k, "hamming_rate"), recall_pct(db, k),
              100 * db.get(k, "fnr"), 100 * db.get(k, "fdr"), tagged_pct(db, k),
              tagged_bases_pct(db, k), read_precision_pct(db, k), read_recall_pct(db, k),
              anchor_precision_pct(db, k), anchor_recall_pct(db, k))
    return t


def t05_main(db, depth):
    t = Table("t05_main_30x", "Headline numbers: phasers on the GIAB truth input, %dx" % depth,
              "The phaser alone: every phaser gets the same unphased truth het SNPs, so the "
              "numbers compare directly. Definitions as in t04_grid.",
              ["reference", "platform", "phaser", "switch_rate_pct", "ng50_mb", "truth_phased_pct",
               "hamming_rate_pct", "reads_tagged_pct", "bases_tagged_pct", "read_precision_snp_pct",
               "read_recall_snp_pct"])
    for ref in st.REFERENCES:
        for p in st.PLATFORMS:
            for ph in st.PHASERS:
                k = (ref, "truth", p, depth, ph)
                if not db.has(k):
                    continue
                t.add(ref_label(ref), plat(p), phaser_label(ph), switch_pct(db, k), n50_mb(db, k),
                      truth_phased_pct(db, k), 100 * db.get(k, "hamming_rate"), tagged_pct(db, k),
                      tagged_bases_pct(db, k), read_precision_pct(db, k), read_recall_pct(db, k))
    return t


def quartiles(values):
    v = sorted(values)
    if len(v) < 2:
        return (v[0], v[0]) if v else (NAN, NAN)
    q = statistics.quantiles(v, n=4, method="inclusive")
    return q[0], q[2]


def paired_p(xs, ys):
    """Two-sided Wilcoxon signed-rank p of log(T2T / GRCh38) = 0."""
    if wilcoxon is None or len(xs) < 6:
        return NAN
    d = [math.log(y / x) for x, y in zip(xs, ys)]
    if all(v == 0 for v in d):
        return 1.0
    try:
        return float(wilcoxon(d).pvalue)
    except ValueError:
        return NAN


def reference_pairs(db, value, platforms, phasers):
    a, b = st.REFERENCES
    pairs = []
    for p in platforms:
        for ph in phasers:
            for d in st.DEPTHS:
                x, y = value(db, (a, "truth", p, d, ph)), value(db, (b, "truth", p, d, ph))
                if finite(x) and finite(y) and x > 0 and y > 0:
                    pairs.append((x, y))
    return pairs


def t06_reference(db, depth):
    a, b = st.REFERENCES
    t = Table("t06a_reference", "GRCh38 against T2T-CHM13, per platform and phaser",
              "Truth input only (the caller is out of the comparison). Values at %dx on both "
              "references and their ratio T2T / GRCh38; then over the five depths: median ratio, "
              "its range, and in how many depths T2T is better (lower switch rate, higher NG50, "
              "higher share phased). The two references also differ in the alignment: provider "
              "BAMs on GRCh38, our minimap2 on T2T (methodology 10.4)." % depth,
              ["platform", "phaser", "metric", "grch38_%dx" % depth, "t2t_%dx" % depth,
               "ratio_%dx" % depth, "median_ratio", "min_ratio", "max_ratio", "t2t_better", "depths"])
    pooled = Table("t06b_reference_pooled", "GRCh38 against T2T-CHM13, pooled",
                   "Ratio T2T / GRCh38 of the same run (truth input), pooled over platforms and "
                   "depths per phaser, and over phasers and depths per technology: median, "
                   "interquartile range, runs where T2T is better, and the two-sided Wilcoxon "
                   "signed-rank p of log ratio = 0 (empty without scipy or under 6 pairs).",
                   ["group", "metric", "pairs", "median_ratio", "q1_ratio", "q3_ratio",
                    "t2t_better", "wilcoxon_p"])
    for p in st.PLATFORMS:
        for ph in st.PHASERS:
            for metric, value, more in METRICS[:3]:
                pairs = reference_pairs(db, value, [p], [ph])
                if not pairs:
                    continue
                x, y = value(db, (a, "truth", p, depth, ph)), value(db, (b, "truth", p, depth, ph))
                r = [yy / xx for xx, yy in pairs]
                better = sum(1 for xx, yy in pairs if (yy > xx if more else yy < xx))
                t.add(plat(p), phaser_label(ph), metric, x, y, ratio(y, x), statistics.median(r),
                      min(r), max(r), better, len(pairs))
    groups = [(phaser_label(ph), st.PLATFORMS, [ph]) for ph in st.PHASERS]
    groups += [("ONT, all phasers", ["ont4k", "ont5k"], st.PHASERS),
               ("PacBio, all phasers", ["hifi", "revio"], st.PHASERS),
               ("All runs", st.PLATFORMS, st.PHASERS)]
    for label, platforms, phasers in groups:
        for metric, value, more in METRICS[:3]:
            pairs = reference_pairs(db, value, platforms, phasers)
            if not pairs:
                continue
            r = [y / x for x, y in pairs]
            q1, q3 = quartiles(r)
            better = sum(1 for x, y in pairs if (y > x if more else y < x))
            pooled.add(label, metric, len(pairs), statistics.median(r), q1, q3, better,
                       paired_p([x for x, _ in pairs], [y for _, y in pairs]))
    return t, pooled


def t07_ranking(db):
    t = Table("t07_ranking", "How often each phaser is the best",
              "For each reference, input and metric: the phasers are ranked within every "
              "platform x depth cell where all of them are computed (rank 1 = best: lowest switch "
              "rate, highest NG50 / share phased / read recall). wins = cells ranked first; "
              "median_rank over the cells.",
              ["reference", "input", "metric", "phaser", "cells", "wins", "median_rank",
               "median_value"])
    for ref in st.REFERENCES:
        for c in st.CALLERS:
            for metric, value, more in METRICS:
                ranks = {ph: [] for ph in st.PHASERS}
                vals = {ph: [] for ph in st.PHASERS}
                for p in st.PLATFORMS:
                    for d in st.DEPTHS:
                        v = {ph: value(db, (ref, c, p, d, ph)) for ph in st.PHASERS}
                        if not all(finite(x) for x in v.values()):
                            continue
                        order = sorted(st.PHASERS, key=lambda ph: -v[ph] if more else v[ph])
                        for i, ph in enumerate(order):
                            ranks[ph].append(i + 1)
                            vals[ph].append(v[ph])
                for ph in st.PHASERS:
                    if ranks[ph]:
                        t.add(ref_label(ref), caller_label(c), metric, phaser_label(ph), len(ranks[ph]),
                              sum(1 for r in ranks[ph] if r == 1), float(statistics.median(ranks[ph])),
                              statistics.median(vals[ph]))
    return t


def t08_coverage(db):
    t = Table("t08_coverage", "Metrics against coverage",
              "Truth input. Value at each depth, relative change from 10x to 30x and from 30x to "
              "50x (percent of the earlier value), and the saturation depth: the lowest depth "
              "whose value is within 5% of the best value over all depths.",
              ["reference", "platform", "phaser", "metric", "10x", "20x", "30x", "40x", "50x",
               "change_10_30_pct", "change_30_50_pct", "saturation_depth"])
    for ref in st.REFERENCES:
        for p in st.PLATFORMS:
            for ph in st.PHASERS:
                for metric, value, more in METRICS[:3]:
                    v = {d: value(db, (ref, "truth", p, d, ph)) for d in st.DEPTHS}
                    shown = {d: x for d, x in v.items() if finite(x)}
                    if not shown:
                        continue
                    best = max(shown.values()) if more else min(shown.values())
                    sat = next((d for d in st.DEPTHS if d in shown and
                                abs(shown[d] - best) <= 0.05 * abs(best)), None)
                    t.add(ref_label(ref), plat(p), phaser_label(ph), metric,
                          *[v[d] for d in st.DEPTHS],
                          100 * ratio(v[30] - v[10], v[10]), 100 * ratio(v[50] - v[30], v[30]), sat)
    return t


def t09_call_vs_phase(db, depth):
    t = Table("t09_call_vs_phase", "What calling and what phasing cost, %dx" % depth,
              "All on one denominator, the benchmark het SNPs of the reference. called: found by "
              "the caller; correct: found and phased correctly; lost_calling = 100 - called; "
              "lost_phasing = called - correct. The switch rate is over the sites this caller "
              "found and this phaser phased, so it is not a caller-against-caller comparison "
              "on its own; truth_phased says how much of the truth it covers.",
              ["reference", "platform", "input", "phaser", "called_pct", "correct_pct",
               "lost_calling_pct", "lost_phasing_pct", "truth_phased_pct", "switch_rate_pct"])
    for ref in st.REFERENCES:
        for p in st.PLATFORMS:
            for c in st.CALLERS:
                for ph in st.PHASERS:
                    k = (ref, c, p, depth, ph)
                    if not db.has(k):
                        continue
                    called = pct(db.get(k, "shared"), db.get(k, "truth_total"))
                    correct = recall_pct(db, k)
                    t.add(ref_label(ref), plat(p), caller_label(c), phaser_label(ph), called, correct,
                          100 - called, called - correct, truth_phased_pct(db, k), switch_pct(db, k))
    return t


def t10_error_types(db, depth):
    t = Table("t10_error_types", "Switch against flip errors, truth input, %dx" % depth,
              "whatshap compare on chr1-22: a flip is a single site phased the wrong way (two "
              "switches in a row, counted once); switch_only are the remaining switches, a "
              "stretch phased the wrong way. Per 100 assessed pairs of consecutive phased sites.",
              ["reference", "platform", "phaser", "switch_only_per100", "flip_per100",
               "errors_per100", "flip_share_pct"])
    for ref in st.REFERENCES:
        for p in st.PLATFORMS:
            for ph in st.PHASERS:
                k = (ref, "truth", p, depth, ph)
                if not db.has(k):
                    continue
                sw, fl = flips(db, k)
                t.add(ref_label(ref), plat(p), phaser_label(ph), sw, fl, sw + fl, pct(fl, sw + fl))
    return t


def t11_reads(db, depth):
    t = Table("t11_reads_30x", "Read-level metrics, %dx, every input" % depth,
              "Reads and bases tagged in blocks of >= 2 phased SNPs on chr1-22; precision and "
              "recall of the tags against the truth tags, blocks oriented by their SNPs or by an "
              "anchor read (definitions in t04_grid).",
              ["reference", "platform", "input", "phaser", "reads_tagged_pct", "bases_tagged_pct",
               "read_precision_snp_pct", "read_recall_snp_pct", "read_precision_anchor_pct",
               "read_recall_anchor_pct"])
    for ref in st.REFERENCES:
        for p in st.PLATFORMS:
            for c in st.CALLERS:
                for ph in st.PHASERS:
                    k = (ref, c, p, depth, ph)
                    if not db.has(k):
                        continue
                    t.add(ref_label(ref), plat(p), caller_label(c), phaser_label(ph), tagged_pct(db, k),
                          tagged_bases_pct(db, k), read_precision_pct(db, k), read_recall_pct(db, k),
                          anchor_precision_pct(db, k), anchor_recall_pct(db, k))
    return t


def stratum_name(s):
    from figures import stratum_label
    return stratum_label(s).replace("\n", " ")


def t12_strata(db, depth):
    a = Table("t12a_strata", "The strata",
              "GIAB genome stratifications v3.6, chr1-22, each cut to the benchmark BED "
              "(regions_statistics.sh). bp_raw: the whole stratum; bp: inside the BED.",
              ["reference", "stratum", "bp_raw_mb", "bp_in_bed_mb", "genome_pct",
               "benchmark_het_snps", "het_snps_per_mb", "median_distance_bp", "within_50bp_pct"])
    b = Table("t12b_strata_switch", "Switch rate per stratum, truth input, %dx" % depth,
              "Switch errors / assessed pairs inside each stratum (regions.sh on chr1-22). "
              "Enrichment: switch errors (block ends) per het SNP in the stratum relative to "
              "the genome (bin/strata_mechanism.sh), where computed.",
              ["reference", "stratum", "platform", "phaser", "benchmark_het_snps", "switches",
               "assessed_pairs", "switch_rate_pct", "switch_enrichment", "block_end_enrichment"])
    for ref in st.REFERENCES:
        stats = db.strata_statistics(ref)
        names = [s for s in db.strata_names() if s in stats]
        for s in names:
            r = stats[s]
            a.add(ref_label(ref), stratum_name(s), r.get("bp_raw", NAN) / 1e6, r.get("bp", NAN) / 1e6,
                  100 * r.get("frac_genome", NAN), r.get("truth_het_snps"), r.get("het_snps_per_mb"),
                  r.get("median_inter_snp_dist"), 100 * r.get("frac_snp_within_50bp", NAN))
        mech = db.mechanism(ref)
        for s in db.strata_names():
            for p in st.PLATFORMS:
                for ph in st.PHASERS:
                    row = db.stratum((ref, "truth", p, depth, ph), s)
                    if not row:
                        continue
                    b.add(ref_label(ref), stratum_name(s), plat(p), phaser_label(ph),
                          row["truth_total"], row["switches"], row["assessed"],
                          pct(row["switches"], row["assessed"]),
                          mech.get((p, depth, "truth", ph, s, "switch_enrichment"), NAN),
                          mech.get((p, depth, "truth", ph, s, "break_enrichment"), NAN))
    return (a, b, t12c_strata_inputs(db, depth), t12d_read_classes(db, depth),
            t12e_strata_reads(db, depth), t12f_strata_features(db, depth))


def t12c_strata_inputs(db, depth):
    t = Table("t12c_strata_inputs", "Phasing per stratum, every input, %dx" % depth,
              "regions.sh: the run's VCF and the truth cut to the stratum (chr1-22). Phased: "
              "phased / all het SNPs of the input in the stratum; switch and Hamming rates over "
              "assessed pairs and covered SNPs; correct: truth het SNPs called and phased right, "
              "of all truth het SNPs in the stratum; FNR / FDR of the input's het SNPs there.",
              ["reference", "stratum", "platform", "input", "phaser", "benchmark_het_snps",
               "phased_pct", "switch_rate_pct", "hamming_rate_pct", "correct_pct", "fnr_pct",
               "fdr_pct"])
    for ref in st.REFERENCES:
        for s in db.strata_names():
            for p in st.PLATFORMS:
                for c in st.CALLERS:
                    for ph in st.PHASERS:
                        row = db.stratum((ref, c, p, depth, ph), s)
                        if not row:
                            continue
                        t.add(ref_label(ref), stratum_name(s), plat(p), caller_label(c),
                              phaser_label(ph), row["truth_total"], pct(row["phased"], row["variants"]),
                              pct(row["switches"], row["assessed"]), pct(row["hamming"], row["covered"]),
                              pct(row["correct"], row["truth_total"]), 100 * row["fnr"],
                              100 * row["fdr"])
    return t


def t12d_read_classes(db, depth):
    from figures_extra import CLASS_ORDER
    t = Table("t12d_read_classes", "Why reads are not tagged, truth input, %dx" % depth,
              "bin/read_classes.sh: every primary alignment on chr1-22 in exactly one class, the "
              "first rule that fits (tagged; MAPQ < 20 in segdups / low mappability; MAPQ < 20 "
              "elsewhere; no het SNP and outside the confident regions; no het SNP; SNPs only in "
              "no block of >= 2 SNPs; SNPs in a block yet no haplotype). Percent of all reads and "
              "of all read bases of the run.",
              ["reference", "platform", "phaser", "reads"]
              + ["reads_%s_pct" % c for c in CLASS_ORDER] + ["bases_%s_pct" % c for c in CLASS_ORDER])
    for ref in st.REFERENCES:
        for p in st.PLATFORMS:
            for ph in st.PHASERS:
                table = db.classes((ref, "truth", p, depth, ph))
                rows = [table.get(("genome", "class", c), {}) for c in CLASS_ORDER]
                reads = sum(r.get("reads", 0) or 0 for r in rows)
                bases = sum(r.get("bases", 0) or 0 for r in rows)
                if not reads:
                    continue
                t.add(ref_label(ref), plat(p), phaser_label(ph), reads,
                      *[pct(r.get("reads", 0) or 0, reads) for r in rows],
                      *[pct(r.get("bases", 0) or 0, bases) for r in rows])
    return t


def t12e_strata_reads(db, depth):
    t = Table("t12e_strata_reads", "Reads per stratum, truth input, %dx" % depth,
              "bin/read_classes.sh: reads with half their alignment or more in the stratum. "
              "Tagged: of those reads; wrong haplotype: tagged reads with a truth label whose "
              "haplotype disagrees with it, blocks oriented by their SNPs; mismatches and indels "
              "per aligned base of all the reads in the stratum.",
              ["reference", "stratum", "platform", "phaser", "reads", "reads_tagged_pct",
               "wrong_haplotype_pct", "mismatch_pct", "indel_pct"])
    for ref in st.REFERENCES:
        for s in ["genome"] + db.strata_names():
            for p in st.PLATFORMS:
                for ph in st.PHASERS:
                    table = db.classes((ref, "truth", p, depth, ph))
                    allr = table.get((s, "all", "all"))
                    if not allr or not allr.get("reads"):
                        continue
                    tagged = table.get((s, "class", "tagged"), {})
                    wrong = NAN
                    if tagged.get("compared"):
                        wrong = 100 - pct(tagged["correct"], tagged["compared"])
                    t.add(ref_label(ref), "Genome" if s == "genome" else stratum_name(s), plat(p),
                          phaser_label(ph), allr["reads"], pct(tagged.get("reads", 0) or 0, allr["reads"]),
                          wrong, 100 * allr.get("mismatch_rate", NAN), 100 * allr.get("indel_rate", NAN))
    return t


def t12f_strata_features(db, depth):
    t = Table("t12f_strata_features", "What makes a stratum hard, %dx" % depth,
              "bin/strata_mechanism.sh, the values of fig11c: |VAF - 0.5| of the truth het SNPs "
              "in the reads; links between two SNPs showing more than two allele patterns (a "
              "paralog sign); mean depth relative to the genome; share of the depth at MAPQ < 20.",
              ["reference", "stratum", "platform", "truth_het_snps", "vaf_dev",
               "paralog_links_pct", "depth_ratio", "lowmapq_depth_pct"])
    for ref in st.REFERENCES:
        mech = db.mechanism(ref)
        for s in db.strata_names():
            for p in st.PLATFORMS:
                get = lambda m: mech.get((p, depth, "-", "-", s, m), NAN)
                if not finite(get("truth_snps")):
                    continue
                t.add(ref_label(ref), stratum_name(s), plat(p), get("truth_snps"), get("vaf_dev"),
                      100 * get("paralog_link_frac"), get("depth_ratio"), 100 * get("lowmapq_frac"))
    return t


def t13_samples(dbs, depth=30, ref="GRCh38"):
    t = Table("t13_samples", "GIAB samples at %dx on GRCh38" % depth,
              "Truth input for the phasing metrics (v5.0q truth for HG002, v4.2.1 for the others, "
              "chr1-22), and the FNR / FDR of each caller in the same cell.",
              ["sample", "platform", "phaser", "switch_rate_pct", "ng50_mb", "truth_phased_pct",
               "read_recall_snp_pct", "clair3_fnr_pct", "deepvariant_fnr_pct", "longcalld_fnr_pct",
               "clair3_fdr_pct", "deepvariant_fdr_pct", "longcalld_fdr_pct"])
    for sample, db in dbs.items():
        for p in st.PLATFORMS:
            calls = [100 * db.call(ref, c, p, depth, "fnr") for c in ("clair3", "deepvariant", "longcalld")]
            calls += [100 * db.call(ref, c, p, depth, "fdr") for c in ("clair3", "deepvariant", "longcalld")]
            for ph in st.PHASERS:
                k = (ref, "truth", p, depth, ph)
                if not db.has(k):
                    continue
                t.add(sample, plat(p), phaser_label(ph), switch_pct(db, k), n50_mb(db, k),
                      truth_phased_pct(db, k), read_recall_pct(db, k), *calls)
    return t


def t14_software(path):
    t = Table("t14_software", "Software",
              "Versions as installed on the servers (conda-meta, file names, --version); "
              "'?' is to be filled in by hand.", ["tool", "version", "where"])
    if path and os.path.exists(path):
        with open(path) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == 3:
                    t.add(*f)
    return t


# --------------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--out", default=None, help="default <root>/tables")
    parser.add_argument("--inputs", default=None, help="from bin/make_tables.sh")
    parser.add_argument("--software", default=None, help="from bin/make_tables.sh")
    parser.add_argument("--depth", type=int, default=30, choices=st.DEPTHS)
    parser.add_argument("--phasers", default=" ".join(st.PHASERS))
    parser.add_argument("--rules", default="auto_min2")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-reads", action="store_true", help="no read statistics in t01")
    parser.add_argument("--skip-truth", action="store_true", help="no t02 (reads the truth VCFs)")
    ctx = parser.parse_args()
    ctx.out = ctx.out or os.path.join(ctx.root, "tables")
    os.makedirs(ctx.out, exist_ok=True)
    st.PHASERS[:] = [p for p in st.PHASERS if p in ctx.phasers.split()]
    if not st.PHASERS:
        parser.error("--phasers names none of %s" % " ".join(st.PHASER_LABEL))

    inputs = read_inputs(ctx.inputs)
    samples = ["HG002"] + sorted(s for s in {s for s, _ in inputs} if s != "HG002"
                                 and os.path.isdir(sample_root(ctx.root, s)))
    db = Metrics(ctx.root, "HG002", st.REGION, readacc_rules=ctx.rules).load(
        st.REFERENCES, st.CALLERS, st.PLATFORMS, st.DEPTHS, st.PHASERS)
    dbs = {"HG002": db}
    for s in samples[1:]:
        dbs[s] = Metrics(ctx.root, s, st.REGION, readacc_rules=ctx.rules).load(
            ["GRCh38"], st.CALLERS, st.PLATFORMS, [30], st.PHASERS)

    tables = [t00_completeness(dbs)]
    tables.append(t01_datasets(ctx.root, inputs, samples, ctx.threads, with_reads=not ctx.skip_reads))
    if not ctx.skip_truth:
        tables.append(t02_truth(inputs, ctx.threads))
    tables += [t03_callers(db), t04_grid(db), t05_main(db, ctx.depth)]
    tables += list(t06_reference(db, ctx.depth))
    tables += [t07_ranking(db), t08_coverage(db), t09_call_vs_phase(db, ctx.depth),
               t10_error_types(db, ctx.depth), t11_reads(db, ctx.depth)]
    tables += list(t12_strata(db, ctx.depth))
    tables += [t13_samples({s: d for s, d in dbs.items()}), t14_software(ctx.software)]

    written = [t for t in tables if t.rows]
    for t in tables:
        print("  %-8s %-24s %d rows" % ("ok" if t.rows else "empty", t.name, len(t.rows)))
    readme = Table("README", "What each table holds", "", ["table", "title", "definition"])
    for t in written:
        readme.add(t.name, t.title, t.note)
    for t in [readme] + written:
        write_tsv(t, ctx.out)
    write_xlsx(os.path.join(ctx.out, "tables.xlsx"),
               [(t.name, t.header, t.rows) for t in [readme] + written])
    print("%d tables -> %s (TSV each, and tables.xlsx)" % (len(written), ctx.out))


if __name__ == "__main__":
    main()
