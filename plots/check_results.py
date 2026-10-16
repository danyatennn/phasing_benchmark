#!/usr/bin/env python3
"""Checks the numbers behind the figures before they go into the paper. Reads
only; nothing in metrics/ or phased/ is changed.

  python plots/check_results.py                  every check
  python plots/check_results.py --only grid,rows
  python plots/check_results.py --samtools /path/to/samtools   for the read check

  grid     why a run of the grid does not count as computed
  rows     the identities every row.tsv must satisfy
  ng50     block NG50 recomputed from the phased VCF against a fixed genome,
           the 22 autosomes, next to what whatshap stats reported
  tagging  which reads the haplotag list counts as untagged: by chromosome,
           and on chr20 by read length and MAPQ
  margin   Margin's share of the truth phased against coverage, and what its
           logs say about the sites it dropped
"""

import argparse
import collections
import gzip
import math
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import style as st
from data import ROW_COLUMNS, Metrics

AUTOSOMES = ["chr%d" % i for i in range(1, 23)]


def title(text):
    print("\n" + "=" * 78 + "\n" + text + "\n" + "=" * 78)


def pct(a, b):
    return 100.0 * a / b if b else float("nan")


# ------------------------------------------------------------------- grid

def check_grid(db, ctx):
    title("grid: runs the figures do not use")
    bad = 0
    for key, state in sorted(db.status.items()):
        if state == "done" or key[0] not in ctx.refs:
            continue
        bad += 1
        path = os.path.join(db.metrics_dir, db.name(*key) + ".row.tsv")
        if state == "locked":
            why = "locked: " + (open(path + ".lock/owner").read().strip()
                                if os.path.exists(path + ".lock/owner") else "no owner file")
        elif not os.path.exists(path):
            why = "no row.tsv"
        else:
            with open(path) as fh:
                lines = fh.read().split("\n")
            lines = [l for l in lines if l]
            first = lines[0].split("\t") if lines else []
            if len(first) != len(ROW_COLUMNS):
                why = "%d fields instead of %d, %d line(s)" % (len(first), len(ROW_COLUMNS), len(lines))
            else:
                said = dict(zip(ROW_COLUMNS, first))
                why = "row says %s %s %s %s %s %s" % tuple(
                    said[c] for c in ("reference", "region", "caller", "platform", "depth", "phaser"))
            why += "\n      " + (lines[0][:200] if lines else "(empty)")
        print("  %-55s %s" % (db.name(*key), why))
    print("  all runs are computed" if not bad else "  %d run(s) left out" % bad)


# ------------------------------------------------------------------- rows

def check_rows(db, ctx):
    """Each line is an identity that holds by construction; a failure is a
    run scored on a truncated file or a column read from the wrong place."""
    title("rows: identities that must hold")
    truth_totals = collections.defaultdict(set)
    failures = 0
    for key, row in sorted(db.rows.items()):
        if key[0] not in ctx.refs:
            continue
        g = row.get
        truth_totals[key[0]].add(g("truth_total"))
        tests = [
            ("hamming <= covered", g("hamming") <= g("covered")),
            ("switches <= assessed pairs", g("switches") <= g("assessed")),
            ("covered <= shared", g("covered") <= g("shared")),
            ("shared <= emitted", g("shared") <= g("emitted")),
            ("shared <= truth_total", g("shared") <= g("truth_total")),
            ("phased <= variants", g("phased") <= g("variants")),
            ("reads_tagged <= reads_total", g("reads_tagged") <= g("reads_total")),
            ("block_n50 > 0", g("block_n50") > 0),
        ]
        if key[1] == "truth":
            # the phaser got the truth sites, so it can only lose them
            tests.append(("truth input: emitted <= truth_total", g("emitted") <= g("truth_total")))
        for label, ok in tests:
            if not ok:
                failures += 1
                print("  FAIL %-45s %s" % (label, db.name(*key)))
    for ref, totals in sorted(truth_totals.items()):
        print("  %-7s truth_total %s" % (ref, ", ".join("%d" % t for t in sorted(totals))))
        if len(totals) != 1:
            failures += 1
            print("  FAIL truth_total differs between runs of %s" % ref)

    # a truth-input run that drops sites: the phaser re-genotyped them, which
    # changes what its switch rate is measured on
    print("\n  truth input, more than 0.1% of the sites not written back by the phaser:")
    for key, row in sorted(db.rows.items()):
        if key[0] in ctx.refs and key[1] == "truth" and row["emitted"] < 0.999 * row["truth_total"]:
            print("    %-50s %.3f%% dropped" % (
                db.name(*key), pct(row["truth_total"] - row["emitted"], row["truth_total"])))
    print("  %d failure(s)" % failures)


# ------------------------------------------------------------------- NG50

def chrom_lengths(root, ref):
    lengths = {}
    with open(os.path.join(root, "data", "chrom_lengths_%s.tsv" % ref)) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2:
                lengths[f[0]] = int(f[1])
    return lengths


def phase_blocks(vcf):
    """(chrom, PS) -> [first, last, n] over phased genotypes, as whatshap stats
    reads them: GT with '|' and a PS."""
    blocks = {}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t", 10)
            keys = f[8].split(":")
            vals = f[9].split(":")
            if "PS" not in keys or "GT" not in keys:
                continue
            gt = vals[keys.index("GT")]
            ps = vals[keys.index("PS")] if keys.index("PS") < len(vals) else "."
            if "|" not in gt or ps in (".", ""):
                continue
            a = gt.split("|")
            if a[0] == a[1]:
                continue
            pos = int(f[1])
            b = blocks.get((f[0], ps))
            if b is None:
                blocks[(f[0], ps)] = [pos, pos, 1]
            else:
                b[0] = min(b[0], pos)
                b[1] = max(b[1], pos)
                b[2] += 1
    return blocks


def ng50(spans, genome):
    acc = 0
    for s in sorted(spans, reverse=True):
        acc += s
        if acc >= genome / 2:
            return s
    return 0


def interleaved(blocks):
    """Blocks whose span overlaps another block's span on the same chromosome:
    their span counts sequence they do not phase."""
    by_chrom = collections.defaultdict(list)
    for (chrom, _), (a, b, n) in blocks.items():
        if n >= 2:
            by_chrom[chrom].append((a, b))
    hits = 0
    for spans in by_chrom.values():
        spans.sort()
        reach = -1
        for a, b in spans:
            if a <= reach:
                hits += 1
            reach = max(reach, b)
    return hits


def check_ng50(db, ctx):
    title("ng50: recomputed from the phased VCFs, %d×, GIAB truth input" % ctx.depth)
    print("  whatshap = block_n50 in row.tsv (whatshap stats --chr-lengths, its own genome size)")
    print("  auto     = NG50 of autosomal blocks of >= 2 SNPs, genome = the 22 autosomes")
    print("  span = last - first phased SNP of a block, as whatshap stats measures it\n")
    for ref in ctx.refs:
        lengths = chrom_lengths(db.root, ref)
        genome = sum(lengths[c] for c in AUTOSOMES if c in lengths)
        print("  %s, autosomes %.1f Mb" % (ref, genome / 1e6))
        print("  %-8s %-10s %10s %10s %8s %9s %12s %s" % (
            "platform", "phaser", "whatshap", "auto", "diff", "blocks", "interleaved",
            "non-autosomal blocks"))
        for platform in st.PLATFORMS:
            for phaser in st.PHASERS:
                key = (ref, "truth", platform, ctx.depth, phaser)
                vcf = os.path.join(db.root, "phased", db.name(*key) + ".vcf.gz")
                if not os.path.exists(vcf):
                    print("  %-8s %-10s no phased VCF" % (platform, phaser))
                    continue
                blocks = phase_blocks(vcf)
                auto = {k: v for k, v in blocks.items() if k[0] in AUTOSOMES and v[2] >= 2}
                other = collections.Counter(k[0] for k, v in blocks.items()
                                            if k[0] not in AUTOSOMES and v[2] >= 2)
                mine = ng50([b - a for a, b, _ in auto.values()], genome) / 1e6
                theirs = db.get(key, "block_n50") / 1e6
                print("  %-8s %-10s %10.3f %10.3f %+7.1f%% %9d %12d %s" % (
                    platform, phaser, theirs, mine, pct(theirs - mine, mine), len(auto),
                    interleaved(auto),
                    " ".join("%s:%d" % kv for kv in sorted(other.items())) or "-"))
        print()


# ---------------------------------------------------------------- tagging

def chrom_class(chrom):
    if chrom in AUTOSOMES:
        return "autosome"
    if chrom in ("chrX", "chrY", "chrM"):
        return chrom
    if chrom in ("", "*", "None", "none"):
        return "unmapped"
    return "other contigs"


def read_tag_list(path):
    """readname -> (tagged?, chromosome of the first line), and how many lines
    the list has: a read can appear once per alignment."""
    reads = {}
    lines = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            lines += 1
            tagged = len(f) > 1 and f[1] != "none"
            chrom = f[3] if len(f) > 3 else ""
            old = reads.get(f[0])
            reads[f[0]] = (tagged or (old[0] if old else False), old[1] if old else chrom)
    return reads, lines


LEN_BINS = [(0, 1000), (1000, 5000), (5000, 10000), (10000, 20000), (20000, 50000), (50000, 10 ** 9)]
MAPQ_BINS = [(0, 1), (1, 20), (20, 60), (60, 256)]


def bin_label(lo, hi, unit=""):
    if hi >= 10 ** 9:
        return ">=%s%s" % (lo // 1000 if unit == "kb" else lo, unit)
    if unit == "kb":
        return "%g-%g kb" % (lo / 1000, hi / 1000)
    return "%d-%d" % (lo, hi - 1)


def check_tagging(db, ctx):
    title("tagging: what 'reads assigned a haplotype' counts, %d×, truth input, WhatsHap" % ctx.depth)
    print("  reads_total = lines of the haplotag list; the figure plots reads_tagged / reads_total\n")
    for ref in ctx.refs:
        print("  " + ref)
        for platform in st.PLATFORMS:
            key = (ref, "truth", platform, ctx.depth, "whatshap")
            path = os.path.join(db.metrics_dir, db.name(*key) + ".haplotag.tsv")
            if not os.path.exists(path):
                print("    %-6s no haplotag list" % platform)
                continue
            reads, lines = read_tag_list(path)
            by = collections.defaultdict(lambda: [0, 0])
            for tagged, chrom in reads.values():
                c = by[chrom_class(chrom)]
                c[0] += 1
                c[1] += tagged
            n = len(reads)
            t = sum(v[1] for v in by.values())
            auto = by["autosome"]
            print("    %-6s lines %d, distinct reads %d; tagged %.2f%% of lines in row.tsv, "
                  "%.2f%% of reads, %.2f%% of autosomal reads" % (
                      platform, lines, n, 100 * db.get(key, "tagging_completeness"),
                      pct(t, n), pct(auto[1], auto[0])))
            for cls in ("autosome", "chrX", "chrY", "chrM", "other contigs", "unmapped"):
                if cls in by:
                    c = by[cls]
                    print("           %-14s %5.2f%% of reads, %6.2f%% of them tagged" % (
                        cls, pct(c[0], n), pct(c[1], c[0])))
            if ctx.samtools:
                by_length(ctx, db, ref, platform, reads)
        print()


def by_length(ctx, db, ref, platform, reads):
    """chr20 primary alignments, joined to the list by read name: untagged by
    read length, by MAPQ, and the share of bases (not reads) that got a tag."""
    bam = os.path.join(db.root, "data", "%s_%s_%s_%dx_wgs" % (platform, ref, st.SAMPLE, ctx.depth),
                       "input.bam")
    if not os.path.exists(bam):
        print("           no BAM %s" % bam)
        return
    lens = collections.defaultdict(lambda: [0, 0])
    mapq = collections.defaultdict(lambda: [0, 0])
    bases = [0, 0]
    proc = subprocess.Popen([ctx.samtools, "view", "-F", "0x904", bam, "chr20"],
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        f = line.split("\t", 11)
        name, q, seq = f[0], int(f[4]), f[9]
        length = len(seq) if seq != "*" else 0
        tagged = reads.get(name, (False, ""))[0]
        for lo, hi in LEN_BINS:
            if lo <= length < hi:
                lens[(lo, hi)][0] += 1
                lens[(lo, hi)][1] += tagged
        for lo, hi in MAPQ_BINS:
            if lo <= q < hi:
                mapq[(lo, hi)][0] += 1
                mapq[(lo, hi)][1] += tagged
        bases[0] += length
        bases[1] += length * tagged
    proc.wait()
    total = sum(v[0] for v in lens.values())
    print("           chr20, %d primary reads; bases in tagged reads %.2f%%, reads tagged %.2f%%" % (
        total, pct(bases[1], bases[0]), pct(sum(v[1] for v in lens.values()), total)))
    print("           %-12s %9s %9s" % ("read length", "of reads", "tagged"))
    for b in LEN_BINS:
        c = lens[b]
        print("           %-12s %8.2f%% %8.2f%%" % (bin_label(*b, unit="kb"), pct(c[0], total), pct(c[1], c[0])))
    print("           %-12s %9s %9s" % ("MAPQ", "of reads", "tagged"))
    for b in MAPQ_BINS:
        c = mapq[b]
        print("           %-12s %8.2f%% %8.2f%%" % (bin_label(*b), pct(c[0], total), pct(c[1], c[0])))


# ----------------------------------------------------------------- margin

MARGIN_LOG = re.compile(r"homozyg|not analy|unanaly|filter|skipp|discard|remov", re.I)


def check_margin(db, ctx):
    title("margin: share of the truth phased against coverage, truth input")
    for ref in ctx.refs:
        print("  %s   truth het SNPs phased (covered / truth_total, %%)" % ref)
        print("  %-8s %-10s %s" % ("platform", "phaser", "  ".join("%9s" % ("%d×" % d) for d in st.DEPTHS)))
        for platform in st.PLATFORMS:
            for phaser in st.PHASERS:
                cells = []
                for d in st.DEPTHS:
                    key = (ref, "truth", platform, d, phaser)
                    cov = pct(db.get(key, "covered"), db.get(key, "truth_total"))
                    cells.append("%9.2f" % cov if math.isfinite(cov) else "%9s" % "-")
                print("  %-8s %-10s %s" % (platform, phaser, "  ".join(cells)))
        print()
        for platform in ("hifi", "revio"):
            for d in (30, 50):
                name = db.name(ref, "truth", platform, d, "margin")
                logf = os.path.join(db.root, "logs", "phase", name + ".log")
                if not os.path.exists(logf):
                    continue
                print("  %s" % logf)
                shown = 0
                with open(logf, errors="replace") as fh:
                    for line in fh:
                        if MARGIN_LOG.search(line) and shown < 12:
                            print("    " + line.rstrip()[:160])
                            shown += 1
        calibration = os.path.join(db.metrics_dir, "calibration.tsv")
        if os.path.exists(calibration):
            print("\n  %s (bin/calibrate_depth.sh, chr20):" % calibration)
            with open(calibration) as fh:
                for line in fh:
                    print("    " + line.rstrip())
        else:
            print("\n  no metrics/calibration.tsv: bin/calibrate_depth.sh has not run")


CHECKS = [("grid", check_grid), ("rows", check_rows), ("ng50", check_ng50),
          ("tagging", check_tagging), ("margin", check_margin)]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=st.ROOT)
    parser.add_argument("--only", default="", help="comma-separated: " + ",".join(n for n, _ in CHECKS))
    parser.add_argument("--depth", type=int, default=30, choices=st.DEPTHS)
    parser.add_argument("--refs", default=" ".join(st.REFERENCES))
    parser.add_argument("--samtools", default="", help="enables the chr20 read-length check")
    ctx = parser.parse_args()
    ctx.refs = [r for r in ctx.refs.split() if r in st.REFERENCES]

    db = Metrics(ctx.root, st.SAMPLE, st.REGION).load(
        st.REFERENCES, st.CALLERS, st.PLATFORMS, st.DEPTHS, st.PHASERS)
    wanted = [w for w in ctx.only.split(",") if w]
    for name, check in CHECKS:
        if not wanted or name in wanted:
            check(db, ctx)


if __name__ == "__main__":
    main()
