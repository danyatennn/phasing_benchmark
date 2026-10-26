#!/usr/bin/env python3
"""Reads that one reference tags and the other does not, for one platform and
coverage: tagged on T2T-CHM13 but not on GRCh38 ("rescued by CHM13") and the
other way round ("lost on CHM13"), per phaser. The numbers behind figR9 of
plots/figures_summary.py.

Reads: those in the read tables of both cells (the T2T cell is drawn from
reads realigned to CHM13, with its own downsampling fraction, so the two read
sets overlap but are not the same; only the shared reads are compared). A read
is tagged by a run if whatshap haplotag gave it H1/H2 in a block of
--min-snps phased het SNPs or more, as everywhere in the paper.

Why a read is untagged on a reference (in this order, as plots/read_classes.py):
  off_autosomes        its primary alignment is not on chr1-22
  lowmapq_paralog      MAPQ < 20, half the alignment or more in segmental
                       duplications or low-mappability regions
  lowmapq_other        MAPQ < 20 elsewhere
  no_snp_outside_bed   no truth het SNP on it, wholly outside the benchmark
  no_snp               no truth het SNP on it
  snps_not_tagged      het SNPs on it, but no H1/H2 in a block that counts
Is the tag right: the read's truth label (whatshap haplotag against the
phased truth) against its label once the block is oriented against the truth
by its SNPs.

Output, with --out <prefix>:
  <prefix>.reads.tsv.gz    every read rescued or lost by a phaser: positions on
                           both references, MAPQ, SNPs, the reason, tier-1
                           strata (half the alignment or more), tag correct
  <prefix>.upset.tsv       reads untagged on GRCh38 by at least one phaser,
                           counted by membership: untagged_GRCh38_<phaser>,
                           tagged_T2T_<phaser> (0/1) -> reads
  <prefix>.regions_T2T.tsv, <prefix>.regions_GRCh38.tsv
                           rescued reads merged into regions (overlapping or
                           closer than --merge-gap): chrom start end reads,
                           reads per phaser, strata the region touches
  <prefix>.summary.tsv     long format: phaser direction category name reads
                           share background_share (share of all shared reads,
                           for the strata)
"""

import argparse
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import IntervalSet, haplotag_list, is_autosome, open_text, phased_het, vcf_records
from read_classes import blocks_and_orientation

MIN_MAPQ = 20
STRATUM_FRACTION = 0.5
REFS = ("GRCh38", "T2T")
REASONS = ["off_autosomes", "lowmapq_paralog", "lowmapq_other", "no_snp_outside_bed", "no_snp",
           "snps_not_tagged"]


def load_truth(path):
    return {(c, p, ref, alt): gt for c, p, ref, alt, gt, ps in vcf_records(path, autosomes_only=True)
            if phased_het(gt) is not None}


def read_rows(path):
    """(read, chrom, start, end, mapq) for every line of a read table"""
    with open_text(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ic, i_s, ie, im = h.index("chrom"), h.index("start"), h.index("end"), h.index("mapq")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            yield f[0], f[ic], int(f[i_s]), int(f[ie]), int(f[im])


def snps_of(path, wanted):
    """read -> het SNPs covered with REF or ALT, for the wanted reads"""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open_text(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ir, ia = h.index("n_ref"), h.index("n_alt")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in wanted:
                out[f[0]] = int(f[ir]) + int(f[ia])
    return out


def truth_of(path, wanted):
    """read -> 0/1 truth haplotype, for the wanted reads"""
    out = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] in wanted and f[0] not in out and f[1] in ("H1", "H2"):
                out[f[0]] = 0 if f[1] == "H1" else 1
    return out


def strata_of(sets, chrom, start, end):
    span = max(1, end - start)
    return [n for n, s in sets if s.overlap(chrom, start, end) >= STRATUM_FRACTION * span]


def reason(pos, nsnp, paralog, bench):
    if pos is None or not is_autosome(pos[0]):
        return "off_autosomes"
    chrom, start, end, mapq = pos
    if mapq < MIN_MAPQ:
        span = max(1, end - start)
        return "lowmapq_paralog" if paralog.overlap(chrom, start, end) >= STRATUM_FRACTION * span \
            else "lowmapq_other"
    if not nsnp:
        return "no_snp" if bench.overlap(chrom, start, end) > 0 else "no_snp_outside_bed"
    return "snps_not_tagged"


def merge(intervals, gap):
    """[(chrom, start, end, phasers)] -> merged [(chrom, start, end, reads, {phaser: reads})]"""
    out = []
    for c, s, e, ph in sorted(intervals):
        if out and out[-1][0] == c and s <= out[-1][2] + gap:
            last = out[-1]
            last[2] = max(last[2], e)
            last[3] += 1
            for p in ph:
                last[4][p] = last[4].get(p, 0) + 1
        else:
            out.append([c, s, e, 1, {p: 1 for p in ph}])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", nargs=4, action="append", required=True,
                    metavar=("REF", "PHASER", "TAGS", "VCF"), help="GRCh38 or T2T, phaser, haplotag list, phased VCF")
    for ref in REFS:
        k = ref.lower()
        ap.add_argument("--reads-%s" % k, required=True, help="read table of the %s cell" % ref)
        ap.add_argument("--alleles-%s" % k, default=None, help="<cell>_truth.reads.tsv.gz on %s" % ref)
        ap.add_argument("--truth-vcf-%s" % k, required=True)
        ap.add_argument("--truth-tags-%s" % k, required=True)
        ap.add_argument("--strata-%s" % k, default=None, help="directory of <ref>_<stratum>.bed")
    ap.add_argument("--stratum", action="append", default=[], help="tier-1 stratum names")
    ap.add_argument("--min-snps", type=int, default=2)
    ap.add_argument("--merge-gap", type=int, default=10000)
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()

    phasers = []
    for ref, phaser, _, _ in args.run:
        if phaser not in phasers:
            phasers.append(phaser)
    runs = {(ref, phaser): (tags, vcf) for ref, phaser, tags, vcf in args.run}
    for ref in REFS:
        for p in phasers:
            if (ref, p) not in runs:
                sys.exit("no %s run for %s" % (ref, p))

    # bits: present on GRCh38, present on T2T, then tagged per (ref, phaser)
    PRES = {"GRCh38": 1, "T2T": 2}
    TAG = {}
    bit = 4
    for ref in REFS:
        for p in phasers:
            TAG[(ref, p)] = bit
            bit <<= 1
    mask = {}
    flips, labels_wanted = {}, set()
    for ref in REFS:
        truth_gt = load_truth(getattr(args, "truth_vcf_%s" % ref.lower()))
        for p in phasers:
            tags, vcf = runs[(ref, p)]
            size, flip = blocks_and_orientation(vcf, truth_gt)
            flips[(ref, p)] = flip
            b = TAG[(ref, p)]
            n = 0
            for read, (label, ps, chrom) in haplotag_list(tags, autosomes_only=True).items():
                if label in ("H1", "H2") and size.get((chrom, ps), 0) >= args.min_snps:
                    mask[read] = mask.get(read, 0) | b
                    n += 1
            print("%s %s: %d reads tagged" % (ref, p, n), file=sys.stderr)
        del truth_gt
    for ref in REFS:
        b = PRES[ref]
        for read, *_ in read_rows(getattr(args, "reads_%s" % ref.lower())):
            mask[read] = mask.get(read, 0) | b
    shared = {r: m for r, m in mask.items() if m & 3 == 3}
    del mask
    print("%d reads in both read tables" % len(shared), file=sys.stderr)

    def tagged(m, ref, p):
        return bool(m & TAG[(ref, p)])

    rescued = {p: {r for r, m in shared.items() if tagged(m, "T2T", p) and not tagged(m, "GRCh38", p)}
               for p in phasers}
    lost = {p: {r for r, m in shared.items() if tagged(m, "GRCh38", p) and not tagged(m, "T2T", p)}
            for p in phasers}
    detail = set().union(*rescued.values(), *lost.values())

    # UpSet membership: reads untagged on GRCh38 by at least one phaser
    combos = {}
    for r, m in shared.items():
        u = tuple(int(not tagged(m, "GRCh38", p)) for p in phasers)
        if not any(u):
            continue
        k = u + tuple(int(tagged(m, "T2T", p)) for p in phasers)
        combos[k] = combos.get(k, 0) + 1
    with open(args.out + ".upset.tsv", "w") as fh:
        fh.write("\t".join(["untagged_GRCh38_" + p for p in phasers] +
                           ["tagged_T2T_" + p for p in phasers] + ["reads"]) + "\n")
        for k, n in sorted(combos.items(), key=lambda kv: -kv[1]):
            fh.write("\t".join(str(x) for x in k) + "\t%d\n" % n)

    # strata, benchmark and paralog sets per reference
    strata, bench, paralog = {}, {}, {}
    for ref in REFS:
        d = getattr(args, "strata_%s" % ref.lower())
        sets = []
        if d:
            for s in args.stratum:
                path = os.path.join(d, "%s_%s.bed" % (ref, s))
                if os.path.exists(path):
                    sets.append((s, IntervalSet.from_bed(path, autosomes_only=True)))
        strata[ref] = sets
        by = dict(sets)
        rows = []
        for n in ("notinalldifficultregions", "alldifficultregions"):
            if n in by:
                for c in by[n].starts:
                    rows += list(zip([c] * len(by[n].starts[c]), by[n].starts[c], by[n].ends[c]))
        bench[ref] = IntervalSet(rows)
        rows = []
        for n in ("segdups", "lowmappabilityall"):
            if n in by:
                for c in by[n].starts:
                    rows += list(zip([c] * len(by[n].starts[c]), by[n].starts[c], by[n].ends[c]))
        paralog[ref] = IntervalSet(rows)

    # positions of the detail reads, and the strata of every shared read (the background)
    pos = {ref: {} for ref in REFS}
    background = {ref: {} for ref in REFS}
    for ref in REFS:
        for read, chrom, start, end, mapq in read_rows(getattr(args, "reads_%s" % ref.lower())):
            if read not in shared:
                continue
            for s in strata_of(strata[ref], chrom, start, end):
                background[ref][s] = background[ref].get(s, 0) + 1
            if read in detail:
                pos[ref][read] = (chrom, start, end, mapq)
    nsnp = {ref: snps_of(getattr(args, "alleles_%s" % ref.lower()), detail) for ref in REFS}
    truth = {ref: truth_of(getattr(args, "truth_tags_%s" % ref.lower()), detail) for ref in REFS}

    # the label and block of the detail reads in every run, for the correctness
    run_tag = {}
    for (ref, p), (tags, _) in runs.items():
        run_tag[(ref, p)] = {r: v for r, v in haplotag_list(tags, autosomes_only=True).items() if r in detail}

    def correct(read, ref, p):
        t = run_tag[(ref, p)].get(read)
        if t is None or read not in truth[ref] or t[0] not in ("H1", "H2"):
            return None
        f = flips[(ref, p)].get((t[2], t[1]))
        if f is None:
            return None
        return int(((0 if t[0] == "H1" else 1) ^ f) == truth[ref][read])

    summary = []
    total = len(shared)
    for direction, sets_, on, off in (("rescued", rescued, "T2T", "GRCh38"), ("lost", lost, "GRCh38", "T2T")):
        for p in phasers:
            reads = sets_[p]
            n = len(reads)
            summary.append((p, direction, "total", "reads", n, n / total if total else 0, ""))
            cnt = {}
            for r in reads:
                k = reason(pos[off].get(r), nsnp[off].get(r, 0), paralog[off], bench[off])
                cnt[k] = cnt.get(k, 0) + 1
            for k in REASONS:
                summary.append((p, direction, "reason_on_" + off, k, cnt.get(k, 0),
                                cnt.get(k, 0) / n if n else 0, ""))
            st_cnt = {}
            for r in reads:
                if r in pos[on]:
                    for s in strata_of(strata[on], *pos[on][r][:3]):
                        st_cnt[s] = st_cnt.get(s, 0) + 1
            for s, _ in strata[on]:
                summary.append((p, direction, "stratum_on_" + on, s, st_cnt.get(s, 0),
                                st_cnt.get(s, 0) / n if n else 0,
                                "%.6f" % (background[on].get(s, 0) / total if total else 0)))
            ok = [correct(r, on, p) for r in reads]
            ok = [x for x in ok if x is not None]
            summary.append((p, direction, "tag_on_" + on, "compared", len(ok), 1.0, ""))
            summary.append((p, direction, "tag_on_" + on, "correct", sum(ok),
                            sum(ok) / len(ok) if ok else 0, ""))
            print("%s %s: %d reads (%.2f%% of shared), %d of %d tags right on %s" % (
                p, direction, n, 100 * n / total if total else 0, sum(ok), len(ok), on), file=sys.stderr)

    with open(args.out + ".summary.tsv", "w") as fh:
        fh.write("phaser\tdirection\tcategory\tname\treads\tshare\tbackground_share\n")
        for row in summary:
            fh.write("%s\t%s\t%s\t%s\t%d\t%.6f\t%s\n" % row)

    with gzip.open(args.out + ".reads.tsv.gz", "wt") as fh:
        cols = ["read", "rescued_by", "lost_by"]
        for ref in REFS:
            k = "g" if ref == "GRCh38" else "t"
            cols += ["%s_%s" % (k, c) for c in ("chrom", "start", "end", "mapq", "snps", "reason", "strata")]
        cols += ["correct_%s_%s" % (ref, p) for ref in REFS for p in phasers]
        fh.write("\t".join(cols) + "\n")
        for r in sorted(detail):
            row = [r, ",".join(p for p in phasers if r in rescued[p]) or "-",
                   ",".join(p for p in phasers if r in lost[p]) or "-"]
            for ref in REFS:
                q = pos[ref].get(r)
                n = nsnp[ref].get(r, 0)
                untagged = not any(tagged(shared[r], ref, p) for p in phasers)
                row += [q[0], str(q[1]), str(q[2]), str(q[3])] if q else ["-"] * 4
                row += [str(n), reason(q, n, paralog[ref], bench[ref]) if untagged else "tagged",
                        (",".join(strata_of(strata[ref], *q[:3])) or "-") if q else "-"]
            for ref in REFS:
                for p in phasers:
                    c = correct(r, ref, p)
                    row.append("NA" if c is None else str(c))
            fh.write("\t".join(row) + "\n")

    # regions of the rescued reads, on both references
    for ref in REFS:
        iv = []
        for r in set().union(*rescued.values()):
            q = pos[ref].get(r)
            if q:
                iv.append((q[0], q[1], q[2], tuple(p for p in phasers if r in rescued[p])))
        regions = merge(iv, args.merge_gap)
        regions.sort(key=lambda x: -x[3])
        with open("%s.regions_%s.tsv" % (args.out, ref), "w") as fh:
            fh.write("\t".join(["chrom", "start", "end", "reads"] + ["reads_" + p for p in phasers] +
                               ["strata"]) + "\n")
            for c, s, e, n, per in regions:
                touched = [name for name, ivs in strata[ref] if ivs.overlap(c, s, e) > 0
                           and name != "notinalldifficultregions"]
                fh.write("\t".join([c, str(s), str(e), str(n)] + [str(per.get(p, 0)) for p in phasers] +
                                   [",".join(touched) or "-"]) + "\n")
    print("written %s.*" % args.out, file=sys.stderr)


if __name__ == "__main__":
    main()
