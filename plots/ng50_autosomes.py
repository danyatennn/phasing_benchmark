#!/usr/bin/env python3
"""Phase block NG50 over the 22 autosomes, one genome size per reference.

whatshap stats --chr-lengths divides by the length of every chromosome the VCF
has a record on. The GRCh38 truth has het SNPs in the chrX PAR, the T2T truth
has none on chrX, so on GRCh38 chrX's 156 Mb join the genome and NG50 comes out
6-15% lower than over the same blocks on T2T - a difference between the truth
sets, not between the references. Here every run is measured against the
autosomes of its reference and nothing else.

A block is what whatshap stats counts: phased het genotypes (GT with '|') that
share a PS, two or more of them; its span is last minus first position.
check_results.py --only ng50 showed this matches whatshap to 0.0% on T2T.

  python plots/ng50_autosomes.py                 every run of the grid, 8 at once
  python plots/ng50_autosomes.py -j 16 --refs T2T
  python plots/ng50_autosomes.py --phasers "lpmod mpwhatshap" --depths 30   meth_run.sh runs

metrics/<name>.ng50.tsv, one line: reference caller platform depth phaser
  ng50 blocks genome_bp. A file newer than its VCF is kept.
"""

import argparse
import gzip
import os
import sys
from multiprocessing import Pool

ROOT = "/autofs/bal36/dten/phasing_benchmark"
SAMPLE, REGION = "HG002", "wgs"
REFERENCES = ["GRCh38", "T2T"]
CALLERS = ["truth", "clair3", "deepvariant", "longcalld"]
PLATFORMS = ["ont4k", "ont5k", "hifi", "revio"]
DEPTHS = [10, 20, 30, 40, 50]
PHASERS = ["whatshap", "longphase", "margin", "hapcut2", "hiphase"]
PACBIO_ONLY = {"hiphase", "mphiphase"}
AUTOSOMES = {"chr%d" % i for i in range(1, 23)}


def genome_size(root, ref):
    total = 0
    with open(os.path.join(root, "data", "chrom_lengths_%s.tsv" % ref)) as fh:
        for line in fh:
            f = line.split("\t")
            if f[0] in AUTOSOMES:
                total += int(f[1])
    return total


def block_spans(vcf):
    blocks = {}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 10)
            if f[0] not in AUTOSOMES:
                continue
            keys = f[8].split(":")
            if "PS" not in keys or "GT" not in keys:
                continue
            vals = f[9].rstrip("\n").split(":")
            i_ps = keys.index("PS")
            gt = vals[keys.index("GT")]
            ps = vals[i_ps] if i_ps < len(vals) else "."
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
                if pos < b[0]:
                    b[0] = pos
                if pos > b[1]:
                    b[1] = pos
                b[2] += 1
    return [b - a for a, b, n in blocks.values() if n >= 2]


def ng50(spans, genome):
    acc = 0
    for s in sorted(spans, reverse=True):
        acc += s
        if acc >= genome / 2:
            return s
    return 0


def one(job):
    root, sample, ref, caller, platform, depth, phaser, genome = job
    name = "%s_%s_%s_%dx_%s_%s_%s" % (platform, ref, sample, depth, REGION, caller, phaser)
    vcf = os.path.join(root, "phased", name + ".vcf.gz")
    out = os.path.join(root, "metrics", name + ".ng50.tsv")
    if not os.path.exists(vcf + ".tbi"):
        return "no VCF   " + name
    if os.path.exists(out) and os.path.getmtime(out) > os.path.getmtime(vcf):
        return None
    spans = block_spans(vcf)
    tmp = "%s.%d" % (out, os.getpid())
    with open(tmp, "w") as fh:
        fh.write("%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\n" % (
            ref, caller, platform, depth, phaser, ng50(spans, genome), len(spans), genome))
    os.replace(tmp, out)
    return "done     " + name


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default=ROOT)
    parser.add_argument("--sample", default=SAMPLE,
                        help="HG002 lives at the root, any other sample under samples/<sample>/")
    parser.add_argument("--refs", default=" ".join(REFERENCES))
    parser.add_argument("--phasers", default=" ".join(PHASERS))
    parser.add_argument("--callers", default=" ".join(CALLERS))
    parser.add_argument("--depths", default=" ".join(str(d) for d in DEPTHS))
    parser.add_argument("-j", type=int, default=8)
    ctx = parser.parse_args()
    phasers, callers = ctx.phasers.split(), ctx.callers.split()
    depths = [int(d) for d in ctx.depths.split()]
    jobs = []
    for ref in ctx.refs.split():
        root = ctx.root if ctx.sample == "HG002" else os.path.join(ctx.root, "samples", ctx.sample)
        if not os.path.exists(os.path.join(root, "data", "chrom_lengths_%s.tsv" % ref)):
            print("%s %s: no chrom_lengths (eval_run.sh writes it) - skipped" % (ctx.sample, ref))
            continue
        genome = genome_size(root, ref)
        print("%s: autosomes %.1f Mb" % (ref, genome / 1e6))
        jobs += [(root, ctx.sample, ref, c, pl, d, ph, genome)
                 for c in callers for pl in PLATFORMS for d in depths for ph in phasers
                 if ph not in PACBIO_ONLY or pl in ("hifi", "revio")]
    with Pool(ctx.j) as pool:
        for msg in pool.imap_unordered(one, jobs):
            if msg:
                print(msg, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
