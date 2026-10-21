#!/usr/bin/env python3
"""Which allele each read carries at each heterozygous SNP of an input VCF, and
what follows from that, for one cell (BAM) and one input (the truth or a
caller's het SNPs):

  <out>.sites.tsv.gz     one line per SNP: depth, reads with REF / ALT / another
                         base / a deletion there, VAF = ALT / (REF + ALT)
  <out>.links.tsv.gz     one line per pair of neighbouring SNPs: reads that
                         cover both with REF or ALT, how many show each of the
                         four allele combinations, and how many combinations
                         are really there (at least 2 reads and 10% of the
                         spanning reads). A het pair shows 2; more points to
                         reads from a paralogous copy mapped here.
  <out>.reads.tsv.gz     one line per read: SNPs it covers, with REF / ALT /
                         other / deletion, mean |VAF - 0.5| over the SNPs it
                         covers with REF or ALT, links it spans and how many of
                         them look paralogous
  <out>.votes.<phaser>.tsv.gz
                         one line per read and phase block of that phaser's
                         VCF: SNPs where the read carries the allele of
                         haplotype 1 of the block (h1), of haplotype 2 (h2)

The alleles come from samtools mpileup with --output-QNAME at the SNP positions
only, every read (MAPQ 0 included, so the reads of a paralog are seen too),
primary alignments only, no base-quality cut, no BAQ. Positions with more than
one record in the input VCF are dropped. Only chr1-22.

  python plots/allele_matrix.py --bam input.bam --sites truth.vcf.gz \\
      --phased whatshap=run1.vcf.gz --phased margin=run2.vcf.gz \\
      --samtools /path/samtools --out reads/<cell>_truth --threads 8
"""

import argparse
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchio import AUTOSOMES, is_autosome, phased_het, vcf_records, write_text

REF, ALT, OTHER, DEL = 0, 1, 2, 3
MIN_PATTERN_READS = 2
MIN_PATTERN_FRACTION = 0.10

SITES_HEADER = "chrom\tpos\tref\talt\tdepth\tn_ref\tn_alt\tn_other\tn_del\tvaf\n"
LINKS_HEADER = "chrom\tpos1\tpos2\tspanning\tc00\tc01\tc10\tc11\tn_patterns\tparalog\n"
READS_HEADER = ("read\tchrom\tn_sites\tn_ref\tn_alt\tn_other\tn_del\tvaf_dev\tn_vaf"
                "\tlinks\tparalog_links\n")
VOTES_HEADER = "read\tchrom\tps\th1\th2\n"


def load_sites(path):
    """chrom -> [(pos, ref, alt)], biallelic SNPs, positions with more than one
    record dropped."""
    per = {}
    count = {}
    for chrom, pos, ref, alt, gt, ps in vcf_records(path, autosomes_only=True):
        if len(ref) != 1 or len(alt) != 1 or alt in (".", "*"):
            continue
        per.setdefault(chrom, []).append((pos, ref, alt))
        count[(chrom, pos)] = count.get((chrom, pos), 0) + 1
    out = {}
    dropped = 0
    for chrom, sites in per.items():
        kept = []
        for pos, ref, alt in sorted(sites):
            if count[(chrom, pos)] > 1:
                dropped += 1
                continue
            kept.append((pos, ref, alt))
        out[chrom] = kept
    return out, dropped


def load_phased(path, chrom, sites_by_pos):
    """pos -> (allele code of haplotype 1, of haplotype 2, PS) for the phased
    het SNPs of one chromosome that match an input site."""
    phase = {}
    for c, pos, ref, alt, gt, ps in vcf_records(path, autosomes_only=True):
        if c != chrom or ps in (".", ""):
            continue
        site = sites_by_pos.get(pos)
        if site is None or site[0] != ref or site[1] != alt:
            continue
        ab = phased_het(gt)
        if ab is None or not set(ab) <= {"0", "1"}:
            continue
        phase[pos] = (int(ab[0]), int(ab[1]), ps)
    return phase


def parse_bases(bases):
    """The read bases of one mpileup column, one entry per read, with the
    start/end markers and inserted/deleted sequences removed."""
    out = []
    i, n = 0, len(bases)
    while i < n:
        c = bases[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            j = i + 1
            while j < n and bases[j].isdigit():
                j += 1
            i = j + int(bases[i + 1:j])
            continue
        out.append(c)
        i += 1
    return out


def allele_code(base, ref, alt):
    if base in "*#":
        return DEL
    b = base.upper()
    if b == ref:
        return REF
    if b == alt:
        return ALT
    return OTHER


def one_chromosome(job):
    chrom, sites, bam, samtools, phased_paths, tmpdir = job
    sites_by_pos = {pos: (ref, alt) for pos, ref, alt in sites}
    index_of = {pos: k for k, (pos, _, _) in enumerate(sites)}
    phases = {name: load_phased(path, chrom, sites_by_pos) for name, path in phased_paths}

    bed = os.path.join(tmpdir, "%s.bed" % chrom)
    with open(bed, "w") as fh:
        for pos, _, _ in sites:
            fh.write("%s\t%d\t%d\n" % (chrom, pos - 1, pos))

    cmd = [samtools, "mpileup", "-B", "-Q", "0", "-q", "0", "-d", "100000",
           "--ff", "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY",
           "--output-QNAME", "-l", bed, "-r", chrom, bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    site_lines, link_lines = [], []
    reads = {}            # read -> [n_sites, n_ref, n_alt, n_other, n_del, vafdev_sum, n_vaf, links, paralog]
    last = {}             # read -> (site index, allele) at the last site it showed REF or ALT
    votes = {name: {} for name in phases}
    bad_columns = 0

    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        pos = int(f[1])
        k = index_of.get(pos)
        if k is None:
            continue
        ref, alt = sites_by_pos[pos]
        depth = int(f[3])
        if depth == 0:
            continue
        bases = parse_bases(f[4])
        names = f[6].split(",")
        if len(bases) != depth or len(names) != depth:
            bad_columns += 1
            continue

        alleles = [allele_code(b, ref, alt) for b in bases]
        n = [0, 0, 0, 0]
        for a in alleles:
            n[a] += 1
        vaf = n[ALT] / (n[REF] + n[ALT]) if n[REF] + n[ALT] else float("nan")
        site_lines.append("%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n" % (
            chrom, pos, ref, alt, depth, n[REF], n[ALT], n[OTHER], n[DEL],
            "nan" if vaf != vaf else "%.4f" % vaf))

        # the link to the previous site is complete once this column is read:
        # every read covering both shows up here
        pattern = [0, 0, 0, 0]
        spanning = []
        for name, a in zip(names, alleles):
            if a > ALT:
                continue
            prev = last.get(name)
            if prev is not None and prev[0] == k - 1:
                pattern[prev[1] * 2 + a] += 1
                spanning.append(name)
        if spanning:
            total = len(spanning)
            real = sum(1 for c in pattern
                       if c >= MIN_PATTERN_READS and c >= MIN_PATTERN_FRACTION * total)
            paralog = 1 if real > 2 else 0
            link_lines.append("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (
                chrom, sites[k - 1][0], pos, total, pattern[0], pattern[1], pattern[2],
                pattern[3], real, paralog))
            for name in spanning:
                acc = reads[name]
                acc[7] += 1
                acc[8] += paralog

        vafdev = abs(vaf - 0.5) if vaf == vaf else None
        for name, a in zip(names, alleles):
            acc = reads.get(name)
            if acc is None:
                acc = reads[name] = [0, 0, 0, 0, 0, 0.0, 0, 0, 0]
            acc[0] += 1
            acc[1 + a] += 1
            if a <= ALT:
                last[name] = (k, a)
                if vafdev is not None:
                    acc[5] += vafdev
                    acc[6] += 1
                for phaser, phase in phases.items():
                    hp = phase.get(pos)
                    if hp is None:
                        continue
                    key = (name, hp[2])
                    v = votes[phaser].get(key)
                    if v is None:
                        v = votes[phaser][key] = [0, 0]
                    if a == hp[0]:
                        v[0] += 1
                    elif a == hp[1]:
                        v[1] += 1

    err = proc.stderr.read()
    if proc.wait() != 0:
        raise RuntimeError("samtools mpileup failed on %s %s: %s" % (bam, chrom, err[-2000:]))
    os.remove(bed)

    read_lines = []
    for name, a in reads.items():
        read_lines.append("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n" % (
            name, chrom, a[0], a[1], a[2], a[3], a[4],
            "%.4f" % (a[5] / a[6]) if a[6] else "nan", a[6], a[7], a[8]))
    vote_lines = {}
    for phaser, v in votes.items():
        vote_lines[phaser] = ["%s\t%s\t%s\t%d\t%d\n" % (name, chrom, ps, h1, h2)
                              for (name, ps), (h1, h2) in v.items()]
    return chrom, site_lines, link_lines, read_lines, vote_lines, bad_columns


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sites", required=True, help="the input VCF: truth or caller het SNPs")
    parser.add_argument("--phased", action="append", default=[],
                        help="phaser=path of a phased VCF over the same input, repeatable")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--out", required=True, help="output prefix")
    parser.add_argument("--threads", type=int, default=8, help="chromosomes at once")
    args = parser.parse_args()

    phased = []
    for item in args.phased:
        name, _, path = item.partition("=")
        if not name or not path or not os.path.exists(path):
            sys.exit("--phased %s: expected phaser=existing.vcf.gz" % item)
        phased.append((name, path))

    sites, dropped = load_sites(args.sites)
    chroms = [c for c in AUTOSOMES if sites.get(c)]
    print("%s: %d SNPs on %d autosomes, %d positions with several records dropped" % (
        args.sites, sum(len(sites[c]) for c in chroms), len(chroms), dropped), file=sys.stderr)

    outs = {
        "sites": write_text(args.out + ".sites.tsv.gz"),
        "links": write_text(args.out + ".links.tsv.gz"),
        "reads": write_text(args.out + ".reads.tsv.gz"),
    }
    for name, _ in phased:
        outs["votes." + name] = write_text("%s.votes.%s.tsv.gz" % (args.out, name))
    outs["sites"].write(SITES_HEADER)
    outs["links"].write(LINKS_HEADER)
    outs["reads"].write(READS_HEADER)
    for name, _ in phased:
        outs["votes." + name].write(VOTES_HEADER)

    tmpdir = tempfile.mkdtemp(prefix="allele_matrix.", dir=os.path.dirname(os.path.abspath(args.out)))
    bad = 0
    try:
        jobs = [(c, sites[c], args.bam, args.samtools, phased, tmpdir) for c in chroms]
        with Pool(max(1, min(args.threads, len(jobs)))) as pool:
            # imap keeps chromosome order and hands each one over as soon as
            # it and those before it are done, so only a few are held at once
            for chrom, s, l, r, v, b in pool.imap(one_chromosome, jobs):
                outs["sites"].write("".join(s))
                outs["links"].write("".join(l))
                outs["reads"].write("".join(r))
                for name, lines in v.items():
                    outs["votes." + name].write("".join(lines))
                bad += b
    except BaseException:
        for o in outs.values():
            o.abandon()
        raise
    finally:
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass
    for o in outs.values():
        o.commit()
    if bad:
        print("WARNING: %d pileup columns did not parse and were skipped" % bad, file=sys.stderr)
    print("%s.*: done" % args.out, file=sys.stderr)


if __name__ == "__main__":
    main()
