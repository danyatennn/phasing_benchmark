"""Readers shared by the analysis scripts: VCFs, BEDs, haplotag lists, read
tables. Plain Python, no pysam, so they run in either environment.

Conventions, the same in every script:
  positions are 1-based in VCFs and haplotag lists, 0-based half-open in BEDs
  and in the read table (start, end), as samtools and pysam give them
  a phase block is (chromosome, PS); a phased het genotype is a|b with a != b
  "autosomes" are chr1-22
"""

import bisect
import gzip
import io
import os
import re

AUTOSOME_RE = re.compile(r"^chr([1-9]|1[0-9]|2[0-2])$")
AUTOSOMES = ["chr%d" % i for i in range(1, 23)]


def is_autosome(chrom):
    return AUTOSOME_RE.match(chrom) is not None


def open_text(path):
    """A text handle for plain or gzip/bgzip files."""
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="ascii", errors="replace")
    return open(path, encoding="ascii", errors="replace")


def write_text(path):
    """Writes through a temporary file; call .commit() on the result to keep it."""
    return AtomicWriter(path)


class AtomicWriter:
    def __init__(self, path):
        self.path = path
        self.tmp = "%s.%d.tmp" % (path, os.getpid())
        if path.endswith(".gz"):
            self.fh = io.TextIOWrapper(gzip.open(self.tmp, "wb", compresslevel=3), encoding="ascii")
        else:
            self.fh = open(self.tmp, "w")

    def write(self, text):
        self.fh.write(text)

    def commit(self):
        self.fh.close()
        os.replace(self.tmp, self.path)

    def abandon(self):
        try:
            self.fh.close()
        finally:
            if os.path.exists(self.tmp):
                os.remove(self.tmp)


# ---------------------------------------------------------------- VCF

def vcf_records(path, autosomes_only=False):
    """(chrom, pos, ref, alt, gt, ps) for every record, from the first sample.
    ref and alt are upper-cased; gt is the raw GT string, ps the raw PS
    ('.' when absent)."""
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 10)
            if len(f) < 10:
                continue
            chrom = f[0]
            if autosomes_only and not is_autosome(chrom):
                continue
            keys = f[8].split(":")
            vals = f[9].split(":")
            gt = vals[keys.index("GT")] if "GT" in keys and keys.index("GT") < len(vals) else "."
            ps = "."
            if "PS" in keys:
                i = keys.index("PS")
                if i < len(vals):
                    ps = vals[i]
            yield chrom, int(f[1]), f[3].upper(), f[4].upper(), gt, ps


def phased_het(gt):
    """(a, b) allele indices of a phased heterozygous GT like 0|1, else None."""
    if "|" not in gt:
        return None
    parts = gt.split("|")
    if len(parts) != 2 or parts[0] == parts[1] or "." in parts:
        return None
    return parts[0], parts[1]


def phase_blocks(path, autosomes_only=True):
    """(chrom, PS) -> sorted positions of the phased het genotypes in it."""
    blocks = {}
    for chrom, pos, ref, alt, gt, ps in vcf_records(path, autosomes_only):
        if ps in (".", "") or phased_het(gt) is None:
            continue
        blocks.setdefault((chrom, ps), []).append(pos)
    for v in blocks.values():
        v.sort()
    return blocks


# ---------------------------------------------------------------- BED

class IntervalSet:
    """Merged intervals per chromosome; overlap and membership by bisect."""

    def __init__(self, intervals=None):
        self.starts = {}
        self.ends = {}
        self.cum = {}
        if intervals is not None:
            self._build(intervals)

    @classmethod
    def from_bed(cls, path, autosomes_only=False):
        rows = []
        with open_text(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.split("\t", 3)
                if len(f) < 3:
                    continue
                if autosomes_only and not is_autosome(f[0]):
                    continue
                rows.append((f[0], int(f[1]), int(f[2])))
        return cls(rows)

    def _build(self, rows):
        per = {}
        for c, s, e in rows:
            if e > s:
                per.setdefault(c, []).append((s, e))
        for c, iv in per.items():
            iv.sort()
            merged = []
            for s, e in iv:
                if merged and s <= merged[-1][1]:
                    if e > merged[-1][1]:
                        merged[-1][1] = e
                else:
                    merged.append([s, e])
            self.starts[c] = [s for s, _ in merged]
            self.ends[c] = [e for _, e in merged]
            acc, cum = 0, [0]
            for s, e in merged:
                acc += e - s
                cum.append(acc)
            self.cum[c] = cum

    def bp(self):
        return sum(v[-1] for v in self.cum.values())

    def contains(self, chrom, pos0):
        """Is the 0-based position inside?"""
        st = self.starts.get(chrom)
        if not st:
            return False
        i = bisect.bisect_right(st, pos0) - 1
        return i >= 0 and pos0 < self.ends[chrom][i]

    def overlap(self, chrom, start, end):
        """Bases of [start, end) inside the set."""
        st = self.starts.get(chrom)
        if not st or end <= start:
            return 0
        en = self.ends[chrom]
        i = bisect.bisect_right(st, start) - 1
        if i < 0 or en[i] <= start:
            i += 1
        j = bisect.bisect_left(st, end) - 1
        if j < i:
            return 0
        total = self.cum[chrom][j + 1] - self.cum[chrom][i]
        total -= max(0, start - st[i])
        total -= max(0, en[j] - end)
        return max(0, total)


# --------------------------------------------------------- haplotag lists

def haplotag_list(path, autosomes_only=True):
    """read -> (label, PS, chrom) from a whatshap --output-haplotag-list file,
    the first line of a read only. label is H1, H2 or none."""
    tags = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] in tags:
                continue
            if autosomes_only and not is_autosome(f[3]):
                continue
            tags[f[0]] = (f[1], f[2], f[3])
    return tags


def truth_labels(path, source, autosomes_only=True, read_chrom=None):
    """read -> 0 or 1, the truth haplotype.
      truthtags  whatshap haplotag list against the phased truth (H1 -> 0)
      matpat     bin/matpat_truth.sh table: read, label (PAT -> 0, MAT -> 1)
      yak        bin/yak_truth.sh table, the same layout
    For matpat and yak the chromosome comes from read_chrom (the read table)."""
    labels = {}
    if source == "truthtags":
        for read, (label, _, _) in haplotag_list(path, autosomes_only).items():
            if label == "H1":
                labels[read] = 0
            elif label == "H2":
                labels[read] = 1
        return labels
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or f[0] in labels:
                continue
            if autosomes_only:
                chrom = read_chrom.get(f[0]) if read_chrom is not None else None
                if chrom is None or not is_autosome(chrom):
                    continue
            if f[1] == "PAT":
                labels[f[0]] = 0
            elif f[1] == "MAT":
                labels[f[0]] = 1
    return labels


# ------------------------------------------------------------- read table

READ_TABLE_COLUMNS = ["read", "chrom", "start", "end", "length", "mapq", "aligned",
                      "softclip", "nm", "ins", "dels", "mean_bq"]


def read_table(path, autosomes_only=True):
    """read -> (chrom, start, end, length, mapq, aligned, softclip, nm, ins,
    dels, mean_bq) from plots/read_table.py output. nm is -1 when the record
    had no NM tag, mean_bq NaN without base qualities."""
    reads = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header != READ_TABLE_COLUMNS:
            raise ValueError("%s: unexpected header %s" % (path, header))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if autosomes_only and not is_autosome(f[1]):
                continue
            reads[f[0]] = (f[1], int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]),
                           int(f[7]), int(f[8]), int(f[9]), int(f[10]), float(f[11]))
    return reads


def read_tsv_header(path):
    """A list of dicts from a TSV with a header line."""
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]
