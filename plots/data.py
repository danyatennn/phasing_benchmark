"""Reads what the pipeline has written so far, for one sample.

The grid is unfinished while this runs, so a run is one of
  done     its row.tsv is whole and says what its name says, and its metrics
           over chr1-22 (bin/autosomes.sh) are there
  no_auto  the row is whole but bin/autosomes.sh has not scored it yet
  locked   <row>.lock exists - some host is writing it right now
  missing  anything else
Only "done" runs reach a figure, so no figure mixes whole-genome numbers with
numbers over chr1-22. Nothing here raises on a bad or absent file.

What a done run carries, on top of its row.tsv (the row's own values are
replaced where the autosome table has the same column):
  <run>.auto.tsv                    switches, hamming, covered, truth_total,
                                    shared, rates, tagging - over chr1-22
  <run>.ng50.tsv                    block_ng50 over chr1-22
  <run>.tagbases_<rules>.tsv        bases_tagged_fraction
  <run>.<readacc rules>.tsv         read_precision, read_recall (blocks oriented
                                    by their SNPs)
  <run>.anchor_<truth>_<rules>.tsv  read_precision_anchor_<truth>, ...
  <run>.readacc_<truth>_<rules>.tsv read_precision_var_<truth>, ... (blocks
                                    oriented by their SNPs, bin/readacc_truths.sh)

HG002 lives at the project root; any other sample under samples/<sample>/.
"""

import math
import os

from style import applies as phaser_applies

NAN = float("nan")

ROW_COLUMNS = [
    "reference", "region", "caller", "platform", "depth", "phaser",
    "variants", "phased", "phased_fraction", "blocks", "block_n50",
    "switches", "switch_rate", "hamming", "hamming_rate", "haplotype_accuracy",
    "covered", "assessed", "truth_total", "emitted", "shared", "correct",
    "fnr", "fdr", "recall_phased", "precision_phased",
    "reads_tagged", "reads_total", "tagging_completeness",
    "reads_correct", "reads_compared", "haplotag_accuracy",
]

# bin/autosomes.sh (plots/autosomes.py)
AUTO_COLUMNS = [
    "reference", "region", "caller", "platform", "depth", "phaser",
    "variants", "phased", "blocks", "switches", "assessed", "hamming", "covered",
    "truth_total", "emitted", "shared", "correct", "switch_rate", "hamming_rate",
    "recall_phased", "fnr", "fdr", "reads_total", "reads_tagged", "tagging_completeness",
    "reads_tagged_min", "tagging_completeness_min",
]

# readacc_run.sh: read-level accuracy with each block oriented by its variants
READACC_COLUMNS = [
    "reference", "region", "caller", "platform", "depth", "phaser",
    "reads_truth_tagged", "reads_tagged", "reads_compared", "reads_correct",
    "reads_unoriented", "read_precision", "read_recall",
    "blocks_oriented", "blocks_flipped",
]

# readacc_anchor.sh (plots/readacc_anchor.py)
ANCHOR_COLUMNS = [
    "reference", "region", "caller", "platform", "depth", "phaser", "truth_source",
    "reads_truth_labelled", "reads_tagged", "reads_scored", "reads_correct",
    "anchors", "blocks_flipped", "read_precision", "read_recall",
]
ANCHOR_TRUTHS = ["truthtags", "matpat", "yak"]

# bin/readacc_truths.sh: read accuracy against each read truth, every block
# oriented by its SNPs as in readacc_run.sh; *_common over the reads that carry
# a label of this truth and of the truth VCF, the numbers to compare truths by
VARTRUTH_COLUMNS = [
    "reference", "region", "caller", "platform", "depth", "phaser", "truth_source",
    "reads_truth_tagged", "reads_tagged", "reads_compared", "reads_correct", "reads_unoriented",
    "read_precision", "read_recall", "reads_truth_common", "reads_compared_common",
    "reads_correct_common", "read_precision_common", "read_recall_common",
]

# plots/ng50_autosomes.py
NG50_COLUMNS = ["reference", "caller", "platform", "depth", "phaser",
                "block_ng50", "blocks_autosomal", "genome_bp"]

# bin/tagged_bases.sh (plots/tagged_bases.py)
TAGBASES_COLUMNS = ["reference", "caller", "platform", "depth", "phaser", "reads_total",
                    "reads_tagged", "bases_total", "bases_tagged", "bases_tagged_fraction"]

# regions.sh, per stratum
STRATA_COLUMNS = [
    "reference", "region", "caller", "platform", "depth", "phaser", "stratum",
    "variants", "phased", "phased_fraction", "blocks", "block_n50",
    "switches", "switch_rate", "hamming", "hamming_rate",
    "covered", "assessed", "truth_total", "emitted", "shared", "correct",
    "fnr", "fdr", "recall_phased", "precision_phased",
]

LABEL_COLUMNS = {"reference", "region", "caller", "platform", "phaser", "stratum",
                 "truth_source", "sample", "group_type", "group", "metric"}

AUTOSOMES = {"chr%d" % i for i in range(1, 23)}

# bin/prepare_strata.sh, TIER1 (chrX_nonPAR is empty once cut to chr1-22)
TIER1_STRATA = [
    "notinalldifficultregions", "alldifficultregions", "segdups",
    "lowmappabilityall", "AllTandemRepeatsandHomopolymers_slop5",
    "satellites_slop5", "MHC", "chrX_nonPAR",
]


def num(text):
    try:
        value = float(text)
    except (TypeError, ValueError):
        return NAN
    return value if math.isfinite(value) else NAN


def ratio(a, b):
    return a / b if b and math.isfinite(a) and math.isfinite(b) else NAN


def _exists(path):
    try:
        return os.path.exists(path)
    except OSError:
        return False


def _read_row(path, columns):
    try:
        with open(path) as fh:
            fields = fh.readline().rstrip("\n").split("\t")
    except (OSError, UnicodeDecodeError):
        return None
    if len(fields) != len(columns):
        return None
    return {c: (v if c in LABEL_COLUMNS else num(v)) for c, v in zip(columns, fields)}


def _read_table(path):
    """A TSV with one header line (a leading '#' is dropped), then rows."""
    try:
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            header[0] = header[0].lstrip("#")
            return [dict(zip(header, line.rstrip("\n").split("\t")))
                    for line in fh if line.strip()]
    except (OSError, UnicodeDecodeError):
        return []


def sample_root(root, sample):
    return root if sample == "HG002" else os.path.join(root, "samples", sample)


class Metrics:
    def __init__(self, root, sample, region, readacc_rules="auto_min2", require_auto=True):
        self.root = root
        self.sample_root = sample_root(root, sample)
        self.metrics_dir = os.path.join(self.sample_root, "metrics")
        self.reads_dir = os.path.join(self.sample_root, "metrics_reads")
        self.strata_dir = os.path.join(self.sample_root, "metrics_strata_auto")
        self.sample = sample
        self.region = region
        self.rules = readacc_rules
        self.require_auto = require_auto
        self.rows = {}
        self.status = {}
        self.calls = {}
        self._compare = {}
        self._strata = {}
        self._classes = {}
        self._mechanism = {}

    def name(self, ref, caller, platform, depth, phaser):
        return "%s_%s_%s_%sx_%s_%s_%s" % (
            platform, ref, self.sample, depth, self.region, caller, phaser)

    def load(self, references, callers, platforms, depths, phasers):
        for ref in references:
            for caller in callers:
                for platform in platforms:
                    for depth in depths:
                        for phaser in phasers:
                            if not phaser_applies(platform, phaser):
                                continue        # e.g. HiPhase on ONT: not a missing run
                            key = (ref, caller, platform, depth, phaser)
                            self.status[key] = self._load_one(key)
        self._load_calls()
        return self

    def _merge(self, row, key, suffix, columns, rename=None):
        ref, caller, platform, depth, phaser = key
        extra = _read_row(os.path.join(self.metrics_dir, self.name(*key) + suffix), columns)
        if extra is None or extra.get("phaser") != phaser or extra.get("caller") != caller:
            return False
        for c in columns[6:] if "region" in columns else columns[5:]:
            row[(rename or {}).get(c, c)] = extra[c]
        return True

    def _load_one(self, key):
        ref, caller, platform, depth, phaser = key
        path = os.path.join(self.metrics_dir, self.name(*key) + ".row.tsv")
        if _exists(path + ".lock"):
            return "locked"
        row = _read_row(path, ROW_COLUMNS)
        if row is None:
            return "missing"
        said = (row["reference"], row["caller"], row["platform"], row["depth"], row["phaser"])
        if said != (ref, caller, platform, float(depth), phaser) or row["region"] != self.region:
            return "missing"
        has_auto = self._merge(row, key, ".auto.tsv", AUTO_COLUMNS)
        if self.require_auto and not has_auto:
            return "no_auto"
        self._merge(row, key, ".%s.tsv" % self._readacc_suffix(), READACC_COLUMNS)
        self._merge(row, key, ".ng50.tsv", NG50_COLUMNS)
        self._merge(row, key, ".tagbases_%s.tsv" % self.rules, TAGBASES_COLUMNS,
                    rename={"reads_total": "tb_reads_total", "reads_tagged": "tb_reads_tagged"})
        for truth in ANCHOR_TRUTHS:
            self._merge(row, key, ".anchor_%s_%s.tsv" % (truth, self.rules), ANCHOR_COLUMNS,
                        rename={c: "%s_anchor_%s" % (c, truth) for c in ANCHOR_COLUMNS})
            self._merge(row, key, ".readacc_%s_%s.tsv" % (truth, self.rules), VARTRUTH_COLUMNS,
                        rename={c: "%s_var_%s" % (c, truth) for c in VARTRUTH_COLUMNS})
        self.rows[key] = row
        return "done"

    def _readacc_suffix(self):
        return "readacc" if self.rules == "all" else "readacc_" + self.rules

    def _load_calls(self):
        """bin/calls_fdr.sh: (ref, caller, platform, depth) -> fnr, fdr, ..."""
        for line in _read_table(os.path.join(self.metrics_dir, "summary_calls.tsv")):
            try:
                key = (line["reference"], line["caller"], line["platform"], int(line["depth"]))
            except (KeyError, ValueError):
                continue
            self.calls[key] = {k: num(v) for k, v in line.items()}

    def call(self, ref, caller, platform, depth, column):
        return self.calls.get((ref, caller, platform, depth), {}).get(column, NAN)

    def get(self, key, column):
        row = self.rows.get(key)
        return row.get(column, NAN) if row is not None else NAN

    def has(self, key):
        return key in self.rows

    def widest(self, ref, caller, platform, depth, phasers):
        """What a caller found does not depend on the phaser, but a phaser may
        drop records. The run that kept the most sites stands for the cell."""
        keys = [(ref, caller, platform, depth, p) for p in phasers]
        keys = [k for k in keys if k in self.rows and math.isfinite(self.rows[k]["emitted"])]
        return max(keys, key=lambda k: self.rows[k]["emitted"]) if keys else None

    def coverage(self, ref, column):
        """How many computed runs of ref have this column."""
        return sum(1 for k, r in self.rows.items()
                   if k[0] == ref and math.isfinite(r.get(column, NAN)))

    def count(self, ref, state):
        return sum(1 for k, s in self.status.items() if k[0] == ref and s == state)

    # per-chromosome whatshap compare, only for runs whose row is done: the row
    # is written last, so its compare.tsv is complete
    def compare(self, key):
        if key not in self.rows:
            return []
        if key not in self._compare:
            path = os.path.join(self.metrics_dir, self.name(*key) + ".compare.tsv")
            self._compare[key] = _read_table(path)
        return self._compare[key]

    def chromosome_switch_rates(self, key):
        """Switch rate of each autosome; a chromosome without a switch is 0."""
        rates = []
        for line in self.compare(key):
            if line.get("chromosome") not in AUTOSOMES:
                continue
            rate = ratio(num(line.get("all_switches")), num(line.get("all_assessed_pairs")))
            if math.isfinite(rate):
                rates.append(rate)
        return rates

    def switch_flip(self, key):
        """whatshap writes all_switchflips as 'switches/flips', where a flip is
        two switches in a row - one wrong site, not a wrong stretch. Summed
        over chr1-22."""
        switches = flips = assessed = 0.0
        seen = False
        for line in self.compare(key):
            if line.get("chromosome") not in AUTOSOMES:
                continue
            parts = str(line.get("all_switchflips", "")).split("/")
            pairs = num(line.get("all_assessed_pairs"))
            if len(parts) != 2 or not math.isfinite(pairs):
                continue
            s, f = num(parts[0]), num(parts[1])
            if not (math.isfinite(s) and math.isfinite(f)):
                continue
            switches += s
            flips += f
            assessed += pairs
            seen = True
        if not seen or assessed <= 0:
            return None
        return switches / assessed, flips / assessed

    def strata_names(self):
        """Tier 1 from regions/strata.tsv when it is there, in the order of
        bin/prepare_strata.sh."""
        listed = set()
        for line in _read_table(os.path.join(self.root, "regions", "strata.tsv")):
            if line.get("tier") == "1" and line.get("stratum"):
                listed.add(line["stratum"])
        if not listed:
            return list(TIER1_STRATA)
        return [s for s in TIER1_STRATA if s in listed] + sorted(listed - set(TIER1_STRATA))

    def stratum(self, key, stratum):
        cache_key = key + (stratum,)
        if cache_key not in self._strata:
            path = os.path.join(
                self.strata_dir, "phase_%s__%s.row.tsv" % (self.name(*key), stratum))
            row = None if _exists(path + ".lock") else _read_row(path, STRATA_COLUMNS)
            if row is not None and row["stratum"] != stratum:
                row = None
            self._strata[cache_key] = row
        return self._strata[cache_key]

    # regions_statistics.sh: what each stratum is, one row per stratum
    def strata_statistics(self, ref):
        path = os.path.join(self.strata_dir, "regions_statistics_%s.tsv" % ref)
        return {line["stratum"]: {k: (v if k in LABEL_COLUMNS or k == "category" else num(v))
                                  for k, v in line.items()}
                for line in _read_table(path) if line.get("stratum")}

    def strata_overlap(self, ref):
        """(a, b) -> Jaccard index of the two strata, both orders."""
        path = os.path.join(self.strata_dir, "jaccard_%s.tsv" % ref)
        pairs = {}
        for line in _read_table(path):
            a, b, j = line.get("stratum_a"), line.get("stratum_b"), num(line.get("jaccard"))
            if a and b and math.isfinite(j):
                pairs[(a, b)] = pairs[(b, a)] = j
        return pairs

    # bin/read_classes.sh: (stratum, group_type, group) -> numbers, per run
    def classes(self, key):
        if key not in self._classes:
            path = os.path.join(self.reads_dir, self.name(*key) + ".classes.tsv")
            table = {}
            for line in _read_table(path):
                k = (line.get("stratum"), line.get("group_type"), line.get("group"))
                table[k] = {c: (v if c in LABEL_COLUMNS else num(v)) for c, v in line.items()}
            self._classes[key] = table
        return self._classes[key]

    # bin/strata_mechanism.sh: (platform, depth, caller, phaser, stratum, metric) -> value
    def mechanism(self, ref):
        if ref not in self._mechanism:
            table = {}
            for line in _read_table(os.path.join(self.strata_dir, "mechanism_%s.tsv" % ref)):
                try:
                    k = (line["platform"], int(line["depth"]), line["caller"], line["phaser"],
                         line["stratum"], line["metric"])
                except (KeyError, ValueError):
                    continue
                table[k] = num(line["value"])
            self._mechanism[ref] = table
        return self._mechanism[ref]

    def strata_coverage(self, ref, platform, depth):
        """bin/strata_coverage.sh: stratum -> (bp, depth_mapq0, depth_mapq20)"""
        cell = "%s_%s_%s_%sx_%s" % (platform, ref, self.sample, depth, self.region)
        out = {}
        for line in _read_table(os.path.join(self.strata_dir, "coverage_%s.tsv" % cell)):
            out[line.get("stratum")] = (num(line.get("bp")), num(line.get("depth_mapq0")),
                                        num(line.get("depth_mapq20")))
        return out
