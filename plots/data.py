"""Reads what eval_run.sh and regions.sh have written so far.

The grid is unfinished while this runs, so a run is one of
  done     its row.tsv is there, whole, and says what its name says
  locked   <row>.lock exists - some host is writing it right now
  missing  anything else, including a row left by an older eval_run.sh
Only "done" runs reach a figure. Nothing here raises on a bad or absent file.
"""

import math
import os

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

# readacc_run.sh: read-level accuracy with each block oriented by its variants
READACC_COLUMNS = [
    "reference", "region", "caller", "platform", "depth", "phaser",
    "reads_truth_tagged", "reads_tagged", "reads_compared", "reads_correct",
    "reads_unoriented", "read_precision", "read_recall",
    "blocks_oriented", "blocks_flipped",
]

STRATA_COLUMNS = [
    "reference", "region", "caller", "platform", "depth", "phaser", "stratum",
    "variants", "phased", "phased_fraction", "blocks", "block_n50",
    "switches", "switch_rate", "hamming", "hamming_rate",
    "covered", "assessed", "truth_total", "emitted", "shared", "correct",
    "fnr", "fdr", "recall_phased", "precision_phased",
]

LABEL_COLUMNS = {"reference", "region", "caller", "platform", "phaser", "stratum"}

AUTOSOMES = {"chr%d" % i for i in range(1, 23)}

# bin/prepare_strata.sh, TIER1
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
    """A whatshap TSV: one header line that starts with '#', then rows."""
    try:
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            header[0] = header[0].lstrip("#")
            return [dict(zip(header, line.rstrip("\n").split("\t")))
                    for line in fh if line.strip()]
    except (OSError, UnicodeDecodeError):
        return []


class Metrics:
    def __init__(self, root, sample, region):
        self.root = root
        self.metrics_dir = os.path.join(root, "metrics")
        self.strata_dir = os.path.join(root, "metrics_strata")
        self.sample = sample
        self.region = region
        self.rows = {}
        self.status = {}
        self._compare = {}
        self._strata = {}

    def name(self, ref, caller, platform, depth, phaser):
        return "%s_%s_%s_%sx_%s_%s_%s" % (
            platform, ref, self.sample, depth, self.region, caller, phaser)

    def load(self, references, callers, platforms, depths, phasers):
        for ref in references:
            for caller in callers:
                for platform in platforms:
                    for depth in depths:
                        for phaser in phasers:
                            key = (ref, caller, platform, depth, phaser)
                            self.status[key] = self._load_one(key)
        return self

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
        # read_precision and read_recall come from readacc_run.sh, when it has
        # scored this run; until then they read as NaN
        acc = _read_row(os.path.join(self.metrics_dir, self.name(*key) + ".readacc.tsv"),
                        READACC_COLUMNS)
        if acc is not None and acc["phaser"] == phaser and acc["caller"] == caller:
            row.update({c: acc[c] for c in READACC_COLUMNS[6:]})
        self.rows[key] = row
        return "done"

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
        two switches in a row - one wrong site, not a wrong stretch."""
        switches = flips = assessed = 0.0
        seen = False
        for line in self.compare(key):
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
