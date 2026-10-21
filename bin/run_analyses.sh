#!/usr/bin/env bash
# Every analysis that reads finished runs, in the order they depend on each
# other, for one sample (SAMPLE, default HG002) and one reference
# (REFERENCE_NAME, default GRCh38). Each step skips what is up to date, so the
# script can be run again after more runs have finished.
#
#   REFERENCE_NAME=T2T bash bin/run_analyses.sh
#   SAMPLE=HG005 bash bin/run_analyses.sh
#
# Steps over the whole grid of the sample: autosome metrics, read accuracy (by
# SNPs and by anchor read), read tables, tagged bases, caller FNR/FDR, NG50,
# strata cut to chr1-22 and the metrics in them. Steps at ANALYSIS_DEPTH only
# (default 30): allele matrices, read classes, switch positions, stratum
# coverage, the mechanism table.

set -euo pipefail

cd /autofs/bal36/dten/phasing_benchmark

ANALYSIS_DEPTH="${ANALYSIS_DEPTH:-30}"
PLOTS_PY="/autofs/bal36/dten/miniforge3/envs/plots/bin/python"
export SAMPLE="${SAMPLE:-HG002}"
export REFERENCE_NAME="${REFERENCE_NAME:-GRCh38}"

step() {
  printf '[%s] ===== %s %s: %s\n' "$(date '+%F %H:%M:%S')" "${SAMPLE}" "${REFERENCE_NAME}" "$*" >&2
  "$@"
}

step bash bin/autosomes.sh
step bash readacc_run.sh
step bash bin/read_tables.sh
step bash bin/tagged_bases.sh
step bash bin/calls_fdr.sh
step bash readacc_anchor.sh
step "${PLOTS_PY}" plots/ng50_autosomes.py --sample "${SAMPLE}" --refs "${REFERENCE_NAME}" -j 16
step bash bin/prepare_strata.sh
step bash regions.sh
step bash regions_statistics.sh
step env PHASING_DEPTHS="${ANALYSIS_DEPTH}" bash bin/allele_matrix.sh
step env PHASING_DEPTHS="${ANALYSIS_DEPTH}" bash bin/read_classes.sh
step env PHASING_DEPTHS="${ANALYSIS_DEPTH}" bash bin/switch_bed.sh
step env PHASING_DEPTHS="${ANALYSIS_DEPTH}" bash bin/strata_coverage.sh
step env PHASING_DEPTHS="${ANALYSIS_DEPTH}" bash bin/strata_mechanism.sh
printf '[%s] ===== %s %s: analyses done\n' "$(date '+%F %H:%M:%S')" "${SAMPLE}" "${REFERENCE_NAME}" >&2
