#!/usr/bin/env bash
# VAF histograms of the het SNPs per cell, genome-wide and per stratum
# (plots/vaf_strata.py), for the VAF panels of plots/figures_regions.py.
#
#   metrics_strata_auto/vaf_<ref>.tsv
#
# Needs bin/allele_matrix.sh (reads/<cell>_<input>.sites.tsv.gz) and the strata
# cut to chr1-22 (bin/prepare_strata.sh). Rebuilt every time; it only reads.
#
#   PHASING_DEPTHS=30 bash bin/vaf_strata.sh
#   PHASING_DEPTHS=30 VAF_INPUTS="truth clair3" REFERENCE_NAME=T2T bash bin/vaf_strata.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

VAF_INPUTS=(${VAF_INPUTS:-truth})
VAF_MIN_READS="${VAF_MIN_READS:-5}"
VAF_THREADS="${VAF_THREADS:-8}"

mkdir -p "${METRICS_STRATA_DIR}" "${TMP_DIR}" "${LOGS_DIR}/regions"
jobs="$(mktemp -p "${TMP_DIR}")"
trap 'rm -f "${jobs}"' EXIT

strata=()
for s in $(STRATA_TIER=1 strata_names); do
  strata+=(--stratum "${s}=$(stratum_bed "${s}")")
done
[[ ${#strata[@]} -gt 0 ]] || die "no strata in ${STRATA_LIST}"

for platform in "${PLATFORMS[@]}"; do
  for depth in "${DEPTHS[@]}"; do
    cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
    for input in "${VAF_INPUTS[@]}"; do
      sites="${READS_DIR}/${cell}_${input}.sites.tsv.gz"
      [[ -s "${sites}" ]] || { warn "${cell} ${input}: no allele matrix - run bin/allele_matrix.sh"; continue; }
      printf '%s\t%s\t%s\t%s\n' "${platform}" "${depth}" "${input}" "${sites}" >>"${jobs}"
    done
  done
done
[[ -s "${jobs}" ]] || die "no allele matrix for ${PLATFORMS[*]} x ${DEPTHS[*]} x ${VAF_INPUTS[*]}"

out="${METRICS_STRATA_DIR}/vaf_${REFERENCE_NAME}.tsv"
log "${SAMPLE} ${REFERENCE_NAME}: VAF histograms of $(wc -l <"${jobs}") cells"
"${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/vaf_strata.py" \
    --reference "${REFERENCE_NAME}" --sample "${SAMPLE}" "${strata[@]}" \
    --min-reads "${VAF_MIN_READS}" --threads "${VAF_THREADS}" --jobs "${jobs}" --out "${out}" \
  2>>"${LOGS_DIR}/regions/vaf_${REFERENCE_NAME}.log" \
  || die "VAF histograms failed, see ${LOGS_DIR}/regions/vaf_${REFERENCE_NAME}.log"
log "VAF: ${out}"
