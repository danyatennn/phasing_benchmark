#!/usr/bin/env bash
# The mechanism table of the difficult regions (plots/strata_mechanism.py): per
# stratum, the allele balance and paralog-like links of the truth SNPs, depth
# and low-MAPQ share, and per run how strongly switch errors and block ends
# gather there.
#
#   metrics_strata_auto/mechanism_<ref>.tsv
#
# Needs: bin/allele_matrix.sh (truth input), bin/strata_coverage.sh,
# bin/switch_bed.sh, and the strata cut to chr1-22 (bin/prepare_strata.sh).
# Rebuilt every time; it only reads.
#
#   PHASING_DEPTHS=30 bash bin/strata_mechanism.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

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
    printf 'cell\t%s\t%s\t%s\t%s\t%s\n' "${platform}" "${depth}" \
        "${READS_DIR}/${cell}_truth.sites.tsv.gz" "${READS_DIR}/${cell}_truth.links.tsv.gz" \
        "${METRICS_STRATA_DIR}/coverage_${cell}.tsv" >>"${jobs}"
    for caller in "${CALLERS[@]}"; do
      for phaser in "${PHASERS[@]}"; do
        name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
        [[ -s "${SWITCHBED_DIR}/${name}.switches.bed" || -e "${SWITCHBED_DIR}/${name}.switches.bed" ]] || continue
        printf 'run\t%s\t%s\t%s\t%s\t%s\t%s\n' "${platform}" "${depth}" "${caller}" "${phaser}" \
            "${SWITCHBED_DIR}/${name}.switches.bed" "${PHASE_DIR}/${name}.vcf.gz" >>"${jobs}"
      done
    done
  done
done

out="${METRICS_STRATA_DIR}/mechanism_${REFERENCE_NAME}.tsv"
log "${SAMPLE} ${REFERENCE_NAME}: mechanism table from $(grep -c '^cell' "${jobs}") cells and $(grep -c '^run' "${jobs}") runs"
"${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/strata_mechanism.py" \
    --reference "${REFERENCE_NAME}" --sample "${SAMPLE}" --truth-vcf "${TRUTH_PHASED}" \
    "${strata[@]}" --min-snps "${READACC_MIN_SNPS}" --jobs "${jobs}" --out "${out}" \
  2>>"${LOGS_DIR}/regions/mechanism_${REFERENCE_NAME}.log" \
  || die "mechanism table failed, see ${LOGS_DIR}/regions/mechanism_${REFERENCE_NAME}.log"
log "mechanism: ${out}"
