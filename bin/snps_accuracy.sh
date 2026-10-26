#!/usr/bin/env bash
# Het SNPs per read against the haplotag accuracy (plots/snps_accuracy.py), for
# every run of every cell, for the figR3 panels of plots/figures_reads.py.
#
#   metrics_reads/<cell>_<input>.snps_accuracy.tsv
#   metrics_reads/summary_snps_accuracy.tsv      all cells, one header
#
# Needs, per cell: the allele matrix of the input (bin/allele_matrix.sh), the
# haplotag lists and the truth haplotag list (eval_run.sh). A table newer than
# all its inputs is kept.
#
#   PHASING_DEPTHS=30 bash bin/snps_accuracy.sh
#   PHASING_DEPTHS=30 PHASING_CALLERS="truth clair3" REFERENCE_NAME=T2T bash bin/snps_accuracy.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

SNPACC_JOBS="${SNPACC_JOBS:-4}"
mkdir -p "${METRICS_READS_DIR}" "${TMP_DIR}" "${LOGS_DIR}/reads"

snpacc_one() {
  set -euo pipefail
  local caller="$1" platform="$2" depth="$3"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local alleles="${READS_DIR}/${cell}_${caller}"
  local truth_tags="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truthtags).tsv"
  local out="${METRICS_READS_DIR}/${cell}_${caller}.snps_accuracy.tsv"
  for f in "${alleles}.reads.tsv.gz" "${truth_tags}"; do
    [[ -s "${f}" ]] || { warn "${cell} ${caller}: missing ${f}"; return 0; }
  done

  local jobs; jobs="$(mktemp -p "${TMP_DIR}")"
  local fresh=1 phaser name tags vcf
  newer_than "${out}" "${alleles}.reads.tsv.gz" && newer_than "${out}" "${truth_tags}" || fresh=0
  for phaser in ${PHASER_LIST}; do
    phaser_applies "${phaser}" "${platform}" || continue
    name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
    tags="${METRICS_DIR}/${name}.haplotag.tsv"
    vcf="${PHASE_DIR}/${name}.vcf.gz"
    [[ -s "${tags}" && -s "${vcf}" && ! -d "${METRICS_DIR}/${name}.row.tsv.lock" ]] || continue
    newer_than "${out}" "${tags}" && newer_than "${out}" "${vcf}" || fresh=0
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${REFERENCE_NAME}" "${SAMPLE}" \
        "${caller}" "${platform}" "${depth}" "${phaser}" "${tags}" "${vcf}" >>"${jobs}"
  done
  if [[ ! -s "${jobs}" || "${fresh}" == 1 ]]; then
    rm -f "${jobs}"
    log "skip ${cell} ${caller}"
    return 0
  fi
  log "${cell} ${caller}: SNPs per read and accuracy for $(wc -l <"${jobs}") runs"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/snps_accuracy.py" \
      --alleles "${alleles}" --truth-tags "${truth_tags}" --truth-vcf "${TRUTH_PHASED}" \
      --min-snps "${READACC_MIN_SNPS}" --jobs "${jobs}" --out "${out}" \
    2>>"${LOGS_DIR}/reads/${cell}_${caller}.snps_accuracy.log" \
    || { rm -f "${jobs}"; die "SNPs per read failed for ${cell} ${caller}, see ${LOGS_DIR}/reads/${cell}_${caller}.snps_accuracy.log"; }
  rm -f "${jobs}"
}

PHASER_LIST="${PHASERS[*]}"
export -f snpacc_one bam_key vcf_name newer_than phaser_applies log warn die
export METRICS_DIR METRICS_READS_DIR PHASE_DIR READS_DIR LOGS_DIR TMP_DIR REFERENCE_NAME \
       REGION_NAME SAMPLE TRUTH_PHASED READACC_MIN_SNPS PLOTS_PYTHON PROJECT_ROOT PHASER_LIST

log "${SAMPLE} ${REFERENCE_NAME}: SNPs per read for ${CALLERS[*]} x ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${SNPACC_JOBS}" snpacc_one {1} {2} {3} \
    ::: "${CALLERS[@]}" ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true

summary="${METRICS_READS_DIR}/summary_snps_accuracy.tsv"
{
  files=("${METRICS_READS_DIR}"/*.snps_accuracy.tsv)
  [[ -e "${files[0]}" ]] && head -1 "${files[0]}"
  for f in "${files[@]}"; do
    [[ -e "${f}" ]] && tail -n +2 "${f}"
  done
  true
} >"${summary}.$$"
mv -f "${summary}.$$" "${summary}"
log "SNPs per read: ${summary}"
