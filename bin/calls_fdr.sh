#!/usr/bin/env bash
# FNR and FDR of every caller's het SNPs, inside the confident regions and on
# chr1-22, positions with several call records left out of both sides
# (plots/calls_fdr.py). The figures of callers use these, not the fnr and fdr
# of the eval_run.sh rows.
#
#   metrics/calls_<platform>_<ref>_<sample>_<depth>x_wgs_<caller>.tsv
#   metrics/summary_calls.tsv
#
#   bash bin/calls_fdr.sh
#   REFERENCE_NAME=T2T bash bin/calls_fdr.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

require_file "${TRUTH_SRC_VCF}" "truth VCF of ${SAMPLE} on ${REFERENCE_NAME}"
require_file "${CONFIDENT_BED}" "confident BED of ${SAMPLE} on ${REFERENCE_NAME}"
mkdir -p "${TMP_DIR}" "${LOGS_DIR}/calls"

jobs="$(mktemp -p "${TMP_DIR}")"
trap 'rm -f "${jobs}"' EXIT
for caller in "${CALLERS[@]}"; do
  [[ "${caller}" == truth ]] && continue
  for platform in "${PLATFORMS[@]}"; do
    for depth in "${DEPTHS[@]}"; do
      cell="$(vcf_name "${platform}" "${depth}" "${caller}")"
      vcf="${VCF_DIR}/${cell}.het_snps.vcf.gz"
      out="${METRICS_DIR}/calls_${cell}.tsv"
      [[ -s "${vcf}" ]] || continue
      newer_than "${out}" "${vcf}" && newer_than "${out}" "${TRUTH_SRC_VCF}" && continue
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${REFERENCE_NAME}" "${SAMPLE}" "${platform}" \
          "${depth}" "${caller}" "${vcf}" "${out}" >>"${jobs}"
    done
  done
done

if [[ -s "${jobs}" ]]; then
  log "${SAMPLE} ${REFERENCE_NAME}: FNR/FDR of $(wc -l <"${jobs}") caller cells"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/calls_fdr.py" \
      --truth-vcf "${TRUTH_SRC_VCF}" --bed "${CONFIDENT_BED}" --jobs "${jobs}" \
    2>>"${LOGS_DIR}/calls/calls_fdr_${REFERENCE_NAME}.log" \
    || die "calls_fdr failed, see ${LOGS_DIR}/calls/calls_fdr_${REFERENCE_NAME}.log"
else
  log "${SAMPLE} ${REFERENCE_NAME}: every caller cell is up to date"
fi

{
  printf 'reference\tsample\tplatform\tdepth\tcaller\ttruth_het_snps\tcalls\tshared\tfn\tfp'
  printf '\tfnr\tfdr\tcalls_outside_bed\tcalls_off_autosomes\tmulti_record_positions\n'
  cat "${METRICS_DIR}"/calls_*.tsv 2>/dev/null || true
} >"${METRICS_DIR}/summary_calls.tsv.$$"
mv -f "${METRICS_DIR}/summary_calls.tsv.$$" "${METRICS_DIR}/summary_calls.tsv"
log "calls: ${METRICS_DIR}/summary_calls.tsv"
