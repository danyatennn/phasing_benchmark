#!/usr/bin/env bash
# Phasing metrics on chr1-22 only for every finished run (plots/autosomes.py),
# from the files eval_run.sh left: nothing is phased or haplotagged again. The
# figures read these, so that every number in the paper is over the autosomes.
#
#   metrics/<run>.auto.tsv        one line per run
#   metrics/summary_auto.tsv      all of them, with a header
#
# A line newer than the run's row, VCF and haplotag list is kept.
#
#   bash bin/autosomes.sh
#   REFERENCE_NAME=T2T bash bin/autosomes.sh
#   PHASING_CALLERS="hifiasm hifiasmtrio verkkotrio" PHASING_PHASERS=dipcall PHASING_DEPTHS="20 30 50" bash bin/autosomes.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

AUTO_THREADS="${AUTO_THREADS:-8}"
mkdir -p "${TMP_DIR}" "${LOGS_DIR}/evaluate"

jobs="$(mktemp -p "${TMP_DIR}")"
trap 'rm -f "${jobs}"' EXIT
for caller in "${CALLERS[@]}"; do
  for platform in "${PLATFORMS[@]}"; do
    for depth in "${DEPTHS[@]}"; do
      for phaser in "${PHASERS[@]}"; do
        name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
        row="${METRICS_DIR}/${name}.row.tsv"
        cmp="${METRICS_DIR}/${name}.compare.tsv"
        stats="${METRICS_DIR}/${name}.stats.tsv"
        tags="${METRICS_DIR}/${name}.haplotag.tsv"
        vcf="${PHASE_DIR}/${name}.vcf.gz"
        out="${METRICS_DIR}/${name}.auto.tsv"
        [[ -s "${row}" && ! -d "${row}.lock" ]] && row_is_valid "${row}" || continue
        [[ -s "${cmp}" && -s "${stats}" && -s "${tags}" && -s "${vcf}" ]] || { warn "${name}: files missing"; continue; }
        newer_than "${out}" "${row}" && newer_than "${out}" "${vcf}" && newer_than "${out}" "${tags}" && continue
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${REFERENCE_NAME}" "${REGION_NAME}" \
            "${caller}" "${platform}" "${depth}" "${phaser}" "${cmp}" "${stats}" "${vcf}" "${tags}" "${out}" >>"${jobs}"
      done
    done
  done
done

if [[ -s "${jobs}" ]]; then
  log "${SAMPLE} ${REFERENCE_NAME}: autosome metrics for $(wc -l <"${jobs}") runs"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/autosomes.py" \
      --truth-vcf "${TRUTH_PHASED}" --min-snps "${READACC_MIN_SNPS}" \
      --jobs "${jobs}" --threads "${AUTO_THREADS}" \
    2>>"${LOGS_DIR}/evaluate/autosomes_${REFERENCE_NAME}.log" \
    || die "autosome metrics failed, see ${LOGS_DIR}/evaluate/autosomes_${REFERENCE_NAME}.log"
else
  log "${SAMPLE} ${REFERENCE_NAME}: every run is up to date"
fi

{
  printf 'reference\tregion\tcaller\tplatform\tdepth\tphaser\tvariants\tphased\tblocks'
  printf '\tswitches\tassessed\thamming\tcovered\ttruth_total\temitted\tshared\tcorrect'
  printf '\tswitch_rate\thamming_rate\trecall_phased\tfnr\tfdr\treads_total\treads_tagged'
  printf '\ttagging_completeness\treads_tagged_min\ttagging_completeness_min\n'
  cat "${METRICS_DIR}"/*.auto.tsv 2>/dev/null || true
} >"${METRICS_DIR}/summary_auto.tsv.$$"
mv -f "${METRICS_DIR}/summary_auto.tsv.$$" "${METRICS_DIR}/summary_auto.tsv"
log "autosomes: ${METRICS_DIR}/summary_auto.tsv"
