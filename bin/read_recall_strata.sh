#!/usr/bin/env bash
# Read-level haplotag recall per genomic region (plots/read_recall_strata.py),
# for every run of every cell: reads on their truth haplotype / reads the truth
# tags, genome-wide and per tier-1 stratum.
#
#   metrics_reads/<run>.recall_strata.tsv
#   metrics_reads/summary_recall_strata.tsv     all runs, one header
#
# Needs per cell: the read table (bin/read_tables.sh), the haplotag lists and
# phased VCFs of its runs and the truth haplotag list (eval_run.sh), the strata
# cut to chr1-22 (bin/prepare_strata.sh). A table newer than all its inputs is
# kept.
#
#   PHASING_DEPTHS=30 bash bin/read_recall_strata.sh
#   PHASING_DEPTHS=30 REFERENCE_NAME=T2T bash bin/read_recall_strata.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

RECALL_JOBS="${RECALL_JOBS:-4}"
mkdir -p "${METRICS_READS_DIR}" "${TMP_DIR}" "${LOGS_DIR}/reads"

STRATA_ARGS=""
for s in $(STRATA_TIER=1 strata_names); do
  STRATA_ARGS="${STRATA_ARGS} --stratum ${s}=$(stratum_bed "${s}")"
done
[[ -n "${STRATA_ARGS}" ]] || warn "no strata in ${STRATA_LIST}: genome-wide rows only"

recall_one() {
  set -euo pipefail
  local caller="$1" platform="$2" depth="$3"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local table="${READS_DIR}/${cell}.reads.tsv.gz"
  local truth_tags="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truthtags).tsv"
  for f in "${table}" "${truth_tags}"; do
    [[ -s "${f}" ]] || { warn "${cell}: missing ${f}"; return 0; }
  done

  local jobs; jobs="$(mktemp -p "${TMP_DIR}")"
  local phaser name tags vcf out
  for phaser in ${PHASER_LIST}; do
    phaser_applies "${phaser}" "${platform}" || continue
    name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
    tags="${METRICS_DIR}/${name}.haplotag.tsv"
    vcf="${PHASE_DIR}/${name}.vcf.gz"
    out="${METRICS_READS_DIR}/${name}.recall_strata.tsv"
    [[ -s "${tags}" && -s "${vcf}" && ! -d "${METRICS_DIR}/${name}.row.tsv.lock" ]] || continue
    if newer_than "${out}" "${tags}" && newer_than "${out}" "${vcf}" && newer_than "${out}" "${table}" \
        && newer_than "${out}" "${truth_tags}"; then
      continue
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${REFERENCE_NAME}" "${SAMPLE}" \
        "${caller}" "${platform}" "${depth}" "${phaser}" "${tags}" "${vcf}" "${out}" >>"${jobs}"
  done
  if [[ ! -s "${jobs}" ]]; then
    rm -f "${jobs}"
    log "skip ${cell} ${caller}"
    return 0
  fi
  log "${cell} ${caller}: read recall per region for $(wc -l <"${jobs}") runs"
  # shellcheck disable=SC2086
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/read_recall_strata.py" \
      --reads-table "${table}" --truth-tags "${truth_tags}" --truth-vcf "${TRUTH_PHASED}" \
      ${STRATA_ARGS} --min-snps "${READACC_MIN_SNPS}" --jobs "${jobs}" \
    2>>"${LOGS_DIR}/reads/${cell}_${caller}.recall_strata.log" \
    || { rm -f "${jobs}"; die "read recall failed for ${cell} ${caller}, see ${LOGS_DIR}/reads/${cell}_${caller}.recall_strata.log"; }
  rm -f "${jobs}"
}

PHASER_LIST="${PHASERS[*]}"
export -f recall_one bam_key vcf_name newer_than phaser_applies log warn die
export METRICS_DIR METRICS_READS_DIR PHASE_DIR READS_DIR LOGS_DIR TMP_DIR REFERENCE_NAME \
       REGION_NAME SAMPLE TRUTH_PHASED STRATA_ARGS READACC_MIN_SNPS PLOTS_PYTHON PROJECT_ROOT PHASER_LIST

log "${SAMPLE} ${REFERENCE_NAME}: read recall per region for ${CALLERS[*]} x ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${RECALL_JOBS}" recall_one {1} {2} {3} \
    ::: "${CALLERS[@]}" ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true

summary="${METRICS_READS_DIR}/summary_recall_strata.tsv"
{
  files=("${METRICS_READS_DIR}"/*.recall_strata.tsv)
  [[ -e "${files[0]}" ]] && head -1 "${files[0]}"
  for f in "${files[@]}"; do
    [[ -e "${f}" ]] && tail -n +2 "${f}"
  done
  true
} >"${summary}.$$"
mv -f "${summary}.$$" "${summary}"
log "read recall per region: ${summary}"
