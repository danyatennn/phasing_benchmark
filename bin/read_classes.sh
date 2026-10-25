#!/usr/bin/env bash
# Classes of untagged reads, confidence of the tagged ones, and read features
# per class, genome-wide and per stratum: plots/read_classes.py for every run
# (caller x phaser) of every cell.
#
#   metrics_reads/<run>.classes.tsv
#   metrics_reads/summary_classes.tsv      all runs, one header
#
# Needs, per cell: the read table (bin/read_tables.sh), the allele matrix of the
# input (bin/allele_matrix.sh), the truth haplotag list (eval_run.sh) and the
# strata cut to the autosomes (bin/prepare_strata.sh with STRATA_AUTOSOMES=1).
# A table newer than all its inputs is kept.
#
#   PHASING_DEPTHS=30 bash bin/read_classes.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

CLASSES_JOBS="${CLASSES_JOBS:-4}"
mkdir -p "${METRICS_READS_DIR}" "${TMP_DIR}" "${LOGS_DIR}/reads"

# the strata of the figures (tier 1), cut to the autosomes; those with nothing
# on the autosomes have no BED and are left out by read_classes.py
STRATA_ARGS=""
for s in $(STRATA_TIER=1 strata_names); do
  STRATA_ARGS="${STRATA_ARGS} --stratum ${s}=$(stratum_bed "${s}")"
done
[[ -n "${STRATA_ARGS}" ]] || warn "no strata in ${STRATA_LIST}: genome-wide rows only"

classes_one() {
  set -euo pipefail
  local caller="$1" platform="$2" depth="$3"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local table="${READS_DIR}/${cell}.reads.tsv.gz"
  local alleles="${READS_DIR}/${cell}_${caller}"
  local truth_tags="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truthtags).tsv"
  for f in "${table}" "${alleles}.reads.tsv.gz" "${truth_tags}"; do
    [[ -s "${f}" ]] || { warn "${cell} ${caller}: missing ${f}"; return 0; }
  done

  local jobs; jobs="$(mktemp -p "${TMP_DIR}")"
  local phaser name tags vcf out
  for phaser in ${PHASER_LIST}; do
    name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
    tags="${METRICS_DIR}/${name}.haplotag.tsv"
    vcf="${PHASE_DIR}/${name}.vcf.gz"
    out="${METRICS_READS_DIR}/${name}.classes.tsv"
    [[ -s "${tags}" && -s "${vcf}" ]] || continue
    if newer_than "${out}" "${tags}" && newer_than "${out}" "${vcf}" && newer_than "${out}" "${table}" \
        && newer_than "${out}" "${alleles}.reads.tsv.gz" && newer_than "${out}" "${alleles}.votes.${phaser}.tsv.gz"; then
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
  log "${cell} ${caller}: read classes for $(wc -l <"${jobs}") runs"
  # shellcheck disable=SC2086
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/read_classes.py" \
      --reads-table "${table}" --alleles "${alleles}" --truth-tags "${truth_tags}" \
      --truth-vcf "${TRUTH_PHASED}" --confident-bed "${CONFIDENT_BED}" \
      ${STRATA_ARGS} --paralog segdups --paralog lowmappabilityall \
      --min-snps "${READACC_MIN_SNPS}" --jobs "${jobs}" \
    2>>"${LOGS_DIR}/reads/${cell}_${caller}.classes.log" \
    || { rm -f "${jobs}"; die "read classes failed for ${cell} ${caller}, see ${LOGS_DIR}/reads/${cell}_${caller}.classes.log"; }
  rm -f "${jobs}"
}

PHASER_LIST="${PHASERS[*]}"
export -f classes_one bam_key vcf_name newer_than log warn die
export METRICS_DIR METRICS_READS_DIR PHASE_DIR READS_DIR LOGS_DIR TMP_DIR REFERENCE_NAME \
       REGION_NAME SAMPLE TRUTH_PHASED CONFIDENT_BED STRATA_ARGS READACC_MIN_SNPS \
       PLOTS_PYTHON PROJECT_ROOT PHASER_LIST

log "${SAMPLE} ${REFERENCE_NAME}: read classes for ${CALLERS[*]} x ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${CLASSES_JOBS}" classes_one {1} {2} {3} \
    ::: "${CALLERS[@]}" ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true

summary="${METRICS_READS_DIR}/summary_classes.tsv"
{
  files=("${METRICS_READS_DIR}"/*.classes.tsv)
  [[ -e "${files[0]}" ]] && head -1 "${files[0]}"
  for f in "${METRICS_READS_DIR}"/*.classes.tsv; do
    [[ -e "${f}" ]] && tail -n +2 "${f}"
  done
  true
} >"${summary}.$$"
mv -f "${summary}.$$" "${summary}"
log "read classes: ${summary}"
