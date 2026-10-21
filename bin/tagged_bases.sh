#!/usr/bin/env bash
# Share of read bases put on a haplotype, for every run of the grid
# (plots/tagged_bases.py), from the cell's read table: chr1-22, blocks of at
# least READACC_MIN_SNPS phased SNPs, the same reads for the truth and the
# callers.
#
#   metrics/<run>.tagbases_<rules>.tsv           <rules>: auto_min2 by default
#   metrics/summary_tagbases_<rules>.tsv
#
# Needs bin/read_tables.sh. A line newer than the run's haplotag list, VCF and
# the read table is kept.
#
#   bash bin/tagged_bases.sh
#   REFERENCE_NAME=T2T bash bin/tagged_bases.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

TAGBASES_JOBS="${TAGBASES_JOBS:-6}"
rules="$(READACC_AUTOSOMES=1 readacc_suffix)"
rules="${rules#readacc_}"
TAGBASES_SUFFIX="tagbases_${rules}"
mkdir -p "${TMP_DIR}" "${LOGS_DIR}/reads"

tagbases_cell() {
  set -euo pipefail
  local platform="$1" depth="$2"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local table="${READS_DIR}/${cell}.reads.tsv.gz"
  [[ -s "${table}" ]] || { warn "${cell}: no read table - run bin/read_tables.sh"; return 0; }
  local jobs; jobs="$(mktemp -p "${TMP_DIR}")"
  local caller phaser name tags vcf out
  for caller in ${CALLER_LIST}; do
    for phaser in ${PHASER_LIST}; do
      name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
      tags="${METRICS_DIR}/${name}.haplotag.tsv"
      vcf="${PHASE_DIR}/${name}.vcf.gz"
      out="${METRICS_DIR}/${name}.${TAGBASES_SUFFIX}.tsv"
      [[ -s "${tags}" && -s "${vcf}" && ! -d "${METRICS_DIR}/${name}.row.tsv.lock" ]] || continue
      newer_than "${out}" "${tags}" && newer_than "${out}" "${vcf}" && newer_than "${out}" "${table}" && continue
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${REFERENCE_NAME}" "${caller}" "${platform}" \
          "${depth}" "${phaser}" "${tags}" "${vcf}" "${out}" >>"${jobs}"
    done
  done
  if [[ ! -s "${jobs}" ]]; then
    rm -f "${jobs}"; log "skip ${cell}"; return 0
  fi
  log "${cell}: tagged bases for $(wc -l <"${jobs}") runs"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/tagged_bases.py" \
      --reads-table "${table}" --min-snps "${READACC_MIN_SNPS}" --jobs "${jobs}" \
    2>>"${LOGS_DIR}/reads/${cell}.tagbases.log" \
    || { rm -f "${jobs}"; die "tagged bases failed for ${cell}, see ${LOGS_DIR}/reads/${cell}.tagbases.log"; }
  rm -f "${jobs}"
}

CALLER_LIST="${CALLERS[*]}"
PHASER_LIST="${PHASERS[*]}"
export -f tagbases_cell bam_key vcf_name newer_than log warn die
export METRICS_DIR PHASE_DIR READS_DIR LOGS_DIR TMP_DIR REFERENCE_NAME REGION_NAME SAMPLE \
       READACC_MIN_SNPS PLOTS_PYTHON PROJECT_ROOT CALLER_LIST PHASER_LIST TAGBASES_SUFFIX

log "${SAMPLE} ${REFERENCE_NAME}: tagged bases (${TAGBASES_SUFFIX}) for ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${TAGBASES_JOBS}" tagbases_cell {1} {2} ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true

{
  printf 'reference\tcaller\tplatform\tdepth\tphaser\treads_total\treads_tagged\tbases_total\tbases_tagged\tbases_tagged_fraction\n'
  cat "${METRICS_DIR}"/*."${TAGBASES_SUFFIX}".tsv 2>/dev/null || true
} >"${METRICS_DIR}/summary_${TAGBASES_SUFFIX}.tsv.$$"
mv -f "${METRICS_DIR}/summary_${TAGBASES_SUFFIX}.tsv.$$" "${METRICS_DIR}/summary_${TAGBASES_SUFFIX}.tsv"
log "tagged bases: ${METRICS_DIR}/summary_${TAGBASES_SUFFIX}.tsv"
