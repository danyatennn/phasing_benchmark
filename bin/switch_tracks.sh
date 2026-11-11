#!/usr/bin/env bash
# Where the phase switches, per run (plots/switch_tracks.py), for the switch
# figures of plots/figures_switch.py:
#
#   metrics_switchtracks/<run>.switch_windows.tsv    switch rate per 1 Mb window
#   metrics_switchtracks/<run>.switch_positions.tsv  every switch
#   metrics_switchtracks/<run>.blocks.tsv            every phase block
#   metrics_switchtracks/<run>.read_depth.tsv.gz     concordant / discordant read depth, 10 kb bins
#
# Needs per cell: the phased VCFs and haplotag lists of its runs and the truth
# haplotag list (eval_run.sh); for the read depth, the read table
# (bin/read_tables.sh). Kept when newer than all of them.
#
#   PHASING_DEPTHS=30 PHASING_CALLERS=truth bash bin/switch_tracks.sh
#   PHASING_DEPTHS=30 PHASING_CALLERS=truth REFERENCE_NAME=T2T bash bin/switch_tracks.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

SWITCHTRACKS_DIR="${SAMPLE_ROOT}/metrics_switchtracks"
TRACKS_JOBS="${TRACKS_JOBS:-4}"
mkdir -p "${SWITCHTRACKS_DIR}" "${TMP_DIR}" "${LOGS_DIR}/switchtracks"
require_file "${CHROM_LENGTHS}" "chromosome lengths (eval_run.sh writes them)"

tracks_one() {
  set -euo pipefail
  local caller="$1" platform="$2" depth="$3"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local table="${READS_DIR}/${cell}.reads.tsv.gz"
  local truth_tags="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truthtags).tsv"
  [[ -s "${truth_tags}" ]] || { warn "${cell}: no truth haplotag list"; return 0; }
  [[ -s "${table}" ]] || warn "${cell}: no read table - no read depth"

  local jobs; jobs="$(mktemp -p "${TMP_DIR}")"
  local phaser name tags vcf out
  for phaser in ${PHASER_LIST}; do
    phaser_applies "${phaser}" "${platform}" || continue
    name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
    tags="${METRICS_DIR}/${name}.haplotag.tsv"
    vcf="${PHASE_DIR}/${name}.vcf.gz"
    out="${SWITCHTRACKS_DIR}/${name}.switch_windows.tsv"
    [[ -s "${tags}" && -s "${vcf}" && ! -d "${METRICS_DIR}/${name}.row.tsv.lock" ]] || continue
    if newer_than "${out}" "${vcf}" && newer_than "${out}" "${tags}" && newer_than "${out}" "${truth_tags}" \
        && { [[ ! -s "${table}" ]] || newer_than "${SWITCHTRACKS_DIR}/${name}.read_depth.tsv.gz" "${table}"; }; then
      continue
    fi
    printf '%s\t%s\t%s\n' "${name}" "${tags}" "${vcf}" >>"${jobs}"
  done
  if [[ ! -s "${jobs}" ]]; then
    rm -f "${jobs}"
    log "skip ${cell} ${caller}"
    return 0
  fi
  log "${cell} ${caller}: switch tracks for $(wc -l <"${jobs}") runs"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/switch_tracks.py" \
      --truth-vcf "${TRUTH_PHASED}" --truth-tags "${truth_tags}" \
      ${table:+--reads-table "${table}"} --chrom-lengths "${CHROM_LENGTHS}" \
      --min-snps "${READACC_MIN_SNPS}" --jobs "${jobs}" --out-dir "${SWITCHTRACKS_DIR}" \
    2>>"${LOGS_DIR}/switchtracks/${cell}_${caller}.log" \
    || { rm -f "${jobs}"; die "switch tracks failed for ${cell} ${caller}, see ${LOGS_DIR}/switchtracks/${cell}_${caller}.log"; }
  rm -f "${jobs}"
}

PHASER_LIST="${PHASERS[*]}"
export -f tracks_one bam_key vcf_name newer_than phaser_applies log warn die
export METRICS_DIR PHASE_DIR READS_DIR LOGS_DIR TMP_DIR SWITCHTRACKS_DIR REFERENCE_NAME \
       REGION_NAME SAMPLE TRUTH_PHASED CHROM_LENGTHS READACC_MIN_SNPS PLOTS_PYTHON \
       PROJECT_ROOT PHASER_LIST

log "${SAMPLE} ${REFERENCE_NAME}: switch tracks for ${CALLERS[*]} x ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${TRACKS_JOBS}" tracks_one {1} {2} {3} \
    ::: "${CALLERS[@]}" ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true
log "switch tracks: ${SWITCHTRACKS_DIR}"
