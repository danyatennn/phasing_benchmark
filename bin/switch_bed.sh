#!/usr/bin/env bash
# Where the switch errors are: whatshap compare --switch-error-bed for every run
# of the grid, against the phased truth, SNVs only - the same comparison as
# eval_run.sh, with the positions written out. One BED line per switch: the
# interval between the two SNPs it falls between.
#
#   metrics_switchbed/<run>.switches.bed
#
# Kept when newer than the run's VCF and the truth.
#
#   PHASING_DEPTHS=30 bash bin/switch_bed.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

SWITCH_JOBS="${SWITCH_JOBS:-12}"
mkdir -p "${SWITCHBED_DIR}" "${LOGS_DIR}/switchbed"

switch_one() {
  set -euo pipefail
  local caller="$1" platform="$2" depth="$3" phaser="$4"
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
  local vcf="${PHASE_DIR}/${name}.vcf.gz"
  local out="${SWITCHBED_DIR}/${name}.switches.bed"
  is_done "${vcf}" || return 0
  newer_than "${out}" "${vcf}" && newer_than "${out}" "${TRUTH_PHASED}" && return 0
  "${WHATSHAP}" compare --sample "${SAMPLE}" --names truth,"${phaser}" --only-snvs \
      --switch-error-bed "${out}.$$" "${TRUTH_PHASED}" "${vcf}" \
    >"${LOGS_DIR}/switchbed/${name}.log" 2>&1 \
    || { rm -f "${out}.$$"; die "whatshap compare failed for ${name}, see ${LOGS_DIR}/switchbed/${name}.log"; }
  # whatshap writes no file when there is no switch at all
  [[ -e "${out}.$$" ]] || : >"${out}.$$"
  mv -f "${out}.$$" "${out}"
  log "${name}: $(wc -l <"${out}") switches"
}

export -f switch_one vcf_name is_done newer_than log warn die
export SWITCHBED_DIR PHASE_DIR LOGS_DIR TRUTH_PHASED WHATSHAP SAMPLE REFERENCE_NAME REGION_NAME

log "${SAMPLE} ${REFERENCE_NAME}: switch positions for ${CALLERS[*]} x ${PLATFORMS[*]} x ${DEPTHS[*]} x ${PHASERS[*]}"
"${PARALLEL}" -j "${SWITCH_JOBS}" switch_one {1} {2} {3} {4} \
    ::: "${CALLERS[@]}" ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" ::: "${PHASERS[@]}" || true
log "switch positions: ${SWITCHBED_DIR}"
