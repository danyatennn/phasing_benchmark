#!/usr/bin/env bash

set -euo pipefail

REFERENCE_NAME="T2T"

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

require_vars SAMPLE SAMTOOLS MOSDEPTH PARALLEL MINIMAP2 REF_CHM13
require_cmds "${SAMTOOLS}" "${MOSDEPTH}" "${PARALLEL}" "${MINIMAP2}"
require_file "${REF_CHM13}" "T2T reference"

mkdir -p "${T2T_SRC_DIR}" "${MMI_DIR}"

minimap2_preset() {
  case "$1" in
    ont4k|ont5k) printf 'lr:hq\n' ;;
    hifi|revio)  printf 'map-hifi\n' ;;
    *) die "no minimap2 preset for platform '$1'" ;;
  esac
}

align_source() {
  local platform="$1"
  local src; src="$(source_bam "${platform}" GRCh38)"
  local out="${T2T_SRC_DIR}/${platform}.bam"
  local preset; preset="$(minimap2_preset "${platform}")"
  local mmi="${MMI_DIR}/chm13.${preset//:/_}.mmi"

  [[ -s "${out}" && -s "${out}.bai" ]] \
    && { log "${platform}: T2T source exists"; return 0; }

  require_file "${src}" "GRCh38 source BAM for ${platform}"

  if [[ ! -s "${mmi}" ]]; then
    log "building minimap2 index for ${preset}"
    "${MINIMAP2}" -t "${ALIGN_THREADS}" -x "${preset}" -d "${mmi}.tmp" "${REF_CHM13}"
    mv -f "${mmi}.tmp" "${mmi}"
  fi

  local rg
  rg="$("${SAMTOOLS}" view -H "${src}" | grep -m1 '^@RG' || true)"
  [[ -n "${rg}" ]] || rg="$(printf '@RG\tID:%s\tSM:%s' "${platform}" "${SAMPLE}")"
  rg="$(printf '%s' "${rg}" | sed 's/\t/\\t/g')"

  log "${platform}: $(basename "${src}") -> T2T (${preset})"
  rm -f "${out}.tmp" "${out}.tmp.bai" "${out}".tmp.tmp.*.bam

  "${SAMTOOLS}" fastq -@ "${SORT_THREADS}" -n -T MM,ML "${src}" \
    | "${MINIMAP2}" -t "${ALIGN_THREADS}" -ax "${preset}" -y -R "${rg}" "${mmi}" - \
    | "${SAMTOOLS}" sort -@ "${SORT_THREADS}" -m "${SORT_MEM}" -o "${out}.tmp" -

  "${SAMTOOLS}" index -@ "${SORT_THREADS}" "${out}.tmp"
  mv -f "${out}.tmp"     "${out}"
  mv -f "${out}.tmp.bai" "${out}.bai"
}

measure_sources() {
  local platform complete=1
  for platform in "${PLATFORMS[@]}"; do
    [[ -s "${DEPTH_TSV}" && -n "$(source_depth "${platform}")" ]] || complete=0
  done
  [[ "${complete}" == 1 ]] \
    && { log "T2T source depth: cached in ${DEPTH_TSV}"; return 0; }

  local tmp="${DATA_DIR}/.source_depth_T2T_tmp"
  rm -rf "${tmp}"; mkdir -p "${tmp}"

  log "measuring T2T source depth: ${PLATFORMS[*]}"
  export -f measure_depth die
  export MOSDEPTH MOSDEPTH_THREADS="${DEPTH_THREADS}"
  "${PARALLEL}" -j "${#PLATFORMS[@]}" \
      "measure_depth ${T2T_SRC_DIR}/{}.bam > ${tmp}/{}" ::: "${PLATFORMS[@]}" \
    || die "T2T source depth measurement failed"

  printf 'platform\tmean_depth\n' >"${DEPTH_TSV}"
  for platform in "${PLATFORMS[@]}"; do
    printf '%s\t%s\n' "${platform}" "$(<"${tmp}/${platform}")" >>"${DEPTH_TSV}"
    log "  ${platform}: $(<"${tmp}/${platform}")x"
  done
  rm -rf "${tmp}"
}

log "aligning ${#PLATFORMS[@]} sources to $(basename "${REF_CHM13}")"
for platform in "${PLATFORMS[@]}"; do
  align_source "${platform}"
done

measure_sources

log "handing the grid to downsample.sh (REFERENCE_NAME=${REFERENCE_NAME})"
REFERENCE_NAME="${REFERENCE_NAME}" bash "${PROJECT_ROOT}/downsample.sh"
