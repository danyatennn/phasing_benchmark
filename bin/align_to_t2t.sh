#!/usr/bin/env bash
# The T2T sources: every source BAM of the sample (its GRCh38 one) turned back
# into reads and aligned to CHM13 with minimap2 (lr:hq for ONT, map-hifi for
# PacBio), into data/t2t_source/<platform>.bam. Only the alignment: the depth
# and the downsampling are the same steps as on GRCh38 (bin/measure_sources.sh,
# downsample.sh), both run by bin/run_reference.sh. A source that exists is kept.
#
#   bash bin/align_to_t2t.sh
#   PHASING_PLATFORMS=ont4k bash bin/align_to_t2t.sh

set -euo pipefail

REFERENCE_NAME="T2T"

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

require_vars SAMPLE SAMTOOLS MINIMAP2 REF_CHM13
require_cmds "${SAMTOOLS}" "${MINIMAP2}"
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

log "${SAMPLE}: aligning ${#PLATFORMS[@]} sources to $(basename "${REF_CHM13}")"
for platform in "${PLATFORMS[@]}"; do
  align_source "${platform}"
done
log "T2T sources: ${T2T_SRC_DIR}"
