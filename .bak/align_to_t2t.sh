#!/usr/bin/env bash

# Realign the GRCh38 grid onto T2T/CHM13, keeping the read sets identical.
#
# Two phases:
#   1. align each source BAM to CHM13 once (4 alignments, not 20)
#   2. subsample those with the fractions the GRCh38 run already recorded
#
# Fractions are computed against the depth of the T2T-aligned source, so every
# cell lands on its nominal depth on T2T. That is the point of doing it this
# way: the depth matches the group, and the read set is consequently NOT the
# same as the matching GRCh38 cell (samtools -s picks reads by name hash, and a
# different fraction means a different hash threshold).
#
# MM/ML stay on the reads; HP/PS/PC do not survive the FASTQ round trip, so the
# T2T BAMs need no strip_tags.sh pass.
#
# Both phases skip finished outputs, so the script can be re-run.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

require_cmds "${SAMTOOLS}" "${MINIMAP2}"
require_file "${REF_CHM13}" "T2T reference"
require_cmds "${MOSDEPTH}" "${PARALLEL}"
mkdir -p "${T2T_SRC_DIR}" "${MMI_DIR}"

minimap2_preset() {
  case "$1" in
    ont4k|ont5k) printf 'lr:hq\n' ;;
    hifi|revio)  printf 'map-hifi\n' ;;
    *) die "no minimap2 preset for platform '$1'" ;;
  esac
}

preset_mmi() {
  printf '%s/chm13.%s.mmi\n' "${MMI_DIR}" "${1//:/_}"
}

# ------------------------------------------------- phase 1: align the sources

# One index per preset: lr:hq and map-hifi produce different indexes.
build_index() {
  local preset="$1"
  local mmi; mmi="$(preset_mmi "${preset}")"

  [[ -s "${mmi}" ]] && { log "index for ${preset}: exists"; return 0; }

  log "building minimap2 index for ${preset}"
  "${MINIMAP2}" -t "${ALIGN_THREADS}" -x "${preset}" -d "${mmi}.tmp" "${REF_CHM13}"
  mv -f "${mmi}.tmp" "${mmi}"
}

align_source() {
  local platform="$1"
  local src; src="$(source_bam "${platform}")"
  local preset; preset="$(minimap2_preset "${platform}")"
  local mmi; mmi="$(preset_mmi "${preset}")"
  local out="${T2T_SRC_DIR}/${platform}.bam"

  [[ -s "${out}" && -s "${out}.bai" ]] \
    && { log "${platform}: T2T source exists"; return 0; }

  require_file "${src}" "source BAM for ${platform}"

  # Carry the source read group over, so SM survives the FASTQ round trip.
  # minimap2 -R wants literal \t, not real tabs.
  local rg
  rg="$("${SAMTOOLS}" view -H "${src}" | grep -m1 '^@RG' || true)"
  [[ -n "${rg}" ]] \
    || rg="$(printf '@RG\tID:%s\tSM:%s' "${platform}" "${SAMPLE}")"
  rg="$(printf '%s' "${rg}" | sed 's/\t/\\t/g')"

  log "${platform}: ${src} -> T2T (${preset})"
  log "  ${rg}"

  rm -f "${out}.tmp" "${out}.tmp.bai"

  "${SAMTOOLS}" fastq -@ "${SORT_THREADS}" -n -T MM,ML "${src}" \
    | "${MINIMAP2}" -t "${ALIGN_THREADS}" -ax "${preset}" -y -R "${rg}" "${mmi}" - \
    | "${SAMTOOLS}" sort -@ "${SORT_THREADS}" -o "${out}.tmp" -

  "${SAMTOOLS}" index -@ "${SORT_THREADS}" "${out}.tmp"

  mv -f "${out}.tmp"     "${out}"
  mv -f "${out}.tmp.bai" "${out}.bai"

  log "${platform}: T2T source ready"
}

# ------------------------------------------- phase 2: depth of the T2T sources

t2t_source_depth() {
  awk -F'\t' -v p="$1" '$1 == p {print $2}' "${T2T_DEPTH_TSV}"
}

measure_t2t_sources() {
  local missing=() platform
  for platform in "${PLATFORMS[@]}"; do
    [[ -s "${T2T_DEPTH_TSV}" ]] && [[ -n "$(t2t_source_depth "${platform}")" ]] \
      && { log "  ${platform}: $(t2t_source_depth "${platform}")x (cached)"; continue; }
    missing+=("${platform}")
  done
  [[ ${#missing[@]} -eq 0 ]] && return 0

  local tmp="${DATA_DIR}/.t2t_source_depth_tmp"
  rm -rf "${tmp}"; mkdir -p "${tmp}"

  # Keep whatever is already measured, re-measure only what is missing.
  for platform in "${PLATFORMS[@]}"; do
    local cached; cached="$([[ -s "${T2T_DEPTH_TSV}" ]] && t2t_source_depth "${platform}" || true)"
    [[ -n "${cached}" ]] && printf '%s\n' "${cached}" >"${tmp}/${platform}"
  done

  log "measuring T2T source depth for: ${missing[*]}"
  MOSDEPTH_THREADS="${DEPTH_THREADS}" "${PARALLEL}" -j "${#missing[@]}" \
      "measure_depth ${T2T_SRC_DIR}/{}.bam > ${tmp}/{}" ::: "${missing[@]}" \
    || die "T2T source depth measurement failed"

  printf 'platform\tmean_depth\n' >"${T2T_DEPTH_TSV}"
  for platform in "${PLATFORMS[@]}"; do
    printf '%s\t%s\n' "${platform}" "$(<"${tmp}/${platform}")" >>"${T2T_DEPTH_TSV}"
    log "  ${platform}: $(<"${tmp}/${platform}")x"
  done
  rm -rf "${tmp}"
}

# ------------------------------------------------- phase 3: subsample the grid

subsample_cell() {
  set -euo pipefail

  local platform="$1" depth="$2"
  local t2t_src="${T2T_SRC_DIR}/${platform}.bam"
  local key; key="$(bam_key "${platform}" "${T2T_NAME}" "${depth}")"
  local out; out="$(sub_bam "${platform}" "${T2T_NAME}" "${depth}")"

  require_file "${t2t_src}" "T2T source for ${platform}"

  local src_depth; src_depth="$(t2t_source_depth "${platform}")"
  [[ -n "${src_depth}" ]] || die "no T2T source depth for ${platform}"

  mkdir -p "$(dirname "${out}")"
  rm -f "${out}.tmp" "${out}.tmp.bai"

  local frac
  frac="$(awk -v target="${depth}" -v source="${src_depth}" \
      'BEGIN { if (source <= 0) exit 1; printf "%.6f", target / source }' \
    || die "${key}: bad T2T source depth '${src_depth}'")"

  if awk -v f="${frac}" 'BEGIN { exit !(f >= 0.999999) }'; then
    warn "${key}: T2T source is ${src_depth}x; using the full BAM"
    frac="1.000000"
    "${SAMTOOLS}" view -@ "${SORT_THREADS}" -b -o "${out}.tmp" "${t2t_src}"
  else
    local sarg
    sarg="$(awk -v seed="${SEED}" -v frac="${frac}" \
        'BEGIN { printf "%d.%06d", seed, int(frac * 1000000 + 0.5) }')"
    log "${key}: T2T source=${src_depth}x target=${depth}x fraction=${frac}"
    "${SAMTOOLS}" view -@ "${SORT_THREADS}" -b -s "${sarg}" -o "${out}.tmp" "${t2t_src}"
  fi

  "${SAMTOOLS}" index -@ "${SORT_THREADS}" "${out}.tmp"
  mv -f "${out}.tmp"     "${out}"
  mv -f "${out}.tmp.bai" "${out}.bai"

  local got; got="$(MOSDEPTH_THREADS="${DEPTH_THREADS}" measure_depth "${out}")"
  log "${key}: achieved ${got}x"

  printf '%s\t%s\t%s\t%s\n' "${platform}" "${depth}" "${frac}" "${got}" \
    >"${T2T_TMP_DEPTH_DIR}/${platform}_${depth}.tsv"
}

export -f subsample_cell measure_depth t2t_source_depth bam_key sub_bam \
          require_file log warn die
export SAMTOOLS MOSDEPTH DATA_DIR SAMPLE SEED SORT_THREADS DEPTH_THREADS \
       T2T_SRC_DIR T2T_NAME T2T_DEPTH_TSV T2T_TMP_DEPTH_DIR

log "phase 1: aligning ${#PLATFORMS[@]} source BAMs to ${REF_CHM13}"
for preset in lr:hq map-hifi; do
  build_index "${preset}"
done
for platform in "${PLATFORMS[@]}"; do
  align_source "${platform}"
done

log "phase 2: T2T source depth"
measure_t2t_sources

log "phase 3: subsampling to ${DEPTHS[*]}x on T2T"
rm -rf "${T2T_TMP_DEPTH_DIR}"
mkdir -p "${T2T_TMP_DEPTH_DIR}"

"${PARALLEL}" \
    --joblog "${DATA_DIR}/parallel_t2t_subsample.log" \
    -j "${PARALLEL_JOBS}" \
    subsample_cell {1} {2} \
    ::: "${PLATFORMS[@]}" \
    ::: "${DEPTHS[@]}"

printf 'platform\tdepth\tfraction\tachieved_depth\n' >"${T2T_ACHIEVED_TSV}"
find "${T2T_TMP_DEPTH_DIR}" -type f -name '*.tsv' -print0 \
  | sort -z | xargs -0 -r cat >>"${T2T_ACHIEVED_TSV}"
rm -rf "${T2T_TMP_DEPTH_DIR}"

log "T2T grid done"
log "  source depth:   ${T2T_DEPTH_TSV}"
log "  achieved depth: ${T2T_ACHIEVED_TSV}"
