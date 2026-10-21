#!/usr/bin/env bash
# Mean depth of each source BAM of the sample
# data/source_depth[_<ref>].tsv
# You can run it as, for example:
#   bash bin/measure_sources.sh                       HG002, GRCh38
#   REFERENCE_NAME=T2T bash bin/measure_sources.sh
#   SAMPLE=HG005 bash bin/measure_sources.sh
# Only the platforms missing from an existing table are measured
# FORCE_DEPTH=1 measures them all again.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

require_cmds "${MOSDEPTH}"
[[ ${#PLATFORMS[@]} -gt 0 ]] || die "${SAMPLE}: no bam_<platform> in ${SAMPLES_TSV}"

todo=()
for platform in "${PLATFORMS[@]}"; do
  if [[ -n "${FORCE_DEPTH:-}" || ! -s "${DEPTH_TSV}" || -z "$(source_depth "${platform}")" ]]; then
    todo+=("${platform}")
  fi
done
if [[ ${#todo[@]} -eq 0 ]]; then
  log "${DEPTH_TSV} has every platform:"
  cat "${DEPTH_TSV}" >&2
  exit 0
fi

mkdir -p "${DATA_DIR}" "${TMP_DIR}"
tmp="$(mktemp -d -p "${TMP_DIR}")"
trap 'rm -rf "${tmp}"' EXIT

pids=()
for platform in "${todo[@]}"; do
  bam="$(source_bam "${platform}")"
  if [[ "${REFERENCE_NAME}" != GRCh38 && ! -s "${bam}" ]]; then
    die "${SAMPLE} ${platform}: no ${REFERENCE_NAME} source ${bam} - run bin/align_to_t2t.sh first"
  fi
  require_file "${bam}" "source BAM for ${SAMPLE} ${platform}"
  require_file "${bam}.bai" "index of the source BAM (downsample.sh needs <bam>.bai)"
  log "${SAMPLE} ${platform}: measuring $(basename "${bam}")"
  ( MOSDEPTH_THREADS="${DEPTH_THREADS}"; measure_depth "${bam}" >"${tmp}/${platform}" ) &
  pids+=($!)
done
for pid in "${pids[@]}"; do
  wait "${pid}" || die "a depth measurement failed"
done

# the lines of the platforms not measured are kept
{
  printf 'platform\tmean_depth\n'
  if [[ -s "${DEPTH_TSV}" ]]; then
    awk -F'\t' -v skip=" ${todo[*]} " 'NR > 1 && index(skip, " " $1 " ") == 0' "${DEPTH_TSV}"
  fi
  for platform in "${todo[@]}"; do
    printf '%s\t%s\n' "${platform}" "$(<"${tmp}/${platform}")"
  done
} >"${DEPTH_TSV}.tmp"
mv -f "${DEPTH_TSV}.tmp" "${DEPTH_TSV}"
log "${DEPTH_TSV}:"
cat "${DEPTH_TSV}" >&2
