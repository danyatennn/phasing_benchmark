#!/usr/bin/env bash
# The read table of every cell of the grid: plots/read_table.py on its BAM.
#
#   reads/<platform>_<ref>_<sample>_<depth>x_wgs.reads.tsv.gz
#
# Used by readacc_anchor.sh, bin/read_classes.sh and plots/tagged_bases.py. A
# table newer than its BAM is kept.
#
#   bash bin/read_tables.sh                                   HG002, GRCh38, all depths
#   REFERENCE_NAME=T2T PHASING_DEPTHS=30 bash bin/read_tables.sh
#   SAMPLE=HG005 bash bin/read_tables.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

READS_JOBS="${READS_JOBS:-6}"
READS_THREADS="${READS_THREADS:-4}"
mkdir -p "${READS_DIR}" "${LOGS_DIR}/reads"

table_one() {
  set -euo pipefail
  local platform="$1" depth="$2"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local out="${READS_DIR}/${cell}.reads.tsv.gz"
  [[ -s "${bam}" ]] || { warn "no BAM for ${cell}"; return 0; }
  newer_than "${out}" "${bam}" && { log "skip ${cell}"; return 0; }

  local lock="${out}.lock"
  mkdir "${lock}" 2>/dev/null || { log "busy elsewhere ${cell}"; return 0; }
  trap "rmdir '${lock}' 2>/dev/null" EXIT INT TERM HUP

  log "${cell}: read table"
  "${PHASING_PYTHON}" "${PROJECT_ROOT}/plots/read_table.py" \
      --bam "${bam}" --out "${out}" --threads "${READS_THREADS}" \
    2>"${LOGS_DIR}/reads/${cell}.reads.log" \
    || die "read table failed for ${cell}, see ${LOGS_DIR}/reads/${cell}.reads.log"
}

export -f table_one bam_key sub_bam newer_than log warn die
export READS_DIR LOGS_DIR DATA_DIR REFERENCE_NAME SAMPLE PHASING_PYTHON PROJECT_ROOT READS_THREADS

log "${SAMPLE} ${REFERENCE_NAME}: read tables for ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${READS_JOBS}" table_one {1} {2} ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true
log "read tables: ${READS_DIR}"
