#!/usr/bin/env bash
# Parental k-mer databases for the trio assemblies, built once from the
# parents' short reads (PAT_READS, MAT_READS in config.sh):
#
#   data/trio/pat.yak, mat.yak                         hifiasm -1 / -2
#   data/trio/pat.k30.only.meryl, mat.k30.only.meryl   verkko --hap-kmers ... trio
#
# Only reads are needed, no VCF and no parental assembly.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

[[ -n "${PAT_READS}" && -n "${MAT_READS}" ]] \
  || die "set PAT_READS and MAT_READS to the parents' FASTQ files (config.sh or the environment)"
for f in ${PAT_READS} ${MAT_READS}; do
  require_file "${f}" "parental reads"
done

export PATH="${TRIO_BIN}:${PATH}"
require_cmds yak meryl
MERQURY="${MERQURY:-$(ls -d "${TRIO_BIN}"/../share/merqury* 2>/dev/null | head -1)}"
[[ -s "${MERQURY}/trio/hapmers.sh" ]] \
  || die "merqury not found in ${TRIO_BIN}/../share - set MERQURY to its directory"
export MERQURY

mkdir -p "${TRIO_DIR}"
cd "${TRIO_DIR}"

# yak counts in two passes, so it reads the same input twice
yak_db() {   # <name> <reads...>
  local name="$1"; shift
  [[ -s "${name}.yak" ]] && { log "skip ${name}.yak"; return 0; }
  log "${name}.yak from $# files"
  yak count -k31 -b37 -t "${TRIO_THREADS}" -o "${name}.yak.tmp" <(cat "$@") <(cat "$@")
  mv -f "${name}.yak.tmp" "${name}.yak"
}

# verkko wants the k-mers counted in homopolymer-compressed space
meryl_db() {   # <name> <reads...>
  local name="$1"; shift
  [[ -d "${name}.k30.meryl" ]] && { log "skip ${name}.k30.meryl"; return 0; }
  log "${name}.k30.meryl from $# files"
  rm -rf "${name}.k30.tmp.meryl"
  meryl count compress k=30 threads="${TRIO_THREADS}" memory="${TRIO_MEMORY}" \
      "$@" output "${name}.k30.tmp.meryl"
  mv "${name}.k30.tmp.meryl" "${name}.k30.meryl"
}

yak_db pat ${PAT_READS}
yak_db mat ${MAT_READS}
meryl_db pat ${PAT_READS}
meryl_db mat ${MAT_READS}

# k-mers found in one parent only; without a child database hapmers.sh names
# them <input>.only.meryl
if [[ -d pat.k30.only.meryl && -d mat.k30.only.meryl ]]; then
  log "skip hapmers"
else
  log "hapmers"
  bash "${MERQURY}/trio/hapmers.sh" mat.k30.meryl pat.k30.meryl >hapmers.log 2>&1 \
    || die "hapmers.sh failed, see ${TRIO_DIR}/hapmers.log"
  [[ -d pat.k30.only.meryl && -d mat.k30.only.meryl ]] \
    || die "hapmers.sh wrote no *.only.meryl, see ${TRIO_DIR}/hapmers.log"
fi

log "trio k-mers ready in ${TRIO_DIR}"
