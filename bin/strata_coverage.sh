#!/usr/bin/env bash
# How deep and how ambiguous the alignments are in each stratum: mosdepth over
# the strata of the figures (cut to chr1-22), once with every read and once
# with MAPQ >= 20, the threshold below which whatshap tags no read.
#
#   metrics_strata_auto/coverage_<cell>.tsv   stratum bp depth_mapq0 depth_mapq20
#
# "confident" is the whole confident BED on chr1-22, the genome-wide reference
# the strata are compared with. Kept when newer than the BAM.
#
#   PHASING_DEPTHS=30 bash bin/strata_coverage.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

COVERAGE_JOBS="${COVERAGE_JOBS:-4}"
COVERAGE_THREADS="${COVERAGE_THREADS:-4}"
mkdir -p "${METRICS_STRATA_DIR}" "${TMP_DIR}"
require_cmds "${MOSDEPTH}"

# one BED with every stratum, named in column 4, sorted
REGIONS_BED="${METRICS_STRATA_DIR}/coverage_regions_${REFERENCE_NAME}.bed"
{
  awk -v OFS='\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ { print $1, $2, $3, "confident" }' "${CONFIDENT_BED}"
  for s in $(STRATA_TIER=1 strata_names); do
    bed="$(stratum_bed "${s}")"
    [[ -s "${bed}" ]] && awk -v OFS='\t' -v n="${s}" '{ print $1, $2, $3, n }' "${bed}"
  done
  true
} | sort -k1,1 -k2,2n -k3,3n >"${REGIONS_BED}.$$"
mv -f "${REGIONS_BED}.$$" "${REGIONS_BED}"

coverage_one() {
  set -euo pipefail
  local platform="$1" depth="$2"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local out="${METRICS_STRATA_DIR}/coverage_${cell}.tsv"
  [[ -s "${bam}" ]] || return 0
  newer_than "${out}" "${bam}" && newer_than "${out}" "${REGIONS_BED}" && return 0
  local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
  local q
  for q in 0 20; do
    "${MOSDEPTH}" -n -t "${COVERAGE_THREADS}" -Q "${q}" -b "${REGIONS_BED}" "${tmp}/q${q}" "${bam}"
  done
  # mean over every interval of a stratum, weighted by length
  {
    printf 'stratum\tbp\tdepth_mapq0\tdepth_mapq20\n'
    paste <(gzip -dc "${tmp}/q0.regions.bed.gz") <(gzip -dc "${tmp}/q20.regions.bed.gz") \
      | awk -F'\t' -v OFS='\t' '
          $1 != $6 || $2 != $7 || $3 != $8 || $4 != $9 { print "regions out of step: " $0 > "/dev/stderr"; bad = 1; exit 1 }
          { len = $3 - $2; bp[$4] += len; d0[$4] += len * $5; d20[$4] += len * $10 }
          END { if (bad) exit 1; for (s in bp) printf "%s\t%d\t%.4f\t%.4f\n", s, bp[s], d0[s] / bp[s], d20[s] / bp[s] }' \
      | sort
  } >"${out}.$$"
  mv -f "${out}.$$" "${out}"
  rm -rf "${tmp}"
  log "${cell}: coverage in $(($(wc -l <"${out}") - 1)) regions"
}

export -f coverage_one bam_key sub_bam newer_than log warn die
export METRICS_STRATA_DIR TMP_DIR DATA_DIR REFERENCE_NAME SAMPLE MOSDEPTH COVERAGE_THREADS REGIONS_BED

log "${SAMPLE} ${REFERENCE_NAME}: stratum coverage for ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${COVERAGE_JOBS}" coverage_one {1} {2} ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true
log "coverage: ${METRICS_STRATA_DIR}/coverage_*.tsv"
