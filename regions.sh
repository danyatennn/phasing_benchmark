#!/usr/bin/env bash
# Metrics inside the GIAB stratification regions, from the VCFs that are already
# on disk. No caller and no phaser is re-run: a stratum is just an intersection.
#
# metrics_strata/calls_<cell>_<caller>__<stratum>.row.tsv   caller in the stratum
# metrics_strata/phase_<run>__<stratum>.row.tsv             phasing in the stratum
# metrics_strata/summary_calls_strata.tsv
# metrics_strata/summary_phase_strata.tsv
#
# Needs bin/prepare_strata.sh to have run for this reference.
#
# Splitting the work:
#   STRATA_TIER=1                      only the headline strata
#   STRATA_ONLY="segdups MHC"          named ones
#   PHASING_PLATFORMS=ont4k            the usual grid variables also narrow it
#   REGIONS_SHARD=2/4                  every 4th job, offset 2 - for four hosts
#
# Read-level metrics are whole-genome only: a read's haplotype comes from every
# variant it covers, so it cannot be restricted to a stratum without a separate
# haplotag pass per stratum, which would be thousands of passes over WGS BAMs.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${METRICS_STRATA_DIR}" "${LOGS_DIR}/regions" "${TMP_DIR}"

# how many of the truth sites in a stratum a VCF got right
shared_with_truth() {
  "${BCFTOOLS}" isec -n=2 -w1 -c none "$1" "$2" \
    | awk '!/^#/ { n++ } END { print n + 0 }'
}

calls_one() {
  set -euo pipefail

  local caller="$1" platform="$2" depth="$3" stratum="$4"
  local cell; cell="$(vcf_name "${platform}" "${depth}" "${caller}")"
  local row="${METRICS_STRATA_DIR}/calls_${cell}__${stratum}.row.tsv"

  [[ -s "${row}" ]] && return 0

  local bed; bed="$(stratum_bed "${stratum}")"
  local truth; truth="$(stratum_truth "${stratum}")"
  [[ -s "${truth}.tbi" ]] || return 0

  local vcf="${VCF_DIR}/${cell}.het_snps.vcf.gz"
  [[ -s "${vcf}" ]] || return 0

  local lock="${row}.lock"
  mkdir "${lock}" 2>/dev/null || return 0
  printf '%s %s %s\n' "$(hostname)" "$$" "$(date '+%F %H:%M')" >"${lock}/owner"
  trap "rm -f '${lock}/owner'; rmdir '${lock}' 2>/dev/null" EXIT INT TERM HUP

  local truth_total; truth_total="$("${BCFTOOLS}" index -n "${truth}")"

  local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
  "${BCFTOOLS}" view -T "${bed}" -Oz -o "${tmp}/in.vcf.gz" "${vcf}"
  "${TABIX}" -f -p vcf "${tmp}/in.vcf.gz"

  local emitted shared
  emitted="$("${BCFTOOLS}" index -n "${tmp}/in.vcf.gz")"
  shared=0
  [[ "${emitted}" -gt 0 && "${truth_total}" -gt 0 ]] \
    && shared="$(shared_with_truth "${truth}" "${tmp}/in.vcf.gz")"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${REFERENCE_NAME}" "${REGION_NAME}" "${caller}" "${platform}" "${depth}" "${stratum}" \
      "${truth_total}" "${emitted}" "${shared}" \
      "$(( emitted - shared ))" "$(( truth_total - shared ))" \
      "$(awk -v s="${shared}" -v m="${emitted}" 'BEGIN{printf "%.6f", (m ? s/m : 0)}')" \
      "$(awk -v s="${shared}" -v t="${truth_total}" 'BEGIN{printf "%.6f", (t ? s/t : 0)}')" \
      "$(awk -v s="${shared}" -v m="${emitted}" -v t="${truth_total}" \
           'BEGIN{p = m ? s/m : 0; r = t ? s/t : 0; printf "%.6f", (p + r ? 2*p*r/(p+r) : 0)}')" \
    >"${row}.$$"
  mv -f "${row}.$$" "${row}"
  rm -rf "${tmp}"
}

phase_one_region() {
  set -euo pipefail

  local caller="$1" platform="$2" depth="$3" phaser="$4" stratum="$5"
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
  local row="${METRICS_STRATA_DIR}/phase_${name}__${stratum}.row.tsv"

  [[ -s "${row}" ]] && return 0

  local bed; bed="$(stratum_bed "${stratum}")"
  local truth; truth="$(stratum_truth "${stratum}")"
  [[ -s "${truth}.tbi" ]] || return 0

  local vcf="${PHASE_DIR}/${name}.vcf.gz"
  [[ -s "${vcf}" ]] || return 0

  local lock="${row}.lock"
  mkdir "${lock}" 2>/dev/null || return 0
  printf '%s %s %s\n' "$(hostname)" "$$" "$(date '+%F %H:%M')" >"${lock}/owner"
  trap "rm -f '${lock}/owner'; rmdir '${lock}' 2>/dev/null" EXIT INT TERM HUP

  local truth_total; truth_total="$("${BCFTOOLS}" index -n "${truth}")"

  local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
  "${BCFTOOLS}" view -T "${bed}" -Oz -o "${tmp}/in.vcf.gz" "${vcf}"
  "${TABIX}" -f -p vcf "${tmp}/in.vcf.gz"

  local emitted; emitted="$("${BCFTOOLS}" index -n "${tmp}/in.vcf.gz")"

  # a stratum can be empty inside the confident BED, and whatshap has nothing to
  # say about an empty VCF - record the zeros so the gap is visible in the table
  local variants=0 phased=0 phased_fraction=0 blocks=0 block_n50=0
  local covered=0 assessed=0 switches=0 hamming=0 shared=0
  if [[ "${emitted}" -gt 0 && "${truth_total}" -gt 0 ]]; then
    local logf="${LOGS_DIR}/regions/${name}__${stratum}.log"
    "${WHATSHAP}" stats --tsv "${tmp}/stats.tsv" --chr-lengths "${CHROM_LENGTHS}" \
        "${tmp}/in.vcf.gz" >/dev/null 2>"${logf}"
    "${WHATSHAP}" compare --sample "${SAMPLE}" --names truth,"${phaser}" \
        --tsv-pairwise "${tmp}/cmp.tsv" --only-snvs \
        "${truth}" "${tmp}/in.vcf.gz" >>"${logf}" 2>&1

    variants="$(stats_col "${tmp}/stats.tsv" variants)"
    phased="$(stats_col "${tmp}/stats.tsv" phased)"
    phased_fraction="$(stats_col "${tmp}/stats.tsv" phased_fraction)"
    blocks="$(stats_col "${tmp}/stats.tsv" blocks)"
    block_n50="$(stats_col "${tmp}/stats.tsv" block_n50)"
    covered="$(cmp_sum "${tmp}/cmp.tsv" covered_variants)"
    assessed="$(cmp_sum "${tmp}/cmp.tsv" all_assessed_pairs)"
    switches="$(cmp_sum "${tmp}/cmp.tsv" all_switches)"
    hamming="$(cmp_sum "${tmp}/cmp.tsv" blockwise_hamming)"
    shared="$(shared_with_truth "${truth}" "${tmp}/in.vcf.gz")"
  fi
  local correct=$(( covered - hamming ))

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${REFERENCE_NAME}" "${REGION_NAME}" "${caller}" "${platform}" "${depth}" "${phaser}" "${stratum}" \
      "${variants}" "${phased}" "${phased_fraction}" "${blocks}" "${block_n50}" \
      "${switches}" \
      "$(awk -v a="${switches}" -v b="${assessed}" 'BEGIN{printf "%.6f", (b ? a/b : 0)}')" \
      "${hamming}" \
      "$(awk -v a="${hamming}" -v b="${covered}" 'BEGIN{printf "%.6f", (b ? a/b : 0)}')" \
      "${covered}" "${assessed}" "${truth_total}" "${emitted}" "${shared}" "${correct}" \
      "$(awk -v s="${shared}" -v t="${truth_total}" 'BEGIN{printf "%.6f", (t ? 1 - s/t : 0)}')" \
      "$(awk -v s="${shared}" -v m="${emitted}"    'BEGIN{printf "%.6f", (m ? 1 - s/m : 0)}')" \
      "$(awk -v c="${correct}" -v t="${truth_total}" 'BEGIN{printf "%.6f", (t ? c/t : 0)}')" \
      "$(awk -v c="${correct}" -v m="${emitted}"     'BEGIN{printf "%.6f", (m ? c/m : 0)}')" \
    >"${row}.$$"
  mv -f "${row}.$$" "${row}"
  rm -rf "${tmp}"
}

# every n-th line, offset i, so REGIONS_SHARD=2/4 on four hosts splits evenly
shard() {
  [[ -z "${REGIONS_SHARD}" ]] && { cat; return 0; }
  awk -v i="${REGIONS_SHARD%/*}" -v n="${REGIONS_SHARD#*/}" 'NR % n == i % n'
}

export -f calls_one phase_one_region shared_with_truth \
          vcf_name stratum_bed stratum_truth stats_col cmp_sum log warn die
export METRICS_STRATA_DIR PREPARED_DIR STRATA_TRUTH_DIR VCF_DIR PHASE_DIR LOGS_DIR \
       TMP_DIR REFERENCE_NAME REGION_NAME SAMPLE CHROM_LENGTHS \
       WHATSHAP BCFTOOLS TABIX

strata="$(strata_names)"
n_strata="$(printf '%s\n' "${strata}" | grep -c . || true)"
log "${REFERENCE_NAME}: ${n_strata} strata"

log "caller performance per stratum"
joblog="${LOGS_DIR}/parallel_regions_calls_${RUN_TAG}.log"
keep_joblog "${joblog}"
for caller in "${CALLERS[@]}"; do
  [[ "${caller}" == truth ]] && continue
  for platform in "${PLATFORMS[@]}"; do
    for depth in "${DEPTHS[@]}"; do
      for stratum in ${strata}; do
        printf '%s\t%s\t%s\t%s\n' "${caller}" "${platform}" "${depth}" "${stratum}"
      done
    done
  done
done | shard | "${PARALLEL}" --joblog "${joblog}" --colsep '\t' -j "${REGIONS_JOBS}" \
    calls_one {1} {2} {3} {4} || true

log "phasing per stratum"
joblog="${LOGS_DIR}/parallel_regions_phase_${RUN_TAG}.log"
keep_joblog "${joblog}"
for caller in "${CALLERS[@]}"; do
  for platform in "${PLATFORMS[@]}"; do
    for depth in "${DEPTHS[@]}"; do
      for phaser in "${PHASERS[@]}"; do
        for stratum in ${strata}; do
          printf '%s\t%s\t%s\t%s\t%s\n' "${caller}" "${platform}" "${depth}" "${phaser}" "${stratum}"
        done
      done
    done
  done
done | shard | "${PARALLEL}" --joblog "${joblog}" --colsep '\t' -j "${REGIONS_JOBS}" \
    phase_one_region {1} {2} {3} {4} {5} || true

{
  printf 'reference\tregion\tcaller\tplatform\tdepth\tstratum'
  printf '\ttruth_total\temitted\tshared\tfp\tfn\tprecision\trecall\tf1\n'
  cat "${METRICS_STRATA_DIR}"/calls_*.row.tsv 2>/dev/null
} >"${METRICS_STRATA_DIR}/summary_calls_strata.tsv.$$"
mv -f "${METRICS_STRATA_DIR}/summary_calls_strata.tsv.$$" "${METRICS_STRATA_DIR}/summary_calls_strata.tsv"

{
  printf 'reference\tregion\tcaller\tplatform\tdepth\tphaser\tstratum'
  printf '\tvariants\tphased\tphased_fraction\tblocks\tblock_n50'
  printf '\tswitches\tswitch_rate\thamming\thamming_rate'
  printf '\tcovered\tassessed\ttruth_total\temitted\tshared\tcorrect'
  printf '\tfnr\tfdr\trecall_phased\tprecision_phased\n'
  cat "${METRICS_STRATA_DIR}"/phase_*.row.tsv 2>/dev/null
} >"${METRICS_STRATA_DIR}/summary_phase_strata.tsv.$$"
mv -f "${METRICS_STRATA_DIR}/summary_phase_strata.tsv.$$" "${METRICS_STRATA_DIR}/summary_phase_strata.tsv"

log "metrics_strata: summary_calls_strata.tsv, summary_phase_strata.tsv"
