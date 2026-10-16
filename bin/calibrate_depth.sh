#!/usr/bin/env bash
# What do the phasers' internal coverage caps cost, and do they buy anything?
#
# whatshap caps coverage in its core DP at 15, and that DP is exponential in the
# cap: every +1 roughly doubles the runtime. margin caps at maxDepth 32 in its
# params and -p overrides it. longphase and HapCUT2 have no cap.
#
# Runs on chr20 at 50x with the truth VCF as phaser input, so the caller is out
# of the picture and only the phaser is being measured. Caps go in increasing
# order with a timeout each, so killing this half way still leaves a curve.
#
# metrics/calibration.tsv
#   phaser cap platform depth seconds covered hamming hamming_rate blocks

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

CAL_PLATFORMS=(${CAL_PLATFORMS:-ont4k hifi})
CAL_DEPTH="${CAL_DEPTH:-50}"
CAL_WHATSHAP_CAPS=(${CAL_WHATSHAP_CAPS:-15 16 17 18 19 20})
CAL_MARGIN_DEPTHS=(${CAL_MARGIN_DEPTHS:-32 40 50})
CAL_TIMEOUT="${CAL_TIMEOUT:-6h}"

truth_phased="${DERIVED_DIR}/chr20.phased.vcf.gz"
truth_unphased="${DERIVED_DIR}/chr20.unphased.vcf.gz"
require_file "${truth_phased}" "chr20 phased truth - run bin/run_chr20.sh once"
require_file "${truth_unphased}" "chr20 unphased truth"

work="${TMP_DIR}/calibrate"
out="${METRICS_DIR}/calibration.tsv"
mkdir -p "${work}" "${METRICS_DIR}" "${LOGS_DIR}/calibrate"

[[ -s "${out}" ]] \
  || printf 'phaser\tcap\tplatform\tdepth\tseconds\tcovered\thamming\thamming_rate\tblocks\n' >"${out}"

# chr20 slice of the BAM, so margin does not walk the whole genome and the two
# phasers are timed on the same input
slice_bam() {
  local platform="$1"
  local bam="${work}/${platform}.chr20.bam"
  [[ -s "${bam}.bai" ]] && { printf '%s\n' "${bam}"; return 0; }
  "${SAMTOOLS}" view -@ "${SORT_THREADS}" -b \
      -o "${bam}" "$(sub_bam "${platform}" "${REFERENCE_NAME}" "${CAL_DEPTH}")" chr20
  "${SAMTOOLS}" index -@ "${SORT_THREADS}" "${bam}"
  printf '%s\n' "${bam}"
}

score() {   # <phased vcf> <tag> -> covered hamming rate blocks
  local vcf="$1" tag="$2"
  local cmp="${work}/${tag}.compare.tsv" stats="${work}/${tag}.stats.tsv"
  "${WHATSHAP}" compare --sample "${SAMPLE}" --names truth,test \
      --tsv-pairwise "${cmp}" --only-snvs \
      "${truth_phased}" "${vcf}" >/dev/null 2>&1
  "${WHATSHAP}" stats --tsv "${stats}" "${vcf}" >/dev/null 2>&1
  local covered hamming
  covered="$(cmp_sum "${cmp}" covered_variants)"
  hamming="$(cmp_sum "${cmp}" blockwise_hamming)"
  printf '%s\t%s\t%s\t%s' "${covered}" "${hamming}" \
    "$(awk -v a="${hamming}" -v b="${covered}" 'BEGIN{printf "%.6f", (b ? a/b : 0)}')" \
    "$(stats_col "${stats}" blocks)"
}

run_one() {   # <phaser> <cap> <platform> <bam>
  local phaser="$1" cap="$2" platform="$3" bam="$4"
  local tag="${phaser}_${cap}_${platform}_${CAL_DEPTH}x"
  local logf="${LOGS_DIR}/calibrate/${tag}.log"

  awk -F'\t' -v p="${phaser}" -v c="${cap}" -v pl="${platform}" -v d="${CAL_DEPTH}" \
      'NR > 1 && $1 == p && $2 == c && $3 == pl && $4 == d { found = 1 }
       END { exit !found }' "${out}" \
    && { log "skip ${tag}"; return 0; }

  log "${tag}"
  local t0 t1
  t0="$(date +%s)"

  case "${phaser}" in
    whatshap)
      timeout "${CAL_TIMEOUT}" "${WHATSHAP}" phase \
          --output "${work}/${tag}.vcf.gz" \
          --reference "${REF}" \
          --sample "${SAMPLE}" \
          --ignore-read-groups \
          --internal-downsampling "${cap}" \
          "${truth_unphased}" "${bam}" \
        >"${logf}" 2>&1 \
        || { warn "${tag}: failed or hit the ${CAL_TIMEOUT} timeout, see ${logf}"; return 0; }
      ;;
    margin)
      "${BCFTOOLS}" view "${truth_unphased}" >"${work}/${tag}.in.vcf"
      timeout "${CAL_TIMEOUT}" "${MARGIN}" phase \
          "${bam}" "${REF}" "${work}/${tag}.in.vcf" "$(margin_params "${platform}")" \
          -t "${PHASE_THREADS}" \
          -p "${cap}" \
          -o "${work}/${tag}" \
          -M \
        >"${logf}" 2>&1 \
        || { warn "${tag}: failed or hit the ${CAL_TIMEOUT} timeout, see ${logf}"; return 0; }
      "${BGZIP}" -f -c "${work}/${tag}.phased.vcf" >"${work}/${tag}.vcf.gz"
      ;;
  esac

  t1="$(date +%s)"
  "${TABIX}" -f -p vcf "${work}/${tag}.vcf.gz"

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${phaser}" "${cap}" "${platform}" "${CAL_DEPTH}" "$((t1 - t0))" \
      "$(score "${work}/${tag}.vcf.gz" "${tag}")" \
    >>"${out}"
  log "  $(tail -1 "${out}")"
}

for platform in "${CAL_PLATFORMS[@]}"; do
  bam="$(slice_bam "${platform}")"
  log "${platform}: chr20 slice $(du -h "${bam}" | cut -f1)"
  for cap in "${CAL_WHATSHAP_CAPS[@]}"; do
    run_one whatshap "${cap}" "${platform}" "${bam}"
  done
  for cap in "${CAL_MARGIN_DEPTHS[@]}"; do
    run_one margin "${cap}" "${platform}" "${bam}"
  done
done

log "calibration: ${out}"
column -t "${out}"
