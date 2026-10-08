#!/usr/bin/env bash
# Metrics for every phased VCF. Per run in metrics/:
#   <name>.stats.tsv     whatshap stats   - blocks, N50, phased fraction
#   <name>.compare.tsv   whatshap compare - switches, Hamming (one row per chromosome)
#   <name>.haplotag.tsv  whatshap haplotag - read -> haplotype
#   <cell>_truthtags.tsv whatshap haplotag against the phased truth, per cell
#   <name>.row.tsv       one summary line
# and metrics/summary.tsv with every row.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${METRICS_DIR}" "${LOGS_DIR}/evaluate"

[[ -s "${CHROM_LENGTHS}" ]] || cut -f1,2 "${REF}.fai" >"${CHROM_LENGTHS}"

eval_one() {
  set -euo pipefail

  local caller="$1" platform="$2" depth="$3" phaser="$4"
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
  local vcf="${PHASE_DIR}/${name}.vcf.gz"
  local row="${METRICS_DIR}/${name}.row.tsv"

  [[ -s "${row}" ]] && { log "skip ${name}"; return 0; }
  [[ -s "${vcf}" ]] || { warn "no phased VCF for ${name}"; return 0; }

  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local logf="${LOGS_DIR}/evaluate/${name}.log"
  local stats="${METRICS_DIR}/${name}.stats.tsv"
  local cmp="${METRICS_DIR}/${name}.compare.tsv"
  local tags="${METRICS_DIR}/${name}.haplotag.tsv"
  # depends on the cell only, so it is shared by all callers and phasers
  local truth_tags="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truthtags).tsv"

  log "evaluating ${name}"
  : >"${logf}"

  [[ -s "${stats}" ]] || "${WHATSHAP}" stats \
      --tsv "${stats}" \
      --chr-lengths "${CHROM_LENGTHS}" \
      "${vcf}" \
    >"${METRICS_DIR}/${name}.stats.txt" 2>>"${logf}"

  # truth first, then the phasing under test
  "${WHATSHAP}" compare \
      --sample "${SAMPLE}" \
      --names truth,"${phaser}" \
      --tsv-pairwise "${cmp}" \
      --only-snvs \
      "${TRUTH_PHASED}" "${vcf}" \
    >>"${logf}" 2>&1

  # the slowest step here, and it only depends on the VCF - keep it so adding
  # columns later means deleting row.tsv, not redoing every haplotag
  if [[ ! -s "${tags}" ]]; then
    "${WHATSHAP}" haplotag \
        --reference "${REF}" \
        --sample "${SAMPLE}" \
        --ignore-read-groups \
        --skip-missing-contigs \
        ${HAPLOTAG_REGION:+--regions ${HAPLOTAG_REGION}} \
        --output-haplotag-list "${tags}.$$" \
        -o /dev/null \
        "${vcf}" "${bam}" \
      >>"${logf}" 2>&1
    mv -f "${tags}.$$" "${tags}"
  fi

  # compare has no aggregate row, so sum numerators and denominators ourselves:
  #   switch_rate  = all_switches      / all_assessed_pairs
  #   hamming_rate = blockwise_hamming / covered_variants
  local covered assessed switches hamming
  covered="$(cmp_sum "${cmp}" covered_variants)"
  assessed="$(cmp_sum "${cmp}" all_assessed_pairs)"
  switches="$(cmp_sum "${cmp}" all_switches)"
  hamming="$(cmp_sum "${cmp}" blockwise_hamming)"

  # True read haplotypes come from the phased truth. 20 cells, not 320 runs.
  if [[ ! -s "${truth_tags}" ]]; then
    "${WHATSHAP}" haplotag \
        --reference "${REF}" \
        --sample "${SAMPLE}" \
        --ignore-read-groups \
        --skip-missing-contigs \
        ${HAPLOTAG_REGION:+--regions ${HAPLOTAG_REGION}} \
        --output-haplotag-list "${truth_tags}.$$" \
        -o /dev/null \
        "${TRUTH_PHASED}" "${bam}" \
      >>"${logf}" 2>&1
    mv -f "${truth_tags}.$$" "${truth_tags}"
  fi

  # haplotagging completeness: reads given a haplotype / reads seen
  local tagging
  tagging="$(awk -F'\t' 'NR > 1 { n++; if ($2 != "none") t++ }
      END { printf "%d\t%d\t%.6f", t + 0, n + 0, (n ? t / n : 0) }' "${tags}")"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${caller}" "${platform}" "${depth}" "${phaser}" \
      "$(stats_col "${stats}" variants)" \
      "$(stats_col "${stats}" phased)" \
      "$(stats_col "${stats}" phased_fraction)" \
      "$(stats_col "${stats}" blocks)" \
      "$(stats_col "${stats}" block_n50)" \
      "${switches}" \
      "$(awk -v a="${switches}" -v b="${assessed}" 'BEGIN{printf "%.6f", (b ? a/b : 0)}')" \
      "${hamming}" \
      "$(awk -v a="${hamming}" -v b="${covered}" 'BEGIN{printf "%.6f", (b ? a/b : 0)}')" \
      "$(awk -v a="${hamming}" -v b="${covered}" 'BEGIN{printf "%.6f", (b ? 1 - a/b : 0)}')" \
      "${tagging}" \
      "$(haplotag_accuracy "${truth_tags}" "${tags}")" \
    >"${row}"

  log "  ${name}: $(cut -f13 "${row}") hamming rate"
}

export -f eval_one vcf_name sub_bam bam_key stats_col cmp_sum haplotag_accuracy \
          log warn die
export METRICS_DIR PHASE_DIR LOGS_DIR DATA_DIR REF SAMPLE REFERENCE_NAME \
       HAPLOTAG_REGION \
       REGION_NAME TRUTH_PHASED CHROM_LENGTHS WHATSHAP

log "evaluating ${#CALLERS[@]} x ${#PLATFORMS[@]} x ${#DEPTHS[@]} x ${#PHASERS[@]} runs"

"${PARALLEL}" \
    --joblog "${LOGS_DIR}/parallel_evaluate.log" \
    -j "${EVAL_JOBS}" \
    eval_one {1} {2} {3} {4} \
    ::: "${CALLERS[@]}" \
    ::: "${PLATFORMS[@]}" \
    ::: "${DEPTHS[@]}" \
    ::: "${PHASERS[@]}"

{
  printf 'caller\tplatform\tdepth\tphaser\tvariants\tphased\tphased_fraction\tblocks\tblock_n50'
  printf '\tswitches\tswitch_rate\thamming\thamming_rate\thaplotype_accuracy'
  printf '\treads_tagged\treads_total\ttagging_completeness'
  printf '\treads_correct\treads_compared\thaplotag_accuracy\n'
  cat "${METRICS_DIR}"/*.row.tsv
} >"${METRICS_DIR}/summary.tsv"

log "metrics: ${METRICS_DIR}/summary.tsv"
