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

# truth_total is the denominator of fnr and recall_phased in every row; an empty
# value would silently turn fnr into 0
TRUTH_TOTAL="$("${BCFTOOLS}" index -n "${TRUTH_PHASED}")" \
  || die "cannot count the truth sites in ${TRUTH_PHASED}"
[[ "${TRUTH_TOTAL}" =~ ^[0-9]+$ && "${TRUTH_TOTAL}" -gt 0 ]] \
  || die "truth has no sites: '${TRUTH_TOTAL}' in ${TRUTH_PHASED}"

if [[ ! -s "${CHROM_LENGTHS}" ]]; then
  cut -f1,2 "${REF}.fai" >"${CHROM_LENGTHS}.$$"
  mv -f "${CHROM_LENGTHS}.$$" "${CHROM_LENGTHS}"
fi

eval_one() {
  set -euo pipefail

  local caller="$1" platform="$2" depth="$3" phaser="$4"
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
  local vcf="${PHASE_DIR}/${name}.vcf.gz"
  local row="${METRICS_DIR}/${name}.row.tsv"

  [[ -s "${row}" ]] && row_is_valid "${row}" && { log "skip ${name}"; return 0; }

  local lock="${row}.lock"
  mkdir "${lock}" 2>/dev/null || { log "busy elsewhere ${name} ($(cat "${lock}/owner" 2>/dev/null))"; return 0; }
  printf '%s %s %s\n' "$(hostname)" "$$" "$(date '+%F %H:%M')" >"${lock}/owner"
  trap "rm -f '${lock}/owner'; rmdir '${lock}' 2>/dev/null" EXIT INT TERM HUP

  # A row is final only if it is whole. One written before these checks
  # existed may not be: it is kept as .bad and the run is computed again.
  if [[ -s "${row}" ]]; then
    row_is_valid "${row}" && { log "skip ${name}"; return 0; }
    warn "${name}: malformed row ($(awk -F'\t' 'NR == 1 { print NF }' "${row}") fields, expected ${ROW_FIELDS}) - recomputing"
    mv -f "${row}" "${row}.bad.$(date +%Y%m%d%H%M%S)"
  fi
  [[ -s "${vcf}" ]] || { warn "no phased VCF for ${name}"; return 0; }

  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local logf="${LOGS_DIR}/evaluate/${name}.log"
  local stats="${METRICS_DIR}/${name}.stats.tsv"
  local cmp="${METRICS_DIR}/${name}.compare.tsv"
  local tags="${METRICS_DIR}/${name}.haplotag.tsv"
  local truth_tags="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truthtags).tsv"

  log "evaluating ${name}"
  : >"${logf}"

  # stats and haplotag are reused only if they were made from this VCF: after
  # re-phasing a run they are older than it and are made again
  if ! newer_than "${stats}" "${vcf}"; then
    "${WHATSHAP}" stats \
        --tsv "${stats}.$$" \
        --chr-lengths "${CHROM_LENGTHS}" \
        "${vcf}" \
      >"${METRICS_DIR}/${name}.stats.txt" 2>>"${logf}"
    mv -f "${stats}.$$" "${stats}"
  fi

  "${WHATSHAP}" compare \
      --sample "${SAMPLE}" \
      --names truth,"${phaser}" \
      --tsv-pairwise "${cmp}" \
      --only-snvs \
      "${TRUTH_PHASED}" "${vcf}" \
    >>"${logf}" 2>&1

  if ! newer_than "${tags}" "${vcf}"; then
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

  #   switch_rate = all_switches / all_assessed_pairs
  #   hamming_rate = blockwise_hamming / covered_variants
  local covered assessed switches hamming
  covered="$(cmp_sum "${cmp}" covered_variants)"
  assessed="$(cmp_sum "${cmp}" all_assessed_pairs)"
  switches="$(cmp_sum "${cmp}" all_switches)"
  hamming="$(cmp_sum "${cmp}" blockwise_hamming)"

  # shared by every run of the cell; made again only if the truth or the BAM
  # changed after it
  if ! newer_than "${truth_tags}" "${TRUTH_PHASED}" || ! newer_than "${truth_tags}" "${bam}"; then
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

  local emitted shared correct
  emitted="$("${BCFTOOLS}" index -n "${vcf}")"
  shared="$("${BCFTOOLS}" isec -n=2 -w1 -c none "${TRUTH_PHASED}" "${vcf}" \
    | awk '!/^#/ { n++ } END { print n + 0 }')"
  correct=$((covered - hamming))

  # Every value is computed into a variable first. A command substitution in the
  # arguments of printf below would fail without stopping the script (set -e
  # does not reach it), and the row would be written with fields missing.
  local tagging accuracy
  tagging="$(awk -F'\t' 'NR > 1 { n++; if ($2 != "none") t++ }
      END { printf "%d\t%d\t%.6f", t + 0, n + 0, (n ? t / n : 0) }' "${tags}")"
  accuracy="$(haplotag_accuracy "${truth_tags}" "${tags}")"

  local variants phased phased_fraction blocks block_n50
  variants="$(stats_col "${stats}" variants)"
  phased="$(stats_col "${stats}" phased)"
  phased_fraction="$(stats_col "${stats}" phased_fraction)"
  blocks="$(stats_col "${stats}" blocks)"
  block_n50="$(stats_col "${stats}" block_n50)"

  local switch_rate hamming_rate haplotype_accuracy fnr fdr recall_phased precision_phased
  switch_rate="$(awk -v a="${switches}" -v b="${assessed}" 'BEGIN{printf "%.6f", (b ? a/b : 0)}')"
  hamming_rate="$(awk -v a="${hamming}" -v b="${covered}" 'BEGIN{printf "%.6f", (b ? a/b : 0)}')"
  haplotype_accuracy="$(awk -v a="${hamming}" -v b="${covered}" 'BEGIN{printf "%.6f", (b ? 1 - a/b : 0)}')"
  fnr="$(awk -v s="${shared}" -v t="${TRUTH_TOTAL}" 'BEGIN{printf "%.6f", (t ? 1 - s/t : 0)}')"
  fdr="$(awk -v s="${shared}" -v m="${emitted}"    'BEGIN{printf "%.6f", (m ? 1 - s/m : 0)}')"
  recall_phased="$(awk -v c="${correct}" -v t="${TRUTH_TOTAL}" 'BEGIN{printf "%.6f", (t ? c/t : 0)}')"
  precision_phased="$(awk -v c="${correct}" -v m="${emitted}"     'BEGIN{printf "%.6f", (m ? c/m : 0)}')"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${REFERENCE_NAME}" "${REGION_NAME}" \
      "${caller}" "${platform}" "${depth}" "${phaser}" \
      "${variants}" "${phased}" "${phased_fraction}" "${blocks}" "${block_n50}" \
      "${switches}" "${switch_rate}" "${hamming}" "${hamming_rate}" "${haplotype_accuracy}" \
      "${covered}" "${assessed}" "${TRUTH_TOTAL}" "${emitted}" "${shared}" "${correct}" \
      "${fnr}" "${fdr}" "${recall_phased}" "${precision_phased}" \
      "${tagging}" \
      "${accuracy}" \
    >"${row}.$$"
  row_is_valid "${row}.$$" \
    || { mv -f "${row}.$$" "${row}.bad.$(date +%Y%m%d%H%M%S)"; die "${name}: malformed row, kept as ${row}.bad.*"; }
  mv -f "${row}.$$" "${row}"

  log "  ${name}: $(cut -f15 "${row}") hamming rate"
}

export -f eval_one vcf_name sub_bam bam_key stats_col cmp_sum haplotag_accuracy \
          newer_than row_is_valid log warn die
export TRUTH_TOTAL BCFTOOLS ROW_FIELDS
export METRICS_DIR PHASE_DIR LOGS_DIR DATA_DIR REF SAMPLE REFERENCE_NAME \
       HAPLOTAG_REGION \
       REGION_NAME TRUTH_PHASED CHROM_LENGTHS WHATSHAP

log "evaluating ${#CALLERS[@]} x ${#PLATFORMS[@]} x ${#DEPTHS[@]} x ${#PHASERS[@]} runs"

joblog="${LOGS_DIR}/parallel_evaluate_${RUN_TAG}.log"
keep_joblog "${joblog}"

"${PARALLEL}" \
    --joblog "${joblog}" \
    -j "${EVAL_JOBS}" \
    eval_one {1} {2} {3} {4} \
    ::: "${CALLERS[@]}" \
    ::: "${PLATFORMS[@]}" \
    ::: "${DEPTHS[@]}" \
    ::: "${PHASERS[@]}" \
  || true

failed="$(awk 'NR > 1 && $7 != 0 { n++ } END { print n + 0 }' "${joblog}")"
[[ "${failed}" -eq 0 ]] || warn "${failed} runs failed, logs in ${LOGS_DIR}/evaluate"

{
  printf 'reference\tregion\tcaller\tplatform\tdepth\tphaser\tvariants\tphased\tphased_fraction\tblocks\tblock_n50'
  printf '\tswitches\tswitch_rate\thamming\thamming_rate\thaplotype_accuracy'
  printf '\tcovered\tassessed\ttruth_total\temitted\tshared\tcorrect'
  printf '\tfnr\tfdr\trecall_phased\tprecision_phased'
  printf '\treads_tagged\treads_total\ttagging_completeness'
  printf '\treads_correct\treads_compared\thaplotag_accuracy\n'
  cat "${METRICS_DIR}"/*.row.tsv
} >"${METRICS_DIR}/summary.tsv.$$"
mv -f "${METRICS_DIR}/summary.tsv.$$" "${METRICS_DIR}/summary.tsv"

log "metrics: ${METRICS_DIR}/summary.tsv"
