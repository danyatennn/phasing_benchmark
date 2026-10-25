#!/usr/bin/env bash
# Experiment A, what eval_run.sh does not measure (plots/variants_eval.py): how
# well the added variants are phased and which reads they let haplotag tag,
# against the same phaser on the SNP-only input. Run after eval_run.sh, which
# writes the haplotag lists, for the experiment-A runs and their truth runs.
#
#   metrics/<run>.variants.tsv          one line per run
#   metrics/summary_variants.tsv        all of them, with a header
#
#   PHASING_DEPTHS=30 bash bin/variants_eval.sh                       truthindel and truthall
#   PHASING_DEPTHS=30 PHASING_CALLERS=truthindel bash bin/variants_eval.sh
#
# A line newer than the run's VCF and both haplotag lists is kept.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

EXP_CALLERS=(${PHASING_CALLERS:-truthindel truthall})
VARIANTS_JOBS="${VARIANTS_JOBS:-6}"
mkdir -p "${METRICS_DIR}" "${LOGS_DIR}/evaluate"

variants_one() {
  set -euo pipefail
  local caller="$1" platform="$2" depth="$3" phaser="$4"
  phaser_applies "${phaser}" "${platform}" || return 0
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
  local vcf="${PHASE_DIR}/${name}.vcf.gz"
  local tags="${METRICS_DIR}/${name}.haplotag.tsv"
  local base_tags="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truth "${phaser}").haplotag.tsv"
  local out="${METRICS_DIR}/${name}.variants.tsv"
  local truth
  case "${caller}" in
    truthindel) truth="${TRUTH_INDEL_PHASED}" ;;
    truthall)   truth="${TRUTH_ALL_PHASED}" ;;
    *) die "variants_eval: '${caller}' is not an experiment-A input" ;;
  esac
  is_done "${vcf}" || return 0
  [[ -s "${tags}" && -s "${base_tags}" ]] || { warn "${name}: haplotag lists missing - run eval_run.sh"; return 0; }
  newer_than "${out}" "${vcf}" && newer_than "${out}" "${tags}" && newer_than "${out}" "${base_tags}" \
    && { log "skip ${name}"; return 0; }

  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/variants_eval.py" \
      --run "${vcf}" --truth-snps "${TRUTH_PHASED}" --truth-extra "${truth}" \
      --tags "${tags}" --base-tags "${base_tags}" --alleles "${READS_DIR}/${cell}_truth.reads.tsv.gz" \
      --label "${REFERENCE_NAME} ${caller} ${platform} ${depth} ${phaser}" --out "${out}" \
    2>>"${LOGS_DIR}/evaluate/${name}.variants.log" \
    || die "variants_eval failed for ${name}, see ${LOGS_DIR}/evaluate/${name}.variants.log"
  tail -1 "${LOGS_DIR}/evaluate/${name}.variants.log" >&2
}

export -f variants_one vcf_name bam_key phaser_applies is_done newer_than log warn die
export METRICS_DIR PHASE_DIR READS_DIR LOGS_DIR REFERENCE_NAME REGION_NAME SAMPLE \
       TRUTH_PHASED TRUTH_INDEL_PHASED TRUTH_ALL_PHASED PLOTS_PYTHON PROJECT_ROOT

log "${SAMPLE} ${REFERENCE_NAME}: added-variant phase and rescued reads for ${EXP_CALLERS[*]} x ${PLATFORMS[*]} x ${DEPTHS[*]} x ${PHASERS[*]}"
"${PARALLEL}" -j "${VARIANTS_JOBS}" variants_one {1} {2} {3} {4} \
    ::: "${EXP_CALLERS[@]}" ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" ::: "${PHASERS[@]}" || true

summary="${METRICS_DIR}/summary_variants.tsv"
{
  printf 'reference\tcaller\tplatform\tdepth\tphaser'
  printf '\tindel_in_input\tindel_phased\tindel_oriented\tindel_correct\tindel_wrong'
  printf '\tother_in_input\tother_phased\tother_oriented\tother_correct\tother_wrong'
  printf '\treads_tagged_base\treads_tagged\treads_rescued\treads_rescued_no_snp\treads_lost\n'
  files=("${METRICS_DIR}"/*.variants.tsv)
  [[ -e "${files[0]}" ]] && cat "${files[@]}"
  true
} >"${summary}.$$"
mv -f "${summary}.$$" "${summary}"
log "variants: ${summary}"
