#!/usr/bin/env bash
# Read-level accuracy with every phase block oriented by its anchor read: the
# longest read of the block that has a truth label (plots/readacc_anchor.py
# explains the rule). readacc_run.sh orients blocks by their SNPs instead;
# both are kept, in different files, to be compared.
#
#   metrics/<run>.anchor_<truth>_<rules>.tsv          one line per run
#   metrics/summary_anchor_<truth>_<rules>.tsv        all of them, with a header
#
#   <truth>  ANCHOR_TRUTH: truthtags (default) - whatshap haplotag against the
#            phased truth; matpat - bin/matpat_truth.sh; yak - bin/yak_truth.sh
#   <rules>  READACC_MIN_SNPS and READACC_AUTOSOMES of config.sh, as in
#            readacc_run.sh: auto_min2 by default
#
# Needs the read table of the cell (bin/read_tables.sh) and, like
# readacc_run.sh, a finished row of the run. A line newer than the run's
# haplotag list, its VCF and the truth labels is kept.
#
#   bash readacc_anchor.sh
#   ANCHOR_TRUTH=matpat PHASING_DEPTHS=30 bash readacc_anchor.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

ANCHOR_TRUTH="${ANCHOR_TRUTH:-truthtags}"
ANCHOR_JOBS="${ANCHOR_JOBS:-6}"
case "${ANCHOR_TRUTH}" in
  truthtags|matpat|yak) ;;
  *) die "ANCHOR_TRUTH must be truthtags, matpat or yak, not '${ANCHOR_TRUTH}'" ;;
esac

rules="$(readacc_suffix)"
rules="${rules#readacc}"
rules="${rules#_}"
[[ -n "${rules}" ]] || rules="all"
ANCHOR_SUFFIX="anchor_${ANCHOR_TRUTH}_${rules}"
mkdir -p "${TMP_DIR}" "${LOGS_DIR}/readacc"

anchor_cell() {
  set -euo pipefail
  local platform="$1" depth="$2"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local table="${READS_DIR}/${cell}.reads.tsv.gz"
  local truth
  case "${ANCHOR_TRUTH}" in
    truthtags) truth="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truthtags).tsv" ;;
    matpat)    truth="${READS_DIR}/${cell}.matpat.tsv.gz" ;;
    yak)       truth="${READS_DIR}/${cell}.yak.tsv.gz" ;;
  esac
  [[ -s "${table}" ]] || { warn "${cell}: no read table - run bin/read_tables.sh"; return 0; }
  [[ -s "${truth}" ]] || { warn "${cell}: no ${ANCHOR_TRUTH} labels ${truth}"; return 0; }

  local jobs; jobs="$(mktemp -p "${TMP_DIR}")"
  local caller phaser name row tags vcf out
  for caller in ${CALLER_LIST}; do
    for phaser in ${PHASER_LIST}; do
      name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
      row="${METRICS_DIR}/${name}.row.tsv"
      tags="${METRICS_DIR}/${name}.haplotag.tsv"
      vcf="${PHASE_DIR}/${name}.vcf.gz"
      out="${METRICS_DIR}/${name}.${ANCHOR_SUFFIX}.tsv"
      [[ -s "${row}" && ! -d "${row}.lock" && -s "${tags}" && -s "${vcf}" ]] || continue
      newer_than "${out}" "${tags}" && newer_than "${out}" "${vcf}" && newer_than "${out}" "${truth}" \
        && newer_than "${out}" "${table}" && continue
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${REFERENCE_NAME}" "${REGION_NAME}" \
          "${caller}" "${platform}" "${depth}" "${phaser}" "${tags}" "${vcf}" "${out}" >>"${jobs}"
    done
  done
  if [[ ! -s "${jobs}" ]]; then
    rm -f "${jobs}"
    log "skip ${cell}"
    return 0
  fi
  log "${cell}: $(wc -l <"${jobs}") runs, anchor orientation, ${ANCHOR_TRUTH} truth"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/readacc_anchor.py" \
      --reads-table "${table}" --truth "${truth}" --truth-source "${ANCHOR_TRUTH}" \
      --min-snps "${READACC_MIN_SNPS}" --autosomes "${READACC_AUTOSOMES}" --jobs "${jobs}" \
    2>>"${LOGS_DIR}/readacc/${cell}.${ANCHOR_SUFFIX}.log" \
    || { rm -f "${jobs}"; die "anchor accuracy failed for ${cell}, see ${LOGS_DIR}/readacc/${cell}.${ANCHOR_SUFFIX}.log"; }
  rm -f "${jobs}"
}

CALLER_LIST="${CALLERS[*]}"
PHASER_LIST="${PHASERS[*]}"
export -f anchor_cell bam_key vcf_name newer_than log warn die
export METRICS_DIR PHASE_DIR READS_DIR LOGS_DIR TMP_DIR REFERENCE_NAME REGION_NAME SAMPLE \
       ANCHOR_TRUTH ANCHOR_SUFFIX READACC_MIN_SNPS READACC_AUTOSOMES PLOTS_PYTHON PROJECT_ROOT \
       CALLER_LIST PHASER_LIST

log "${SAMPLE} ${REFERENCE_NAME}: anchor-read accuracy (${ANCHOR_SUFFIX}) for ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${ANCHOR_JOBS}" anchor_cell {1} {2} ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true

{
  printf 'reference\tregion\tcaller\tplatform\tdepth\tphaser\ttruth_source'
  printf '\treads_truth_labelled\treads_tagged\treads_scored\treads_correct'
  printf '\tanchors\tblocks_flipped\tread_precision\tread_recall\n'
  cat "${METRICS_DIR}"/*."${ANCHOR_SUFFIX}".tsv 2>/dev/null || true
} >"${METRICS_DIR}/summary_${ANCHOR_SUFFIX}.tsv.$$"
mv -f "${METRICS_DIR}/summary_${ANCHOR_SUFFIX}.tsv.$$" "${METRICS_DIR}/summary_${ANCHOR_SUFFIX}.tsv"
log "anchor accuracy: ${METRICS_DIR}/summary_${ANCHOR_SUFFIX}.tsv"
