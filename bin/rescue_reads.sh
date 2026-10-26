#!/usr/bin/env bash
# Reads tagged on T2T-CHM13 but not on GRCh38, and the other way round, per
# platform (plots/rescue_reads.py), for figR9 of plots/figures_summary.py:
#
#   metrics_reads/rescue/<platform>_<depth>x_<input>.{reads.tsv.gz,upset.tsv,summary.tsv,
#                                                     regions_T2T.tsv,regions_GRCh38.tsv}
#
# Needs, on both references: the read tables (bin/read_tables.sh), the allele
# matrices of the truth input (bin/allele_matrix.sh), the haplotag lists and
# phased VCFs of the phasers, and the truth haplotag lists (eval_run.sh).
# Holds every read name of a cell in memory: run it on a host with room.
#
#   bash bin/rescue_reads.sh
#   RESCUE_PHASERS="whatshap longphase hapcut2" PHASING_PLATFORMS=ont5k bash bin/rescue_reads.sh

set -euo pipefail

export REFERENCE_NAME=GRCh38
source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

RESCUE_PHASERS=(${RESCUE_PHASERS:-whatshap longphase})
RESCUE_DEPTH="${RESCUE_DEPTH:-30}"
RESCUE_CALLER="${RESCUE_CALLER:-truth}"
OUT_DIR="${METRICS_READS_DIR}/rescue"
mkdir -p "${OUT_DIR}" "${LOGS_DIR}/reads"

truth_vcf() { [[ "$1" == GRCh38 ]] && echo "${DERIVED_DIR}/truth.het_snps.phased.vcf.gz" \
                                    || echo "${DERIVED_DIR}/truth_$1.het_snps.phased.vcf.gz"; }
cell_of() { printf '%s_%s_%s_%sx_%s\n' "$1" "$2" "${SAMPLE}" "${RESCUE_DEPTH}" "${REGION_NAME}"; }

strata=()
for s in $(STRATA_TIER=1 strata_names); do strata+=(--stratum "${s}"); done

for platform in "${PLATFORMS[@]}"; do
  out="${OUT_DIR}/${platform}_${RESCUE_DEPTH}x_${RESCUE_CALLER}"
  args=()
  missing=0
  for ref in GRCh38 T2T; do
    k="$(tr '[:upper:]' '[:lower:]' <<<"${ref}")"
    cell="$(cell_of "${platform}" "${ref}")"
    args+=(--reads-"${k}" "${READS_DIR}/${cell}.reads.tsv.gz"
           --alleles-"${k}" "${READS_DIR}/${cell}_truth.reads.tsv.gz"
           --truth-vcf-"${k}" "$(truth_vcf "${ref}")"
           --truth-tags-"${k}" "${METRICS_DIR}/${cell}_truthtags.tsv"
           --strata-"${k}" "${PREPARED_DIR}")
    for f in "${READS_DIR}/${cell}.reads.tsv.gz" "$(truth_vcf "${ref}")" "${METRICS_DIR}/${cell}_truthtags.tsv"; do
      [[ -s "${f}" ]] || { warn "${platform}: missing ${f}"; missing=1; }
    done
    for phaser in "${RESCUE_PHASERS[@]}"; do
      name="${cell}_${RESCUE_CALLER}_${phaser}"
      tags="${METRICS_DIR}/${name}.haplotag.tsv"
      vcf="${PHASE_DIR}/${name}.vcf.gz"
      [[ -s "${tags}" && -s "${vcf}" ]] || { warn "${platform}: no ${ref} ${phaser} run"; missing=1; }
      args+=(--run "${ref}" "${phaser}" "${tags}" "${vcf}")
      newer_than "${out}.summary.tsv" "${tags}" && newer_than "${out}.summary.tsv" "${vcf}" || rm -f "${out}.summary.tsv"
    done
  done
  [[ "${missing}" == 0 ]] || { warn "skip ${platform}"; continue; }
  [[ -s "${out}.summary.tsv" ]] && { log "skip ${platform}: up to date"; continue; }
  log "${platform} ${RESCUE_DEPTH}x: reads tagged on one reference only (${RESCUE_PHASERS[*]})"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/rescue_reads.py" "${args[@]}" "${strata[@]}" \
      --min-snps "${READACC_MIN_SNPS}" --out "${out}" \
    2>>"${LOGS_DIR}/reads/rescue_${platform}_${RESCUE_DEPTH}x.log" \
    || die "rescue analysis failed for ${platform}, see ${LOGS_DIR}/reads/rescue_${platform}_${RESCUE_DEPTH}x.log"
  tail -n "$(( 2 * ${#RESCUE_PHASERS[@]} + 1 ))" "${LOGS_DIR}/reads/rescue_${platform}_${RESCUE_DEPTH}x.log" >&2
done
log "rescue: ${OUT_DIR}"
