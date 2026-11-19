#!/usr/bin/env bash
# The state of each read on GRCh38 against T2T-CHM13 (plots/read_transitions.py),
# per platform, for figR20 and figR21 of plots/figures_summary.py:
#
#   metrics_reads/transitions/<platform>_<depth>x_<input>.tsv   phaser grch38 t2t region reads
#
# region: both, t2t, grch38 or other - in all difficult regions of both
# references, of T2T-CHM13 only, of GRCh38 only or of neither, each by the read's
# alignment on that reference (bin/prepare_strata.sh). A table of the earlier
# layout (difficult/other by T2T-CHM13 alone) is set aside as .t2t_regions and
# made again.
#
# Needs, on both references: the read tables (bin/read_tables.sh), the truth
# allele matrices (bin/allele_matrix.sh), the consensus parental labels
# (bin/consensus_truth.sh, either reference), the truth haplotag lists and the
# haplotag lists and phased VCFs of the phasers (eval_run.sh). HG002 only (the
# Q100 labels). Holds every read name of a cell in memory: run it on a host
# with room.
#
#   bash bin/read_transitions.sh
#   TRANSITION_PHASERS="whatshap longphase" PHASING_PLATFORMS=ont5k bash bin/read_transitions.sh

set -euo pipefail

export REFERENCE_NAME=GRCh38
source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

[[ "${SAMPLE}" == HG002 ]] || die "the parental labels are HG002's; there are none for ${SAMPLE}"
TRANSITION_PHASERS=(${TRANSITION_PHASERS:-whatshap})
TRANSITION_DEPTH="${TRANSITION_DEPTH:-30}"
TRANSITION_CALLER="${TRANSITION_CALLER:-truth}"
OUT_DIR="${METRICS_READS_DIR}/transitions"
mkdir -p "${OUT_DIR}" "${LOGS_DIR}/reads"

truth_vcf() { [[ "$1" == GRCh38 ]] && echo "${DERIVED_DIR}/truth.het_snps.phased.vcf.gz" \
                                    || echo "${DERIVED_DIR}/truth_$1.het_snps.phased.vcf.gz"; }
bench_bed() { [[ "$1" == GRCh38 ]] && echo "${HG002_GRCH38_V5Q_BED}" || echo "${HG002_CHM13_V5Q_BED}"; }

for platform in "${PLATFORMS[@]}"; do
  out="${OUT_DIR}/${platform}_${TRANSITION_DEPTH}x_${TRANSITION_CALLER}.tsv"
  args=()
  missing=0
  for ref in GRCh38 T2T; do
    k="$(tr '[:upper:]' '[:lower:]' <<<"${ref}")"
    cell="$(bam_key "${platform}" "${ref}" "${TRANSITION_DEPTH}")"
    labels="${READS_DIR}/${cell}.consensus.tsv.gz"
    [[ -s "${labels}" ]] && args+=(--labels "${labels}")
    args+=(--reads-"${k}" "${READS_DIR}/${cell}.reads.tsv.gz"
           --alleles-"${k}" "${READS_DIR}/${cell}_truth.reads.tsv.gz"
           --truth-vcf-"${k}" "$(truth_vcf "${ref}")"
           --truth-tags-"${k}" "${METRICS_DIR}/${cell}_truthtags.tsv"
           --bench-"${k}" "$(bench_bed "${ref}")")
    for f in "${READS_DIR}/${cell}.reads.tsv.gz" "${READS_DIR}/${cell}_truth.reads.tsv.gz" \
             "$(truth_vcf "${ref}")" "${METRICS_DIR}/${cell}_truthtags.tsv" "$(bench_bed "${ref}")"; do
      [[ -s "${f}" ]] || { warn "${platform}: missing ${f}"; missing=1; }
    done
    for phaser in "${TRANSITION_PHASERS[@]}"; do
      name="${cell}_${TRANSITION_CALLER}_${phaser}"
      tags="${METRICS_DIR}/${name}.haplotag.tsv"
      vcf="${PHASE_DIR}/${name}.vcf.gz"
      [[ -s "${tags}" && -s "${vcf}" ]] || { warn "${platform}: no ${ref} ${phaser} run"; missing=1; }
      args+=(--run "${ref}" "${phaser}" "${tags}" "${vcf}")
      newer_than "${out}" "${tags}" && newer_than "${out}" "${vcf}" || rm -f "${out}"
    done
  done
  [[ " ${args[*]} " == *" --labels "* ]] || { warn "${platform}: no reads/<cell>.consensus.tsv.gz - run bin/consensus_truth.sh"; missing=1; }
  difficult="$(REFERENCE_NAME=T2T stratum_bed alldifficultregions)"
  difficult_grch38="$(REFERENCE_NAME=GRCh38 stratum_bed alldifficultregions)"
  [[ -s "${difficult}" ]] || { warn "${platform}: missing ${difficult}"; missing=1; }
  [[ -s "${difficult_grch38}" ]] || { warn "${platform}: missing ${difficult_grch38}"; missing=1; }
  [[ "${missing}" == 0 ]] || { warn "skip ${platform}"; continue; }
  # a table of the earlier layout, without the region column, is made again
  [[ -s "${out}" ]] && ! head -1 "${out}" | grep -q region && rm -f "${out}"
  if [[ -s "${out}" ]] && awk -F'\t' 'NR > 1 && $4 == "difficult" { f = 1; exit } END { exit !f }' "${out}"; then
    mv -f "${out}" "${out}.t2t_regions"
  fi
  [[ -s "${out}" ]] && { log "skip ${platform}: up to date"; continue; }
  log "${platform} ${TRANSITION_DEPTH}x: read states on GRCh38 and T2T-CHM13 (${TRANSITION_PHASERS[*]})"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/read_transitions.py" "${args[@]}" \
      --difficult "${difficult}" --difficult-grch38 "${difficult_grch38}" --min-snps "${READACC_MIN_SNPS}" --out "${out}" \
    2>>"${LOGS_DIR}/reads/transitions_${platform}_${TRANSITION_DEPTH}x.log" \
    || die "read transitions failed for ${platform}, see ${LOGS_DIR}/reads/transitions_${platform}_${TRANSITION_DEPTH}x.log"
  tail -n "$(( ${#TRANSITION_PHASERS[@]} + 2 ))" "${LOGS_DIR}/reads/transitions_${platform}_${TRANSITION_DEPTH}x.log" >&2
done
log "read transitions: ${OUT_DIR}"
