#!/usr/bin/env bash
# The allele matrix of every cell and phaser input (plots/allele_matrix.py):
# which allele each read carries at each het SNP of the input, and the votes
# of each read for the two haplotypes of every phaser's blocks.
#
#   reads/<cell>_<caller>.sites.tsv.gz, .links.tsv.gz, .reads.tsv.gz,
#   reads/<cell>_<caller>.votes.<phaser>.tsv.gz
#
# <caller> is truth (the phased truth's sites) or a caller (its het SNPs). Kept
# when newer than the BAM, the input VCF and every phased VCF it covers.
#
#   PHASING_DEPTHS=30 bash bin/allele_matrix.sh
#   PHASING_DEPTHS=30 PHASING_CALLERS=truth bash bin/allele_matrix.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

ALLELE_JOBS="${ALLELE_JOBS:-4}"
ALLELE_THREADS="${ALLELE_THREADS:-8}"
mkdir -p "${READS_DIR}" "${LOGS_DIR}/reads"

matrix_one() {
  set -euo pipefail
  local caller="$1" platform="$2" depth="$3"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local sites
  if [[ "${caller}" == truth ]]; then
    sites="${TRUTH_PHASED}"
  else
    sites="${VCF_DIR}/$(vcf_name "${platform}" "${depth}" "${caller}").het_snps.vcf.gz"
  fi
  local out="${READS_DIR}/${cell}_${caller}"
  [[ -s "${bam}" && -s "${sites}" ]] || { warn "${cell} ${caller}: no BAM or no input VCF"; return 0; }

  local phased=() fresh=1 phaser vcf
  newer_than "${out}.reads.tsv.gz" "${bam}" && newer_than "${out}.reads.tsv.gz" "${sites}" || fresh=0
  for phaser in ${PHASER_LIST}; do
    vcf="${PHASE_DIR}/$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}").vcf.gz"
    is_done "${vcf}" || continue
    phased+=(--phased "${phaser}=${vcf}")
    newer_than "${out}.votes.${phaser}.tsv.gz" "${vcf}" || fresh=0
  done
  [[ "${fresh}" == 1 ]] && { log "skip ${cell} ${caller}"; return 0; }

  local lock="${out}.lock"
  mkdir "${lock}" 2>/dev/null || { log "busy elsewhere ${cell} ${caller}"; return 0; }
  trap "rmdir '${lock}' 2>/dev/null" EXIT INT TERM HUP

  log "${cell} ${caller}: allele matrix, $(( ${#phased[@]} / 2 )) phasers"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/allele_matrix.py" \
      --bam "${bam}" --sites "${sites}" ${phased[@]+"${phased[@]}"} \
      --samtools "${SAMTOOLS}" --out "${out}" --threads "${ALLELE_THREADS}" \
    2>"${LOGS_DIR}/reads/${cell}_${caller}.alleles.log" \
    || die "allele matrix failed for ${cell} ${caller}, see ${LOGS_DIR}/reads/${cell}_${caller}.alleles.log"
}

PHASER_LIST="${PHASERS[*]}"
export -f matrix_one bam_key sub_bam vcf_name is_done newer_than log warn die
export READS_DIR LOGS_DIR DATA_DIR VCF_DIR PHASE_DIR REFERENCE_NAME REGION_NAME SAMPLE \
       TRUTH_PHASED PLOTS_PYTHON PROJECT_ROOT SAMTOOLS ALLELE_THREADS PHASER_LIST

log "${SAMPLE} ${REFERENCE_NAME}: allele matrices for ${CALLERS[*]} x ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${ALLELE_JOBS}" matrix_one {1} {2} {3} \
    ::: "${CALLERS[@]}" ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true
log "allele matrices: ${READS_DIR}"
