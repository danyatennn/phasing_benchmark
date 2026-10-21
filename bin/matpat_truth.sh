#!/usr/bin/env bash
# Parental-haplotype truth for the reads of HG002 cells: every read aligned to
# the diploid HG002 Q100 assembly (both parental copies of each chromosome,
# Q100_FASTA in config.sh) and labelled with the copy it fits better
# (plots/matpat_assign.py). It does not depend on any VCF or on whatshap.
#
#   reads/<cell>.matpat.tsv.gz       read label(MAT/PAT/none) chrom as_mat as_pat
#
# The minimap2 index is built once per preset next to the FASTA: lr:hq for ONT,
# map-hifi for PacBio, as for the T2T alignments. A table newer than its BAM is
# kept. Then: ANCHOR_TRUTH=matpat bash readacc_anchor.sh
#
#   PHASING_DEPTHS=30 bash bin/matpat_truth.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

[[ "${SAMPLE}" == HG002 ]] || die "the Q100 assembly is HG002's; there is none for ${SAMPLE}"
require_file "${Q100_FASTA}" "diploid HG002 Q100 assembly (download it, or set Q100_FASTA)"
MATPAT_THREADS="${MATPAT_THREADS:-32}"
MATPAT_MIN_DIFF="${MATPAT_MIN_DIFF:-10}"
mkdir -p "${READS_DIR}" "${LOGS_DIR}/reads"

preset_of() {
  case "$1" in
    ont4k|ont5k) printf 'lr:hq\n' ;;
    hifi|revio)  printf 'map-hifi\n' ;;
  esac
}

for platform in "${PLATFORMS[@]}"; do
  for depth in "${DEPTHS[@]}"; do
    cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
    bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
    out="${READS_DIR}/${cell}.matpat.tsv.gz"
    [[ -s "${bam}" ]] || { warn "no BAM for ${cell}"; continue; }
    newer_than "${out}" "${bam}" && newer_than "${out}" "${Q100_FASTA}" && { log "skip ${cell}"; continue; }

    preset="$(preset_of "${platform}")"
    mmi="${Q100_FASTA%.gz}"
    mmi="${mmi%.fa*}.${preset//:/_}.mmi"
    if [[ ! -s "${mmi}" ]]; then
      log "minimap2 index of $(basename "${Q100_FASTA}") for ${preset}"
      "${MINIMAP2}" -t "${MATPAT_THREADS}" -x "${preset}" -d "${mmi}.tmp" "${Q100_FASTA}"
      mv -f "${mmi}.tmp" "${mmi}"
    fi

    log "${cell}: reads -> Q100 (${preset})"
    logf="${LOGS_DIR}/reads/${cell}.matpat.log"
    "${SAMTOOLS}" fastq -@ 4 -n -F 0x900 "${bam}" 2>"${logf}" \
      | "${MINIMAP2}" -t "${MATPAT_THREADS}" -c --secondary=yes -N 10 -p 0.5 -x "${preset}" "${mmi}" - 2>>"${logf}" \
      | "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/matpat_assign.py" --out "${out}" --min-diff "${MATPAT_MIN_DIFF}" 2>>"${logf}" \
      || die "matpat truth failed for ${cell}, see ${logf}"
    tail -2 "${logf}" >&2
  done
done
log "matpat truth: ${READS_DIR}/*.matpat.tsv.gz"
