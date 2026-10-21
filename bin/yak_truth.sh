#!/usr/bin/env bash
# Parental-haplotype truth for the reads from the parents' k-mers: yak triobin
# with pat.yak and mat.yak of trio_kmers.sh. A read is paternal or maternal if
# it holds k-mers of one parent only; a second, alignment-free check of the
# Q100 truth (bin/matpat_truth.sh), and the only one for HG005.
#
#   reads/<cell>.yak.tsv.gz      read label(PAT/MAT/none) and yak's counts
#
# Kept when newer than the BAM and the k-mer databases.
#
#   PHASING_DEPTHS=30 bash bin/yak_truth.sh
#   SAMPLE=HG005 bash bin/yak_truth.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

YAK_THREADS="${YAK_THREADS:-16}"
export PATH="${TRIO_BIN}:${PATH}"
require_cmds yak
require_file "${TRIO_DIR}/pat.yak" "paternal k-mers - run trio_kmers.sh"
require_file "${TRIO_DIR}/mat.yak" "maternal k-mers - run trio_kmers.sh"
mkdir -p "${READS_DIR}" "${LOGS_DIR}/reads"

for platform in "${PLATFORMS[@]}"; do
  for depth in "${DEPTHS[@]}"; do
    cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
    bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
    out="${READS_DIR}/${cell}.yak.tsv.gz"
    [[ -s "${bam}" ]] || { warn "no BAM for ${cell}"; continue; }
    newer_than "${out}" "${bam}" && newer_than "${out}" "${TRIO_DIR}/pat.yak" \
      && newer_than "${out}" "${TRIO_DIR}/mat.yak" && { log "skip ${cell}"; continue; }

    log "${cell}: yak triobin"
    logf="${LOGS_DIR}/reads/${cell}.yak.log"
    # yak triobin: read name, then p (paternal), m (maternal), a (both), 0 (no
    # parental k-mer), then its counts; p and m become PAT and MAT
    "${SAMTOOLS}" fastq -@ 4 -n -F 0x900 "${bam}" 2>"${logf}" \
      | yak triobin -t "${YAK_THREADS}" "${TRIO_DIR}/pat.yak" "${TRIO_DIR}/mat.yak" /dev/stdin 2>>"${logf}" \
      | awk -F'\t' -v OFS='\t' 'BEGIN { print "#read", "label", "yak" }
          { l = ($2 == "p") ? "PAT" : ($2 == "m") ? "MAT" : "none"; $2 = l OFS $2; print }' \
      | "${BGZIP}" -c >"${out}.tmp" \
      || die "yak triobin failed for ${cell}, see ${logf}"
    mv -f "${out}.tmp" "${out}"
    log "${cell}: $(gzip -dc "${out}" | awk -F'\t' 'NR > 1 { n[$2]++ } END { printf "PAT %d, MAT %d, none %d", n["PAT"], n["MAT"], n["none"] }')"
  done
done
log "yak truth: ${READS_DIR}/*.yak.tsv.gz"
