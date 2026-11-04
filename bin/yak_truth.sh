#!/usr/bin/env bash
# Parental-haplotype truth for the reads from the parents' k-mers: yak triobin
# with pat.yak and mat.yak of trio_kmers.sh. A read is paternal or maternal if
# it holds k-mers of one parent only; a second, alignment-free check of the
# Q100 truth (bin/matpat_truth.sh), and the only one for HG005.
#
#   reads/<cell>.yak.tsv.gz      read label(PAT/MAT/none) and yak's counts
#   reads/<cell>.kmer.tsv.gz     read label(PAT/MAT/none) by our own rule on
#                                yak's counts: PAT or MAT if that parent holds
#                                at least KMER_MIN_FRAC (0.8) of the read's
#                                parent-specific k-mers (pat + mat) and at
#                                least KMER_MIN_COUNT (3) of them; READ_TRUTH=kmer
#
# Kept when newer than the BAM and the k-mer databases; the kmer table when
# newer than the yak table.
#
#   PHASING_DEPTHS=30 bash bin/yak_truth.sh
#   SAMPLE=HG005 bash bin/yak_truth.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

YAK_THREADS="${YAK_THREADS:-16}"
KMER_MIN_FRAC="${KMER_MIN_FRAC:-0.8}"
KMER_MIN_COUNT="${KMER_MIN_COUNT:-3}"
export PATH="${TRIO_BIN}:${PATH}"
require_cmds yak
require_file "${TRIO_DIR}/pat.yak" "paternal k-mers - run trio_kmers.sh"
require_file "${TRIO_DIR}/mat.yak" "maternal k-mers - run trio_kmers.sh"
mkdir -p "${READS_DIR}" "${LOGS_DIR}/reads"

# kmer_labels <yak table> <out> <cell> - PAT/MAT/none by the share of the
# read's parent-specific k-mers. Columns 4 and 5 of the yak table are taken as
# yak's paternal and maternal counts; that is checked against yak's own p/m
# calls before anything is written. A change of KMER_MIN_* needs the kmer
# tables removed, they are kept while newer than the yak table.
kmer_labels() {
  local yak="$1" out="$2" cell="$3" agree
  agree="$(gzip -dc "${yak}" | awk -F'\t' '!/^#/ && ($3 == "p" || $3 == "m") {
      n++; if (($3 == "p" && $4 + 0 > $5 + 0) || ($3 == "m" && $5 + 0 > $4 + 0)) a++ }
      END { printf "%.4f", n ? a / n : 0 }')"
  awk -v a="${agree}" 'BEGIN { exit !(a >= 0.95) }' \
    || die "${cell}: columns 4-5 of ${yak} do not look like yak's paternal/maternal counts (they agree with its p/m calls on ${agree} of the reads); first line: $(gzip -dc "${yak}" | awk '!/^#/ { print; exit }')"
  gzip -dc "${yak}" \
    | awk -F'\t' -v OFS='\t' -v f="${KMER_MIN_FRAC}" -v c="${KMER_MIN_COUNT}" '
        BEGIN { print "#read", "label", "pat", "mat" }
        /^#/ { next }
        { p = $4 + 0; m = $5 + 0; l = "none"
          if (p >= c && p >= f * (p + m)) l = "PAT"
          else if (m >= c && m >= f * (p + m)) l = "MAT"
          print $1, l, p, m }' \
    | "${BGZIP}" -c >"${out}.tmp" \
    || die "kmer labels failed for ${cell}"
  mv -f "${out}.tmp" "${out}"
  log "${cell}: kmer labels (>= ${KMER_MIN_FRAC} of the parent-specific k-mers, >= ${KMER_MIN_COUNT} of them, yak p/m agree ${agree}): $(gzip -dc "${out}" | awk -F'\t' 'NR > 1 { n[$2]++ } END { printf "PAT %d, MAT %d, none %d", n["PAT"], n["MAT"], n["none"] }')"
}

for platform in "${PLATFORMS[@]}"; do
  for depth in "${DEPTHS[@]}"; do
    cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
    bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
    out="${READS_DIR}/${cell}.yak.tsv.gz"
    [[ -s "${bam}" ]] || { warn "no BAM for ${cell}"; continue; }
    kmer="${READS_DIR}/${cell}.kmer.tsv.gz"
    if newer_than "${out}" "${bam}" && newer_than "${out}" "${TRIO_DIR}/pat.yak" \
        && newer_than "${out}" "${TRIO_DIR}/mat.yak"; then
      newer_than "${kmer}" "${out}" && { log "skip ${cell}"; continue; }
      kmer_labels "${out}" "${kmer}" "${cell}"
      continue
    fi

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
    kmer_labels "${out}" "${kmer}" "${cell}"
  done
done
log "yak truth: ${READS_DIR}/*.yak.tsv.gz"
