#!/usr/bin/env bash
# The phaser inputs of experiment A, from the truth VCF inside the confident
# BED, phased records only (as bin/make_truth.sh). Both keep the truth het
# SNPs of the main grid exactly (TRUTH_PHASED) and add:
#
#   truthindel  biallelic het indels
#   truthall    every het record: 1/2 sites split into two records
#               (bcftools norm -m -any), MNPs and complex variants atomized
#               into SNPs and indels (-a); a part that both haplotypes share
#               comes out twice and is dropped (plots/truth_variants.py)
#
# Indels are left-aligned and trimmed against REF (bcftools norm -f). Each file
# is sorted once; the phased and the unphased one have the same records in the
# same order, which matters where two records share a position (whatshap takes
# the first). Output:
#   derived/truth[_ref].{indel,all}.phased.vcf.gz     to score the added variants
#   derived/truth[_ref].{indel,all}.unphased.vcf.gz   the phasers' input
#
#   bash bin/make_truth_variants.sh                 HG002, GRCh38
#   REFERENCE_NAME=T2T bash bin/make_truth_variants.sh
#   SAMPLE=HG005 bash bin/make_truth_variants.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

require_file "${TRUTH_PHASED}" "truth het SNPs - run bin/make_truth.sh first"
require_file "${CONFIDENT_BED}" "confident BED of ${SAMPLE} on ${REFERENCE_NAME}"
require_file "${TRUTH_SRC_VCF}" "truth VCF of ${SAMPLE} on ${REFERENCE_NAME}"
mkdir -p "${DERIVED_DIR}" "${TMP_DIR}" "${LOGS_DIR}"
normlog="${LOGS_DIR}/make_truth_variants_${REFERENCE_NAME}.log"

# extra_records <bcftools view filters...> - the source records inside the BED,
# phased, GT only
extra_records() {
  "${BCFTOOLS}" view -T "${CONFIDENT_BED}" -p "$@" "${TRUTH_SRC_VCF}" \
    | "${BCFTOOLS}" annotate -x 'INFO,^FORMAT/GT'
}

# todo <name> <phased> <unphased> - is the input still to be built?
todo() {
  if [[ -z "${FORCE_TRUTH:-}" ]] && is_done "$2" && is_done "$3" && newer_than "$2" "${TRUTH_PHASED}"; then
    log "${REFERENCE_NAME}: $1 already built (FORCE_TRUTH=1 to rebuild)"
    return 1
  fi
}

# build <name> <phased out> <unphased out> - reads the extra records on stdin
build() {
  local name="$1" phased="$2" unphased="$3"
  local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
  cat >"${tmp}/extra.vcf"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/truth_variants.py" \
      --snps "${TRUTH_PHASED}" --extra "${tmp}/extra.vcf" --sample "${SAMPLE}" --out "${tmp}/merged.vcf"
  "${BCFTOOLS}" sort -T "${tmp}" -Oz -o "${phased}.tmp" "${tmp}/merged.vcf"
  "${TABIX}" -f -p vcf "${phased}.tmp"
  "${BCFTOOLS}" +setGT "${phased}.tmp" -- -t a -n u | "${BGZIP}" -c >"${unphased}.tmp"
  "${TABIX}" -f -p vcf "${unphased}.tmp"

  # the SNPs of the main grid, all of them, and no phase left in the input
  local snps leaked
  snps="$("${BCFTOOLS}" isec -n=2 -w1 -c none "${phased}.tmp" "${TRUTH_PHASED}" | grep -vc '^#' || true)"
  [[ "${snps}" -eq "$("${BCFTOOLS}" index -n "${TRUTH_PHASED}")" ]] \
    || die "${name}: ${snps} of the $("${BCFTOOLS}" index -n "${TRUTH_PHASED}") truth SNPs kept"
  leaked="$("${BCFTOOLS}" query -f '[%GT]\n' "${unphased}.tmp" | grep -c '|' || true)"
  [[ "${leaked}" -eq 0 ]] || die "${name}: phase leaked into the input (${leaked} genotypes)"

  for f in "${phased}" "${unphased}"; do
    mv -f "${f}.tmp.tbi" "${f}.tbi"
    mv -f "${f}.tmp" "${f}"
  done
  rm -rf "${tmp}"
  log "${REFERENCE_NAME} ${name}: $("${BCFTOOLS}" index -n "${phased}") records" \
      "($("${BCFTOOLS}" view -H -v indels "${phased}" | wc -l) indels," \
      "$("${BCFTOOLS}" query -f '%CHROM\t%POS\n' "${phased}" | uniq -d | wc -l) positions with two records)"
}

log "${SAMPLE} ${REFERENCE_NAME}: experiment A inputs from $(basename "${TRUTH_SRC_VCF}")"

if todo truthindel "${TRUTH_INDEL_PHASED}" "${TRUTH_INDEL_UNPHASED}"; then
  extra_records -m2 -M2 -v indels -g het \
    | "${BCFTOOLS}" norm -f "${REF}" -c w 2>>"${normlog}" \
    | build truthindel "${TRUTH_INDEL_PHASED}" "${TRUTH_INDEL_UNPHASED}"
fi

if todo truthall "${TRUTH_ALL_PHASED}" "${TRUTH_ALL_UNPHASED}"; then
  extra_records -g het \
    | "${BCFTOOLS}" norm -m -any 2>>"${normlog}" \
    | "${BCFTOOLS}" norm -a --atom-overlaps . 2>>"${normlog}" \
    | "${BCFTOOLS}" norm -f "${REF}" -c w 2>>"${normlog}" \
    | build truthall "${TRUTH_ALL_PHASED}" "${TRUTH_ALL_UNPHASED}"
fi
