set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

if [[ -z "${FORCE_TRUTH:-}" ]] && is_done "${TRUTH_PHASED}" && is_done "${TRUTH_UNPHASED}"; then
  log "${REFERENCE_NAME}: truth already built (FORCE_TRUTH=1 to rebuild)"
  exit 0
fi

require_file "${CONFIDENT_BED}" "confident-regions BED for ${REFERENCE_NAME}"
require_file "${TRUTH_SRC_VCF}" "truth VCF for ${REFERENCE_NAME}"

mkdir -p "${DERIVED_DIR}" "${TMP_DIR}"

log "${REFERENCE_NAME}: heterozygous biallelic SNPs from $(basename "${TRUTH_SRC_VCF}")"
"${BCFTOOLS}" view \
    -T "${CONFIDENT_BED}" \
    -p \
    -m2 -M2 -v snps -g het \
    "${TRUTH_SRC_VCF}" \
| "${BCFTOOLS}" annotate -x 'INFO,^FORMAT/GT' -Oz -o "${TRUTH_PHASED}"
"${TABIX}" -f -p vcf "${TRUTH_PHASED}"

if [[ "$("${BCFTOOLS}" query -l "${TRUTH_PHASED}")" != "${SAMPLE}" ]]; then
  log "renaming sample $("${BCFTOOLS}" query -l "${TRUTH_PHASED}") -> ${SAMPLE}"
  printf '%s\n' "${SAMPLE}" >"${DERIVED_DIR}/sample${_ref_sfx}.txt"
  "${BCFTOOLS}" reheader -s "${DERIVED_DIR}/sample${_ref_sfx}.txt" \
      -o "${TRUTH_PHASED}.tmp" "${TRUTH_PHASED}"
  mv -f "${TRUTH_PHASED}.tmp" "${TRUTH_PHASED}"
  "${TABIX}" -f -p vcf "${TRUTH_PHASED}"
fi

log "stripping phase to build the phaser input"
"${BCFTOOLS}" +setGT "${TRUTH_PHASED}" -- -t a -n u \
| "${BCFTOOLS}" view -Oz -o "${TRUTH_UNPHASED}"
"${TABIX}" -f -p vcf "${TRUTH_UNPHASED}"

# validation
n_phased=$("${BCFTOOLS}" index -n "${TRUTH_PHASED}")
n_unphased=$("${BCFTOOLS}" index -n "${TRUTH_UNPHASED}")
[[ "${n_phased}" == "${n_unphased}" ]] \
  || die "site sets differ: ${n_phased} phased vs ${n_unphased} unphased"

gt_tmp="$(mktemp -p "${TMP_DIR}")"
trap 'rm -f "${gt_tmp}"' EXIT

"${BCFTOOLS}" query -f '[%GT]\n' "${TRUTH_UNPHASED}" >"${gt_tmp}"
leaked=$(grep -c '|' "${gt_tmp}" || true)
[[ "${leaked}" -eq 0 ]] \
  || die "phase leaked into the phaser input: ${leaked} genotypes still contain '|'"

"${BCFTOOLS}" query -f '[%GT]\n' "${TRUTH_PHASED}" >"${gt_tmp}"
kept=$(grep -c '|' "${gt_tmp}" || true)
[[ "${kept}" -eq "${n_phased}" ]] \
  || die "evaluation reference is not fully phased: ${kept}/${n_phased}"

log "${REFERENCE_NAME} truth ready: ${n_phased} heterozygous biallelic SNPs"
log "  phased:   ${TRUTH_PHASED}"
log "  unphased: ${TRUTH_UNPHASED}"
