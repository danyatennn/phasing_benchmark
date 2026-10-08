#!/usr/bin/env bash

# Produces two VCFs over exactly the same sites:
# truth.het_snps.phased.vcf.gz - phase kept
# truth.het_snps.unphased.vcf.gz - phase removed
set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

require_vars SAMPLE BCFTOOLS TABIX MOSDEPTH \
  HG002_GRCH38_V5Q_BED HG002_GRCH38_V5Q_TRUTH_VCF \
  ONT_HG002_V14_5KHZ_NO_ALT_BAM ONT_HG002_dorado_NO_ALT_BAM \
  PB_HG002_20K_NO_ALT_BAM HG002_REVIO_BAM
require_cmds "${BCFTOOLS}" "${TABIX}" "${MOSDEPTH}"
require_file "${HG002_GRCH38_V5Q_BED}" "confident-regions BED"
require_file "${HG002_GRCH38_V5Q_TRUTH_VCF}" "truth VCF"

mkdir -p "${DERIVED_DIR}"
mkdir -p "${DATA_DIR}"

log "selecting heterozygous biallelic SNPs"
"${BCFTOOLS}" view \
    -T "${HG002_GRCH38_V5Q_BED}" \
    -p \
    -m2 -M2 -v snps -g het \
    "${HG002_GRCH38_V5Q_TRUTH_VCF}" \
| "${BCFTOOLS}" annotate -x 'INFO,^FORMAT/GT' -Oz -o "${TRUTH_PHASED}"
"${TABIX}" -f -p vcf "${TRUTH_PHASED}"


# unphased phaser input
log "stripping phase to build the phaser input"
"${BCFTOOLS}" +setGT "${TRUTH_PHASED}" -- -t a -n u \
| "${BCFTOOLS}" view -Oz -o "${TRUTH_UNPHASED}"
"${TABIX}" -f -p vcf "${TRUTH_UNPHASED}"

# validation
n_phased=$("${BCFTOOLS}" index -n "${TRUTH_PHASED}")
n_unphased=$("${BCFTOOLS}" index -n "${TRUTH_UNPHASED}")

[[ "${n_phased}" == "${n_unphased}" ]] \
  || die "site sets differ: ${n_phased} phased vs ${n_unphased} unphased"

mkdir -p "${TMP_DIR}"
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

log "truth ready: ${n_phased} heterozygous biallelic SNPs"
log "  phased truth: ${TRUTH_PHASED}"
log "  unphased (unphased truth vcf input): ${TRUTH_UNPHASED}"

# source depth, one mosdepth per platform, all platforms at once
SRC_DEPTH_TMP="${DATA_DIR}/.source_depth_tmp"
rm -rf "${SRC_DEPTH_TMP}"
mkdir -p "${SRC_DEPTH_TMP}"

declare -A depth_pid=()

for platform in "${PLATFORMS[@]}"; do
  bam="$(source_bam "${platform}")"
  require_file "${bam}" "source BAM for ${platform}"
  require_file "${bam}.bai" "source BAM index for ${platform}"
  log "measuring source depth for ${platform}"
  (
    MOSDEPTH_THREADS="${DEPTH_THREADS}"
    measure_depth "${bam}" >"${SRC_DEPTH_TMP}/${platform}"
  ) &
  depth_pid["${platform}"]=$!
done

failed=()
for platform in "${PLATFORMS[@]}"; do
  if wait "${depth_pid[${platform}]}"; then
    log "  ${platform}: $(<"${SRC_DEPTH_TMP}/${platform}")x"
  else
    warn "  ${platform}: depth measurement failed"
    failed+=("${platform}")
  fi
done

[[ ${#failed[@]} -eq 0 ]] \
  || die "source depth measurement failed for: ${failed[*]}"

printf 'platform\tmean_depth\n' >"${DEPTH_TSV}"
for platform in "${PLATFORMS[@]}"; do
  printf '%s\t%s\n' "${platform}" "$(<"${SRC_DEPTH_TMP}/${platform}")" \
    >>"${DEPTH_TSV}"
done

rm -rf "${SRC_DEPTH_TMP}"
log "source depth table: ${DEPTH_TSV}"