#!/usr/bin/env bash
# Clair3 over the grid
# vcf/<platform>_<ref>_<sample>_<depth>x_<region>_clair3.vcf.gz - all calls
# vcf/<platform>_<ref>_<sample>_<depth>x_<region>_clair3.het_snps.vcf.gz - phaser input

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${VCF_DIR}" "${CLAIR3_DIR}" "${LOGS_DIR}/clair3"

call_one() {
  local platform="$1" depth="$2"
  local name; name="$(vcf_name "${platform}" "${depth}" clair3)"
  local out="${VCF_DIR}/${name}.vcf.gz"
  local het="${VCF_DIR}/${name}.het_snps.vcf.gz"

  is_done "${het}" && { log "skip ${name}"; return 0; }

  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local work="${CLAIR3_DIR}/${name}"
  local logf="${LOGS_DIR}/clair3/${name}.log"

  rm -rf "${work}"
  mkdir -p "${work}"

  log "${name}: $(clair3_model "${platform}") on GPU"

  "${SING}" exec --nv --cleanenv --env TMPDIR=/tmp \
      -B "${DATA_DIR},$(dirname "${REF}"),$(dirname "${CONFIDENT_BED}"),${CLAIR3_DIR}" \
      "${CLAIR3_SIF}" \
      /opt/bin/run_clair3.sh \
        --bam_fn="${bam}" \
        --ref_fn="${REF}" \
        --threads="${CLAIR3_THREADS}" \
        --platform="$(clair3_platform "${platform}")" \
        --model_path="/opt/models/$(clair3_model "${platform}")" \
        --sample_name="${SAMPLE}" \
        --bed_fn="${CONFIDENT_BED}" \
        --output="${work}" \
        --use_gpu \
      >"${logf}" 2>&1 \
    || die "clair3 failed for ${name}, see ${logf}"

  cp "${work}/merge_output.vcf.gz" "${out}"
  "${TABIX}" -f -p vcf "${out}"

  "${BCFTOOLS}" view -f PASS -m2 -M2 -v snps -g het -Oz -o "${het}" "${out}"
  "${TABIX}" -f -p vcf "${het}"

  log "${name}: $("${BCFTOOLS}" index -n "${het}") het SNPs"
}

log "clair3: ${#PLATFORMS[@]} platforms x ${#DEPTHS[@]} depths on ${REFERENCE_NAME}"
for_each_bam call_one
log "clair3 done"
