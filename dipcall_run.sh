#!/usr/bin/env bash
# Assembly -> phased VCF with dipcall, so eval_run.sh can score assemblies
# with the same metrics as the read-based phasers.
#
# phased/<platform>_<ref>_<sample>_<depth>x_<region>_hifiasm_dipcall.vcf.gz
#
# dipcall also writes <prefix>.dip.bed: the regions where both haplotypes have
# exactly one alignment. That is the "phase came from continuous assembly"
# filter - kept alongside the VCF.
#
# Score these with:
#   PHASING_CALLERS=hifiasm PHASING_PHASERS=dipcall bash eval_run.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${PHASE_DIR}" "${LOGS_DIR}/dipcall"

dipcall_one() {
  local platform="$1" depth="$2"
  local asm="${ASM_DIR}/$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")_hifiasm"
  local name; name="$(vcf_name "${platform}" "${depth}" hifiasm dipcall)"
  local out="${PHASE_DIR}/${name}.vcf.gz"
  local work="${PHASE_DIR}/${name}.dipcall"
  local logf="${LOGS_DIR}/dipcall/${name}.log"

  is_done "${out}" && { log "skip ${name}"; return 0; }

  require_file "${asm}/hap1.fa" "hap1 for ${name}"
  require_file "${asm}/hap2.fa" "hap2 for ${name}"

  rm -rf "${work}"
  mkdir -p "${work}"

  log "${name}: dipcall against $(basename "${REF}")"

  local par=()
  [[ -n "${DIPCALL_PAR}" ]] && par=(-x "${DIPCALL_PAR}")

  "${DIPCALL_BIN}/run-dipcall" \
      ${par[@]+"${par[@]}"} \
      -t "${DIPCALL_THREADS}" \
      "${work}/dip" \
      "${REF}" \
      "${asm}/hap1.fa" \
      "${asm}/hap2.fa" \
    >"${work}/dip.mak" 2>"${logf}" \
    || die "run-dipcall failed for ${name}, see ${logf}"

  make -j2 -f "${work}/dip.mak" >>"${logf}" 2>&1 \
    || die "dipcall make failed for ${name}, see ${logf}"

  require_file "${work}/dip.dip.vcf.gz" "dipcall VCF for ${name}"

  cp "${work}/dip.dip.vcf.gz" "${out}"
  "${TABIX}" -f -p vcf "${out}"
  cp "${work}/dip.dip.bed" "${PHASE_DIR}/${name}.dip.bed"

  log "${name}: $("${BCFTOOLS}" index -n "${out}") variants, confident regions in ${name}.dip.bed"
}

log "dipcall: ${#ASM_PLATFORMS[@]} platforms x ${#ASM_DEPTHS[@]} depths against ${REFERENCE_NAME}"
for platform in "${ASM_PLATFORMS[@]}"; do
  for depth in "${ASM_DEPTHS[@]}"; do
    dipcall_one "${platform}" "${depth}"
  done
done
log "dipcall done"
