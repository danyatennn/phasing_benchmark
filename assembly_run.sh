#!/usr/bin/env bash

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${ASM_DIR}" "${LOGS_DIR}/asm"

hifiasm_preset() {
  case "$1" in
    ont4k|ont5k) printf -- '--ont\n' ;;
    hifi|revio)  printf '\n' ;;
  esac
}

assemble_one() {
  local platform="$1" depth="$2"
  local name; name="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")_hifiasm"
  local out="${ASM_DIR}/${name}"
  local logf="${LOGS_DIR}/asm/${name}.log"

  [[ -s "${out}/hap1.fa" && -s "${out}/hap2.fa" ]] \
    && { log "skip ${name}"; return 0; }

  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  require_file "${bam}" "cell BAM for ${name}"
  mkdir -p "${out}"

  local fq="${out}/reads.fq.gz"
  if [[ ! -s "${fq}" ]]; then
    log "${name}: BAM -> FASTQ${ASM_REGION:+ (${ASM_REGION} only)}"
    "${SAMTOOLS}" view -u -@ "${SORT_THREADS}" "${bam}" ${ASM_REGION} \
      | "${SAMTOOLS}" fastq -@ "${SORT_THREADS}" -n - \
      | "${BGZIP}" -@ "${SORT_THREADS}" -c >"${fq}.tmp"
    mv -f "${fq}.tmp" "${fq}"
  fi

  log "${name}: hifiasm $(hifiasm_preset "${platform}") on ${ASM_THREADS} threads"
  "${HIFIASM}" \
      -o "${out}/asm" \
      -t "${ASM_THREADS}" \
      $(hifiasm_preset "${platform}") \
      "${fq}" \
    >"${logf}" 2>&1 \
    || die "hifiasm failed for ${name}, see ${logf}"

  local h1 h2
  h1="$(ls "${out}"/*hap1*p_ctg.gfa 2>/dev/null | head -1)"
  h2="$(ls "${out}"/*hap2*p_ctg.gfa 2>/dev/null | head -1)"
  [[ -n "${h1}" && -n "${h2}" ]] \
    || die "${name}: no hap1/hap2 p_ctg.gfa in ${out}, see ${logf}"

  log "${name}: $(basename "${h1}") and $(basename "${h2}") -> FASTA"
  awk '/^S/ { print ">" $2; print $3 }' "${h1}" >"${out}/hap1.fa"
  awk '/^S/ { print ">" $2; print $3 }' "${h2}" >"${out}/hap2.fa"

  rm -f "${fq}"
  log "${name}: hap1 $(grep -c '^>' "${out}/hap1.fa") contigs, hap2 $(grep -c '^>' "${out}/hap2.fa") contigs"
}

log "assembling ${#ASM_PLATFORMS[@]} platforms x ${#ASM_DEPTHS[@]} depths"
for platform in "${ASM_PLATFORMS[@]}"; do
  for depth in "${ASM_DEPTHS[@]}"; do
    assemble_one "${platform}" "${depth}"
  done
done
log "assembly done"
