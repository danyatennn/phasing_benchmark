#!/usr/bin/env bash
# Haplotype-resolved assemblies with hifiasm, one per cell.
# asm/<platform>_<ref>_<sample>_<depth>x_hifiasm/hap1.fa  and  hap2.fa
#
# Assembly is reference-free, so the reads come straight out of the cell BAM.
# hifiasm splits haplotypes on graph bubbles: phase is continuous inside a
# contig, but which contig is "hap1" is arbitrary. dipcall handles that by
# setting PS per contig.
#
# Verkko is not here: it only phases with parental Meryl k-mers, Hi-C or
# Pore-C, none of which we have.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${ASM_DIR}" "${LOGS_DIR}/asm"

# hifiasm preset: ONT simplex needs --ont, HiFi is the default
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
    log "${name}: BAM -> FASTQ"
    "${SAMTOOLS}" fastq -@ "${SORT_THREADS}" -n "${bam}" \
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

  # output naming moved across hifiasm versions, so find the two haplotypes
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
