#!/usr/bin/env bash
# hifiasm over the assembly grid, ASM_PLATFORMS x ASM_DEPTHS.
#
#   ASSEMBLER=hifiasm      reads only: hap1/hap2 agree within a contig, not
#                          from one contig to the next
#   ASSEMBLER=hifiasmtrio  phased with the parental k-mers of trio_kmers.sh
#
# asm/<platform>_<ref>_<sample>_<depth>x_wgs_<assembler>/hap1.fa, hap2.fa

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

ASSEMBLER="${ASSEMBLER:-hifiasm}"
case "${ASSEMBLER}" in
  hifiasm)
    trio=()
    ;;
  hifiasmtrio)
    require_file "${TRIO_DIR}/pat.yak" "paternal k-mers - run trio_kmers.sh first"
    require_file "${TRIO_DIR}/mat.yak" "maternal k-mers - run trio_kmers.sh first"
    trio=(-1 "${TRIO_DIR}/pat.yak" -2 "${TRIO_DIR}/mat.yak")
    ;;
  *)
    die "ASSEMBLER must be hifiasm or hifiasmtrio, not '${ASSEMBLER}'"
    ;;
esac

mkdir -p "${ASM_DIR}" "${LOGS_DIR}/asm"

hifiasm_preset() {
  case "$1" in
    ont4k|ont5k) printf -- '--ont\n' ;;
    hifi|revio)  printf '\n' ;;
  esac
}

assemble_one() {
  local platform="$1" depth="$2"
  local out; out="$(asm_dir "${platform}" "${depth}" "${ASSEMBLER}")"
  local name; name="$(basename "${out}")"
  local logf="${LOGS_DIR}/asm/${name}.log"

  [[ -s "${out}/hap1.fa" && -s "${out}/hap2.fa" ]] \
    && { log "skip ${name}"; return 0; }

  local bam; bam="$(sub_bam "${platform}" "${ASM_READS_REF}" "${depth}")"
  require_file "${bam}" "cell BAM for ${name}"
  mkdir -p "${out}"

  local fq="${out}/reads.fq.gz"
  if [[ ! -s "${fq}" ]]; then
    log "${name}: BAM -> FASTQ${ASM_REGION:+ (${ASM_REGION} only)}"
    bam_to_fastq "${bam}" "${fq}"
  fi

  log "${name}: hifiasm $(hifiasm_preset "${platform}")${trio[*]+ with parental k-mers} on ${ASM_THREADS} threads"
  "${HIFIASM}" \
      -o "${out}/asm" \
      -t "${ASM_THREADS}" \
      $(hifiasm_preset "${platform}") \
      ${trio[@]+"${trio[@]}"} \
      "${fq}" \
    >"${logf}" 2>&1 \
    || die "hifiasm failed for ${name}, see ${logf}"

  # reads only: asm.bp.hap1.p_ctg.gfa; with parental k-mers: asm.dip.hap1.p_ctg.gfa
  local h1 h2
  h1="$(ls "${out}"/asm.*.hap1.p_ctg.gfa 2>/dev/null | head -1)"
  h2="$(ls "${out}"/asm.*.hap2.p_ctg.gfa 2>/dev/null | head -1)"
  [[ -n "${h1}" && -n "${h2}" ]] \
    || die "${name}: no hap1/hap2 p_ctg.gfa in ${out}, see ${logf}"

  log "${name}: $(basename "${h1}") and $(basename "${h2}") -> FASTA"
  awk '/^S/ { print ">" $2; print $3 }' "${h1}" >"${out}/hap1.fa"
  awk '/^S/ { print ">" $2; print $3 }' "${h2}" >"${out}/hap2.fa"

  rm -f "${fq}"
  log "${name}: hap1 $(grep -c '^>' "${out}/hap1.fa") contigs, hap2 $(grep -c '^>' "${out}/hap2.fa") contigs"
}

log "${ASSEMBLER}: ${#ASM_PLATFORMS[@]} platforms x ${#ASM_DEPTHS[@]} depths"
for platform in "${ASM_PLATFORMS[@]}"; do
  for depth in "${ASM_DEPTHS[@]}"; do
    assemble_one "${platform}" "${depth}"
  done
done
log "assembly done"
