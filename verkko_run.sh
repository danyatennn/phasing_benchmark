#!/usr/bin/env bash
# verkko phased with the parental k-mers of trio_kmers.sh, over the assembly
# grid. PacBio cells only: verkko needs HiFi-grade reads, and our ONT reads are
# raw simplex. The reads are the cell's alone, as for every other method.
#
# asm/<platform>_<ref>_<sample>_<depth>x_wgs_verkkotrio/hap1.fa, hap2.fa
#
# verkko picks up where it stopped when rerun on the same directory.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

for db in "${TRIO_DIR}/pat.k30.only.meryl" "${TRIO_DIR}/mat.k30.only.meryl"; do
  [[ -d "${db}" ]] || die "missing ${db} - run trio_kmers.sh first"
done
[[ -x "${VERKKO_BIN}/verkko" ]] || die "verkko not found in ${VERKKO_BIN}"
export PATH="${VERKKO_BIN}:${PATH}"

mkdir -p "${ASM_DIR}" "${LOGS_DIR}/asm"

verkko_one() {
  local platform="$1" depth="$2"
  local out; out="$(asm_dir "${platform}" "${depth}" verkkotrio)"
  local name; name="$(basename "${out}")"
  local logf="${LOGS_DIR}/asm/${name}.log"

  case "${platform}" in
    hifi|revio) ;;
    *) log "skip ${name}: verkko needs HiFi-grade reads"; return 0 ;;
  esac

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

  log "${name}: verkko on ${VERKKO_THREADS} CPUs, ${VERKKO_MEMORY} GB"
  "${VERKKO_BIN}/verkko" \
      -d "${out}/verkko" \
      --hifi "${fq}" \
      --hap-kmers "${TRIO_DIR}/pat.k30.only.meryl" "${TRIO_DIR}/mat.k30.only.meryl" trio \
      --local-cpus "${VERKKO_THREADS}" \
      --local-memory "${VERKKO_MEMORY}" \
    >>"${logf}" 2>&1 \
    || die "verkko failed for ${name}, see ${logf}"

  local h
  for h in 1 2; do
    require_file "${out}/verkko/assembly.haplotype${h}.fasta" "verkko haplotype ${h} for ${name}"
    cp "${out}/verkko/assembly.haplotype${h}.fasta" "${out}/hap${h}.fa"
  done

  rm -f "${fq}"
  log "${name}: hap1 $(grep -c '^>' "${out}/hap1.fa") contigs, hap2 $(grep -c '^>' "${out}/hap2.fa") contigs"
}

log "verkkotrio: ${#ASM_PLATFORMS[@]} platforms x ${#ASM_DEPTHS[@]} depths"
for platform in "${ASM_PLATFORMS[@]}"; do
  for depth in "${ASM_DEPTHS[@]}"; do
    verkko_one "${platform}" "${depth}"
  done
done
log "verkko done"
