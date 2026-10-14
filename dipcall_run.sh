#!/usr/bin/env bash
# Assemblies -> phased VCFs that eval_run.sh scores like any read-based phaser.
#
#   phased/<platform>_<ref>_<sample>_<depth>x_<region>_<assembler>_dipcall.vcf.gz
#
# dipcall aligns hap1 and hap2 to the reference and calls where they differ,
# GT hap1|hap2. It writes no phase sets, and whatshap would read a whole
# chromosome as one block. A phase set has to be what the assembler actually
# vouches for, so switch errors are counted under the same rule as for a
# read-based phaser's blocks:
#   hifiasm      reads alone: past a contig end hap1/hap2 are assigned at
#                random, so a set is a run of SNPs on one pair of contigs
#   *trio        hap1 is paternal on every contig, so a set is a chromosome
#
# Kept, so that the VCF is comparable with the read-based ones:
#   biallelic het SNPs, GT only
#   autosomes - no PAR file is needed, and hap1/hap2 mean nothing on X and Y
#     of a male sample when the assembly was not phased with the parents
#   confident BED x dip.bed - dip.bed is where each haplotype has exactly one
#     long alignment, i.e. where the assembly can be trusted at all
#
# Score with:
#   PHASING_CALLERS="hifiasm hifiasmtrio verkkotrio" PHASING_PHASERS=dipcall \
#   PHASING_DEPTHS="20 30 50" bash eval_run.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${PHASE_DIR}" "${LOGS_DIR}/dipcall"

dipcall_one() {
  local assembler="$1" platform="$2" depth="$3"
  local asm; asm="$(asm_dir "${platform}" "${depth}" "${assembler}")"
  local name; name="$(vcf_name "${platform}" "${depth}" "${assembler}" dipcall)"
  local out="${PHASE_DIR}/${name}.vcf.gz"
  local work="${PHASE_DIR}/${name}.dipcall"
  local pre="${work}/dip"
  local logf="${LOGS_DIR}/dipcall/${name}.log"

  is_done "${out}" && { log "skip ${name}"; return 0; }
  [[ -s "${asm}/hap1.fa" && -s "${asm}/hap2.fa" ]] \
    || { log "skip ${name}: no assembly yet"; return 0; }

  # dipcall itself is the slow part; a finished run is reused
  if [[ ! -s "${pre}.dip.vcf.gz" || ! -s "${pre}.dip.bed" ]]; then
    rm -rf "${work}"
    mkdir -p "${work}"
    log "${name}: dipcall against $(basename "${REF}")"
    "${DIPCALL_BIN}/run-dipcall" -t "${DIPCALL_THREADS}" \
        "${pre}" "${REF}" "${asm}/hap1.fa" "${asm}/hap2.fa" \
      >"${work}/dip.mak" 2>"${logf}" \
      || die "run-dipcall failed for ${name}, see ${logf}"
    make -j2 -f "${work}/dip.mak" >>"${logf}" 2>&1 \
      || die "dipcall make failed for ${name}, see ${logf}"
  fi

  # where the comparison happens: autosomes, confident BED, dip.bed
  grep -E '^chr([1-9]|1[0-9]|2[0-2])[[:space:]]' "${pre}.dip.bed" \
    | sort -k1,1 -k2,2n >"${pre}.blk.dip.bed"
  bed_intersect "${CONFIDENT_BED}" "${pre}.blk.dip.bed" >"${pre}.blk.region.bed"

  "${BCFTOOLS}" view -T "${pre}.blk.region.bed" -m2 -M2 -v snps -g het "${pre}.dip.vcf.gz" \
    | "${BCFTOOLS}" annotate -x 'INFO,^FORMAT/GT' >"${pre}.blk.snps.vcf"

  # which contig each haplotype has where: its long alignments (>= 50 kb,
  # mapQ >= 5), the same ones dipcall builds dip.bed from
  local h
  for h in 1 2; do
    gzip -dc "${pre}.hap${h}.paf.gz" \
      | awk -F'\t' -v OFS='\t' '$12 >= 5 && $9 - $8 >= 50000 { print $6, $8, $9, $1 }' \
      | sort -k1,1 -k2,2n >"${pre}.blk.hap${h}.bed"
  done

  local span=""
  [[ "${assembler}" == *trio ]] && span=chromosome
  assembly_phase_sets "${pre}.blk.hap1.bed" "${pre}.blk.hap2.bed" "${pre}.blk.snps.vcf" ${span} \
    >"${pre}.blk.ps.tsv"
  [[ -s "${pre}.blk.ps.tsv" ]] || die "${name}: no SNP kept, see ${work}"

  add_phase_sets "${pre}.blk.ps.tsv" "${pre}.blk.snps.vcf" "${SAMPLE}" \
    | "${BGZIP}" -c >"${out}.tmp"
  "${TABIX}" -f -p vcf "${out}.tmp"
  mv -f "${out}.tmp.tbi" "${out}.tbi"
  mv -f "${out}.tmp" "${out}"
  cp "${pre}.dip.bed" "${PHASE_DIR}/${name}.dip.bed"

  log "${name}: $("${BCFTOOLS}" index -n "${out}") het SNPs in $(cut -f3 "${pre}.blk.ps.tsv" | sort -u | wc -l) phase sets"
}

log "dipcall: ${ASM_CALLERS[*]} x ${ASM_PLATFORMS[*]} x ${ASM_DEPTHS[*]} against ${REFERENCE_NAME}"
for assembler in "${ASM_CALLERS[@]}"; do
  for platform in "${ASM_PLATFORMS[@]}"; do
    for depth in "${ASM_DEPTHS[@]}"; do
      dipcall_one "${assembler}" "${platform}" "${depth}"
    done
  done
done
log "dipcall done"
