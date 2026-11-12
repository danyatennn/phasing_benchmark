#!/usr/bin/env bash
# WhatsHap haplotag against LongPhase haplotag on one chromosome: the phased VCF
# of each 30x run (truth input by default) and the cell's BAM, both cut to the
# chromosome; LongPhase haplotag with its defaults. WhatsHap's tags are the
# grid's own (metrics/<run>.haplotag.tsv, eval_run.sh). How the reads differ:
# plots/haplotag_compare.py.
#
#   bash bin/haplotag_compare.sh                         GRCh38, chr20
#   HAPLOTAG_CHROM=chr1 REFERENCE_NAME=T2T bash bin/haplotag_compare.sh
#   PHASING_CALLERS="truth clair3" PHASING_PLATFORMS=hifi bash bin/haplotag_compare.sh
#
# -> metrics_haplotag/haplotag_compare_<ref>_<chrom>.tsv

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

CHROM="${HAPLOTAG_CHROM:-chr20}"
DEPTH=30
CMP_CALLERS="${PHASING_CALLERS:-truth}"
OUT_DIR="${SAMPLE_ROOT}/metrics_haplotag"
out="${OUT_DIR}/haplotag_compare_${REFERENCE_NAME}_${CHROM}.tsv"
mkdir -p "${OUT_DIR}/work"
rm -f "${out}"

for platform in "${PLATFORMS[@]}"; do
  w="${OUT_DIR}/work/$(bam_key "${platform}" "${REFERENCE_NAME}" "${DEPTH}")_${CHROM}"
  mkdir -p "${w}"
  bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${DEPTH}")"
  [[ -s "${w}/in.bam.bai" ]] || { "${SAMTOOLS}" view -b -o "${w}/in.bam" "${bam}" "${CHROM}" && "${SAMTOOLS}" index "${w}/in.bam"; }
  "${SAMTOOLS}" view -F 0x904 "${w}/in.bam" | awk -v OFS='\t' '{ print $1, $5, length($10) }' >"${w}/reads.tsv"
  for caller in ${CMP_CALLERS}; do
    for phaser in "${PHASERS[@]}"; do
      phaser_applies "${phaser}" "${platform}" || continue
      name="$(vcf_name "${platform}" "${DEPTH}" "${caller}" "${phaser}")"
      vcf="${PHASE_DIR}/${name}.vcf.gz"
      tags="${METRICS_DIR}/${name}.haplotag.tsv"
      [[ -s "${vcf}" && -s "${tags}" ]] || { warn "${name}: no phased VCF or WhatsHap tags, skipped"; continue; }
      log "${name} ${CHROM}: longphase haplotag"
      "${BCFTOOLS}" view -r "${CHROM}" -Oz -o "${w}/${phaser}_${caller}.vcf.gz" "${vcf}"
      "${TABIX}" -f -p vcf "${w}/${phaser}_${caller}.vcf.gz"
      "${LONGPHASE}" haplotag -r "${REF}" -s "${w}/${phaser}_${caller}.vcf.gz" -b "${w}/in.bam" \
          -t "${PHASE_THREADS}" -o "${w}/lp_${phaser}_${caller}" >"${w}/lp_${phaser}_${caller}.log" 2>&1 \
        || die "longphase haplotag failed, see ${w}/lp_${phaser}_${caller}.log"
      "${SAMTOOLS}" view -F 0x904 "${w}/lp_${phaser}_${caller}.bam" \
        | awk -v OFS='\t' '{ hp = ""; ps = ""
            for (i = 12; i <= NF; i++) { if ($i ~ /^HP:i:/) hp = substr($i, 6); if ($i ~ /^PS:i:/) ps = substr($i, 6) }
            print $1, hp, ps }' >"${w}/lp_${phaser}_${caller}.tsv"
      "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/haplotag_compare.py" --reads "${w}/reads.tsv" \
          --whatshap "${tags}" --longphase "${w}/lp_${phaser}_${caller}.tsv" --chrom "${CHROM}" \
          --label "${platform}_${caller}_${phaser}" --out "${out}"
      rm -f "${w}/lp_${phaser}_${caller}.bam"*
    done
  done
done
log "-> ${out}"
