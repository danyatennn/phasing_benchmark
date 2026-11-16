#!/usr/bin/env bash
# CMRG with the files as GIAB released them, no confident BED anywhere:
#   GRCh38   HG002_GRCh38_CMRG_smallvar_v1.00.{vcf.gz,bed}, our 30x GRCh38 BAMs
#   CHM13v1  HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.{vcf.gz,bed}, the reads
#            of our 30x T2T-CHM13 v2.0 BAMs aligned again to T2T-CHM13 v1.0
#            (the same reads, minimap2 and presets as for v2.0)
# Phaser input, CMRG_INPUT (het SNPs of chr1-22, phase removed, no confident BED):
#   cmrg   (default) the CMRG VCF itself inside the CMRG BED, on both references alike
#   truth  GRCh38: GIAB v5.0q, all of it; CHM13v1: v5.0q has no v1.0 release,
#          so the Q100 assembly (what v5.0q is built from) by dipcall against
#          v1.0, inside dip.bed
# Scored per gene against the CMRG VCF (plots/cmrg_eval.py): phased, in one
# block, correct end to end, switches, Hamming; next to each gene the share of
# it in the confident BED of GRCh38 and of T2T-CHM13 (regions/genes/CMRG/genes.tsv).
#
#   bash bin/cmrg_dev.sh ref                                download v1.0, index, compare with v2.0
#   PHASING_PLATFORMS="ont4k ont5k" bash bin/cmrg_dev.sh align    (CHM13v1 only)
#   CMRG_REF=CHM13v1 bash bin/cmrg_dev.sh truth            (minutes; CMRG_INPUT=truth on CHM13v1: dipcall, hours)
#   CMRG_REF=GRCh38 CMRG_JOBS=6 bash bin/cmrg_dev.sh phase truth input, every platform x phaser
#   CMRG_REF=GRCh38 bash bin/cmrg_dev.sh eval
#   CMRG_REF=GRCh38 bash bin/cmrg_dev.sh tags              reads of the CMRG regions: SNPs covered, tags (benchmark, every run)
#   bash bin/cmrg_dev.sh reads                             read recall, consensus labels, reads anchored on v1.0 (after tags on both)
#
# Written to cmrg_dev/<ref>/bam/ (CHM13v1) and cmrg_dev/<ref>/<input>/: truth/, phased/, logs/, eval/.

set -euo pipefail

STEP="${1:?step: ref, align, truth, phase, eval, tags or reads}"
CMRG_REF="${CMRG_REF:-GRCh38}"
CMRG_INPUT="${CMRG_INPUT:-cmrg}"
[[ "${CMRG_INPUT}" == cmrg || "${CMRG_INPUT}" == truth ]] || { echo "CMRG_INPUT must be cmrg or truth" >&2; exit 1; }
case "${CMRG_REF}" in
  GRCh38)  export REFERENCE_NAME=GRCh38 ;;
  CHM13v1) export REFERENCE_NAME=T2T ;;      # config.sh of T2T: its 30x BAMs are the reads to align
  *) echo "CMRG_REF must be GRCh38 or CHM13v1" >&2; exit 1 ;;
esac

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

[[ "${SAMPLE}" == HG002 ]] || die "CMRG is HG002 only"
CMRG_SRC="/autofs/bal33/zxzheng/data/cmrg"
V1_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.0.fasta.gz"
V1_FA="${PROJECT_ROOT}/ref/chm13.draft_v1.0.fa"
GENES_TSV="${PROJECT_ROOT}/regions/genes/CMRG/genes.tsv"
GENES_GRCH38="${PROJECT_ROOT}/data/annotation/cmrg/GRCh38_CMRG_benchmark_gene_coordinates.bed"
# GIAB, CHM13v1.0/SupplementaryFiles: the GRCh38 equivalent of each row of the
# v1.0 CMRG BED, INTERVAL_<n> = its row n (lengths agree for 2615 of 2616 rows)
GENES_EQUIV="${PROJECT_ROOT}/data/annotation/cmrg/HG002_CHM13_CMRG_smallvar_v1.00_GRCh38-equiv-regions_draft.bed"
OUT="${PROJECT_ROOT}/cmrg_dev"
B="${OUT}/${CMRG_REF}/bam"
D="${OUT}/${CMRG_REF}/${CMRG_INPUT}"
DEPTH=30
if [[ "${CMRG_REF}" == GRCh38 ]]; then
  CREF="${REF_GRCH38}"
  CVCF="${CMRG_SRC}/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz"
  CBED="${CMRG_SRC}/HG002_GRCh38_CMRG_smallvar_v1.00.bed"
else
  CREF="${V1_FA}"
  CVCF="${CMRG_SRC}/HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.vcf.gz"
  CBED="${CMRG_SRC}/HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.bed"
fi
TRUTH_P="${D}/truth/truth.het_snps.phased.vcf.gz"
TRUTH_U="${D}/truth/truth.het_snps.unphased.vcf.gz"
mkdir -p "${B}" "${D}"/{truth,phased,logs/phase,eval}
AUTO_RE='^chr([1-9]|1[0-9]|2[0-2])[[:space:]]'

case "${STEP}" in

ref)
  if [[ ! -s "${V1_FA}.fai" ]]; then
    log "downloading ${V1_URL}"
    curl -fL "${V1_URL}" | gzip -dc >"${V1_FA}.tmp"
    mv -f "${V1_FA}.tmp" "${V1_FA}"
    "${SAMTOOLS}" faidx "${V1_FA}"
  fi
  # which chromosomes of v1.0 are the very sequence of v2.0 (same coordinates)
  log "v1.0 against v2.0, chromosome by chromosome (sequence md5)"
  for c in $(cut -f1 "${V1_FA}.fai"); do
    a="$("${SAMTOOLS}" faidx "${V1_FA}" "${c}" | tail -n +2 | tr -d '\n' | md5sum | cut -c1-12)"
    b="$("${SAMTOOLS}" faidx "${REF_CHM13}" "${c}" 2>/dev/null | tail -n +2 | tr -d '\n' | md5sum | cut -c1-12)"
    printf '%s\t%s\t%s\n' "${c}" "$(awk -v c="${c}" '$1 == c { print $2 }' "${V1_FA}.fai")" \
      "$([[ "${a}" == "${b}" ]] && echo same || echo differs)"
  done | tee "${OUT}/v1_vs_v2.tsv"
  ;;

align)
  [[ "${CMRG_REF}" == CHM13v1 ]] || die "align is for CMRG_REF=CHM13v1 (GRCh38 uses the grid's BAMs)"
  require_file "${V1_FA}.fai" "T2T-CHM13 v1.0 (run the ref step)"
  for p in "${PLATFORMS[@]}"; do
    out="${B}/${p}.bam"
    [[ -s "${out}.bai" ]] && { log "${p}: ${out} exists"; continue; }
    src="$(sub_bam "${p}" T2T "${DEPTH}")"
    require_file "${src}" "30x T2T BAM of ${p}"
    case "${p}" in ont4k|ont5k) preset=lr:hq ;; *) preset=map-hifi ;; esac
    mmi="${B}/chm13v1.${preset//:/_}.mmi"
    [[ -s "${mmi}" ]] || { "${MINIMAP2}" -t "${ALIGN_THREADS}" -x "${preset}" -d "${mmi}.tmp" "${V1_FA}" && mv -f "${mmi}.tmp" "${mmi}"; }
    rg="$("${SAMTOOLS}" view -H "${src}" | grep -m1 '^@RG' | sed 's/\t/\\t/g')"
    log "${p}: $(basename "${src}") -> T2T-CHM13 v1.0 (${preset})"
    "${SAMTOOLS}" fastq -@ "${SORT_THREADS}" -n -T MM,ML "${src}" \
      | "${MINIMAP2}" -t "${ALIGN_THREADS}" -ax "${preset}" -y ${rg:+-R "${rg}"} "${mmi}" - \
      | "${SAMTOOLS}" sort -@ "${SORT_THREADS}" -m "${SORT_MEM}" -o "${out}.tmp" -
    "${SAMTOOLS}" index -@ "${SORT_THREADS}" "${out}.tmp"
    mv -f "${out}.tmp" "${out}"
    mv -f "${out}.tmp.bai" "${out}.bai"
  done
  ;;

truth)
  w="${D}/truth/work"
  mkdir -p "${w}"
  if [[ "${CMRG_INPUT}" == cmrg ]]; then
    # inside the CMRG BED: the v1.0 draft VCF covers the whole genome, the
    # GRCh38 one the CMRG regions only, so both are cut to the benchmark BED
    log "${CMRG_REF}: het SNPs of $(basename "${CVCF}") inside $(basename "${CBED}"), chr1-22"
    grep -E "${AUTO_RE}" "${CBED}" >"${w}/cmrg.auto.bed"
    "${BCFTOOLS}" view -T "${w}/cmrg.auto.bed" -m2 -M2 -v snps -g het "${CVCF}" \
      | "${BCFTOOLS}" annotate -x 'INFO,^FORMAT/GT' -Oz -o "${w}/raw.vcf.gz"
  elif [[ "${CMRG_REF}" == GRCh38 ]]; then
    log "GRCh38: phased het SNPs of $(basename "${TRUTH_SRC_VCF}"), chr1-22, no BED"
    "${BCFTOOLS}" view -t "$(printf 'chr%s,' {1..22} | sed 's/,$//')" -p -m2 -M2 -v snps -g het "${TRUTH_SRC_VCF}" \
      | "${BCFTOOLS}" annotate -x 'INFO,^FORMAT/GT' -Oz -o "${w}/raw.vcf.gz"
  else
    require_file "${V1_FA}.fai" "T2T-CHM13 v1.0 (run the ref step)"
    if [[ ! -s "${w}/dip.dip.vcf.gz" ]]; then
      log "CHM13v1: Q100 paternal and maternal copies -> dipcall against v1.0"
      gzip -dc "${Q100_FASTA}" | awk -v p="${w}/pat.fa" -v m="${w}/mat.fa" \
        '/^>/ { out = ($1 ~ /_PATERNAL$/) ? p : m } { print > out }'
      "${DIPCALL_BIN}/run-dipcall" -t "${DIPCALL_THREADS}" "${w}/dip" "${V1_FA}" "${w}/pat.fa" "${w}/mat.fa" >"${w}/dip.mak"
      make -j2 -f "${w}/dip.mak" >"${D}/logs/dipcall.log" 2>&1 || die "dipcall failed, see ${D}/logs/dipcall.log"
    fi
    grep -E "${AUTO_RE}" "${w}/dip.dip.bed" >"${w}/dip.auto.bed"
    "${BCFTOOLS}" view -T "${w}/dip.auto.bed" -m2 -M2 -v snps -g het "${w}/dip.dip.vcf.gz" \
      | "${BCFTOOLS}" annotate -x 'INFO,^FORMAT/GT' \
      | awk -F'\t' -v OFS='\t' '/^#/ { print; next } { $4 = toupper($4); $5 = toupper($5); print }' \
      | "${BGZIP}" -c >"${w}/raw.vcf.gz"
  fi
  printf '%s\n' "${SAMPLE}" >"${w}/sample.txt"
  "${BCFTOOLS}" reheader -s "${w}/sample.txt" -o "${TRUTH_P}" "${w}/raw.vcf.gz"
  "${TABIX}" -f -p vcf "${TRUTH_P}"
  "${BCFTOOLS}" +setGT "${TRUTH_P}" -- -t a -n u | "${BCFTOOLS}" view -Oz -o "${TRUTH_U}"
  "${TABIX}" -f -p vcf "${TRUTH_U}"
  leaked="$("${BCFTOOLS}" query -f '[%GT]\n' "${TRUTH_U}" | grep -c '|' || true)"
  [[ "${leaked}" -eq 0 ]] || die "phase left in the phaser input: ${leaked}"
  log "${CMRG_REF}: $("${BCFTOOLS}" index -n "${TRUTH_U}") het SNPs in the phaser input ${TRUTH_U}"
  ;;

phase)
  require_file "${TRUTH_U}" "phaser input (run the truth step)"
  eval "$(sed -n '/^phase_one() {$/,/^}$/p' "${PROJECT_ROOT}/phase_run.sh")"
  # names say the reference of this experiment; the BAMs of CHM13v1 are ours
  REFERENCE_NAME="${CMRG_REF}"
  REF="${CREF}"
  if [[ "${CMRG_REF}" == CHM13v1 ]]; then
    sub_bam() { printf '%s/%s.bam\n' "${B}" "$1"; }
  else
    sub_bam() { printf '%s/%s/input.bam\n' "${DATA_DIR}" "$(bam_key "$1" GRCh38 "$3")"; }
  fi
  PHASE_DIR="${D}/phased" LOGS_DIR="${D}/logs" TRUTH_UNPHASED="${TRUTH_U}"
  export -f phase_one vcf_name sub_bam bam_key longphase_flag margin_params \
            hapcut2_flag phaser_applies is_done log warn die
  export PHASE_DIR VCF_DIR LOGS_DIR DATA_DIR REF SAMPLE REFERENCE_NAME REGION_NAME B \
         TRUTH_UNPHASED WHATSHAP LONGPHASE MARGIN MARGIN_PARAMS_DIR HAPCUT2_BIN HIPHASE \
         BCFTOOLS BGZIP TABIX PHASE_THREADS TMP_DIR WHATSHAP_DOWNSAMPLING MARGIN_DEPTH
  log "${CMRG_REF}: truth input x ${PLATFORMS[*]} x ${PHASERS[*]}, ${CMRG_JOBS:-4} at once"
  "${PARALLEL}" -j "${CMRG_JOBS:-4}" phase_one truth {1} "${DEPTH}" {2} ::: "${PLATFORMS[@]}" ::: "${PHASERS[@]}"
  ;;

eval)
  genes="${D}/eval/genes.bed"
  if [[ "${CMRG_REF}" == GRCh38 ]]; then
    awk -F'\t' -v OFS='\t' '{ print $1, $2, $3, $4 }' "${GENES_GRCH38}" >"${genes}"
  else
    # a gene on v1.0: the rows of the v1.0 CMRG BED whose GRCh38 equivalent
    # (GIAB) overlaps the gene on GRCh38, from the first to the last of them, on
    # the chromosome with most of their bases
    require_file "${GENES_EQUIV}" "GRCh38 equivalent of the v1.0 CMRG BED (GIAB CHM13v1.0 SupplementaryFiles)"
    log "placing the CMRG genes on v1.0 by the GRCh38 equivalent of the v1.0 CMRG BED"
    awk -F'\t' -v OFS='\t' '
      FILENAME == ARGV[1] { gc[FNR] = $1; gs[FNR] = $2; ge[FNR] = $3; gn[FNR] = $4; ng = FNR; next }
      FILENAME == ARGV[2] { vc[FNR] = $1; vs[FNR] = $2; ve[FNR] = $3; next }
      { n = substr($4, 10) + 0
        if (!(n in vc)) { print "no row " n " in the v1.0 BED: " $0 > "/dev/stderr"; next }
        for (i = 1; i <= ng; i++) if (gc[i] == $1 && $2 < ge[i] && $3 > gs[i]) {
          k = gn[i] SUBSEP vc[n]
          if (!(k in lo) || vs[n] < lo[k]) lo[k] = vs[n]
          if (ve[n] > hi[k]) hi[k] = ve[n]
          bp[k] += ve[n] - vs[n]
          if ($2 < gs[i] || $3 > ge[i]) over[gn[i]]++
        } }
      END {
        for (k in bp) { split(k, x, SUBSEP); chroms[x[1]]++
          if (bp[k] > best[x[1]]) { best[x[1]] = bp[k]; row[x[1]] = x[2] OFS lo[k] OFS hi[k] OFS x[1] } }
        for (g in row) print row[g]
        for (g in chroms) if (chroms[g] > 1) print "  " g ": rows on " chroms[g] " chromosomes of v1.0, the largest kept" > "/dev/stderr"
        for (g in over) print "  " g ": " over[g] " row(s) reaching past the gene on GRCh38, kept whole" > "/dev/stderr"
      }' "${GENES_GRCH38}" "${CBED}" "${GENES_EQUIV}" | sort -k1,1 -k2,2n >"${genes}"
    log "$(wc -l <"${genes}") of $(wc -l <"${GENES_GRCH38}") genes placed on v1.0"
  fi
  runs=()
  for f in "${D}"/phased/*_truth_*.vcf.gz; do
    b="$(basename "${f}" .vcf.gz)"
    runs+=(--run "${b%%_*}_${b##*_}=${f}")
  done
  [[ ${#runs[@]} -gt 0 ]] || die "no phased VCFs in ${D}/phased"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/cmrg_eval.py" --truth "${CVCF}" --bed "${CBED}" --genes "${genes}" \
      --conf "${GENES_TSV}" "${runs[@]}" --out "${D}/eval/cmrg_${CMRG_REF}"
  log "-> ${D}/eval/cmrg_${CMRG_REF}.tsv, .genes.tsv"
  ;;

tags)
  # the reads of the CMRG regions (primary alignments, chr1-22): where they
  # align, the benchmark SNPs each covers with REF or ALT, and their tags by
  # the benchmark and by every run - for the read recall of the reads step
  require_file "${D}/eval/genes.bed" "gene BED (run the eval step)"
  R="${D}/reads"
  mkdir -p "${R}"
  grep -E "${AUTO_RE}" "${CBED}" >"${R}/cmrg.auto.bed"
  for p in "${PLATFORMS[@]}"; do
    if [[ "${CMRG_REF}" == CHM13v1 ]]; then bam="${B}/${p}.bam"; else bam="$(sub_bam "${p}" GRCh38 "${DEPTH}")"; fi
    require_file "${bam}" "30x BAM of ${p}"
    if [[ ! -s "${R}/${p}.bam.bai" ]]; then
      "${SAMTOOLS}" view -b -F 0x904 -M -L "${R}/cmrg.auto.bed" -o "${R}/${p}.bam" "${bam}"
      "${SAMTOOLS}" index "${R}/${p}.bam"
    fi
    # name, chrom, start, end of each alignment (end from the CIGAR)
    "${SAMTOOLS}" view "${R}/${p}.bam" | awk -F'\t' -v OFS='\t' '{ n = 0; c = $6
        while (match(c, /[0-9]+[MDN=X]/)) { n += substr(c, RSTART, RLENGTH - 1); c = substr(c, RSTART + RLENGTH) }
        print $1, $3, $4 - 1, $4 - 1 + n }' >"${R}/${p}.aln.tsv"
    "${SAMTOOLS}" view "${R}/${p}.bam" | "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/read_snp_count.py" \
        --vcf "${TRUTH_P}" >"${R}/${p}.snps.tsv"
    log "${CMRG_REF} ${p}: $(wc -l <"${R}/${p}.aln.tsv") reads in the CMRG regions, haplotagging"
    "${WHATSHAP}" haplotag --reference "${CREF}" --sample "${SAMPLE}" --ignore-read-groups --skip-missing-contigs \
        --output-haplotag-list "${R}/${p}.truth.tsv" -o /dev/null "${TRUTH_P}" "${R}/${p}.bam" >"${R}/${p}.truth.log" 2>&1
    for vcf in "${D}"/phased/"${p}"_*_truth_*.vcf.gz; do
      name="$(basename "${vcf}" .vcf.gz)"
      "${WHATSHAP}" haplotag --reference "${CREF}" --sample "${SAMPLE}" --ignore-read-groups --skip-missing-contigs \
          --output-haplotag-list "${R}/${name}.tags.tsv" -o /dev/null "${vcf}" "${R}/${p}.bam" >"${R}/${name}.log" 2>&1 \
        || warn "haplotag failed for ${name}, see ${R}/${name}.log"
    done
  done
  ;;

reads)
  # both references at once, anchored on T2T-CHM13 v1.0 (its reads are those of
  # the T2T-CHM13 30x cell, realigned): their consensus labels, and a primary
  # alignment in the GRCh38 30x cell too (its read table)
  lab=()
  for p in "${PLATFORMS[@]}"; do
    f="${READS_DIR}/$(bam_key "${p}" T2T "${DEPTH}").consensus.tsv.gz"
    q="${READS_DIR}/$(bam_key "${p}" GRCh38 "${DEPTH}").reads.tsv.gz"
    if [[ -s "${f}" && -s "${q}" ]]; then lab+=(--labels "${p}=${f}" --require "${p}=${q}"); else warn "${p}: no ${f} or ${q}"; fi
  done
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/cmrg_reads.py" "${lab[@]}" \
      --anchor "CHM13v1=${OUT}/CHM13v1/${CMRG_INPUT}" --anchor-bed "${CMRG_SRC}/HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.bed" \
      --other "GRCh38=${OUT}/GRCh38/${CMRG_INPUT}" --out "${OUT}/cmrg_reads_${CMRG_INPUT}"
  ;;

*) die "unknown step ${STEP}: ref, align, truth, phase, eval, tags or reads" ;;
esac
