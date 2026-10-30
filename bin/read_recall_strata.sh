#!/usr/bin/env bash
# Read-level haplotag recall per genomic region (plots/read_recall_strata.py),
# for every run of every cell: reads on their truth haplotype / reads the truth
# tags with READ_RECALL_MIN_SNPS truth het SNPs or more, in reads and in read
# bases, genome-wide and per tier-1 stratum.
#
#   metrics_reads/<run>.recall_strata<sfx>.tsv
#   metrics_reads/summary_recall_strata<sfx>.tsv     all runs, one header
#
# READ_TRUTH picks the read labels (sfx: empty for truthtags, _<truth> else):
#   truthtags  whatshap haplotag against the truth VCF of the run's reference;
#              regions, read table and allele matrix of that reference
#   matpat     parental copy of the Q100 assembly (bin/matpat_truth.sh)
#   yak        parental k-mers (bin/yak_truth.sh)
# matpat and yak depend on no reference: the reads are the same on both, so
# every run, GRCh38 or T2T, is scored on the reads, labels, chr1-22, het SNP
# count and strata of REGION_REF (T2T by default) - the same denominator and
# the same regions on both references. Blocks are oriented by their SNPs
# against the truth VCF of the run's reference, as readacc_run.sh; the parental
# labels are put on its haplotypes per chromosome through the cell's truth
# haplotag list, as readacc_truths.sh.
#
# Needs per cell: the read table (bin/read_tables.sh), the truth allele matrix
# (PHASING_CALLERS=truth bash bin/allele_matrix.sh), the haplotag lists and
# phased VCFs of its runs and the truth haplotag list (eval_run.sh), the strata
# cut to chr1-22 (bin/prepare_strata.sh). A table newer than all its inputs is
# kept.
#
#   PHASING_DEPTHS=30 bash bin/read_recall_strata.sh
#   PHASING_DEPTHS=30 REFERENCE_NAME=T2T bash bin/read_recall_strata.sh
#   PHASING_DEPTHS=30 READ_TRUTH=matpat REFERENCE_NAME=GRCh38 bash bin/read_recall_strata.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

RECALL_JOBS="${RECALL_JOBS:-4}"
READ_RECALL_MIN_SNPS="${READ_RECALL_MIN_SNPS:-2}"
READ_TRUTH="${READ_TRUTH:-truthtags}"
case "${READ_TRUTH}" in
  truthtags)  REGION_REF="${REFERENCE_NAME}"; RECALL_SFX="" ;;
  matpat|yak) REGION_REF="${REGION_REF:-T2T}"; RECALL_SFX="_${READ_TRUTH}" ;;
  *) die "READ_TRUTH must be truthtags, matpat or yak, not ${READ_TRUTH}" ;;
esac
mkdir -p "${METRICS_READS_DIR}" "${TMP_DIR}" "${LOGS_DIR}/reads"

STRATA_ARGS=""
for s in $(STRATA_TIER=1 strata_names); do
  STRATA_ARGS="${STRATA_ARGS} --stratum ${s}=$(REFERENCE_NAME="${REGION_REF}" stratum_bed "${s}")"
done
[[ -n "${STRATA_ARGS}" ]] || warn "no strata in ${STRATA_LIST}: genome-wide rows only"

recall_one() {
  set -euo pipefail
  local caller="$1" platform="$2" depth="$3"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  # reads, alleles and labels of the cell on REGION_REF (the cell itself for truthtags)
  local rcell; rcell="$(bam_key "${platform}" "${REGION_REF}" "${depth}")"
  local table="${READS_DIR}/${rcell}.reads.tsv.gz"
  local alleles="${READS_DIR}/${rcell}_truth.reads.tsv.gz"
  local truth_tags="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truthtags).tsv"
  # parental labels: only reads with a primary alignment in the cell on both
  # references; the truth haplotag list of this cell maps PAT/MAT on the truth
  local require="" t ref_tags="${truth_tags}"
  if [[ "${READ_TRUTH}" != truthtags ]]; then
    [[ -s "${ref_tags}" ]] || { warn "${cell}: missing ${ref_tags}"; return 0; }
    require="--truth-haplotags ${ref_tags}"
    truth_tags="${READS_DIR}/${rcell}.${READ_TRUTH}.tsv.gz"
    for t in GRCh38 T2T; do
      t="${READS_DIR}/$(bam_key "${platform}" "${t}" "${depth}").reads.tsv.gz"
      [[ -s "${t}" ]] || { warn "${cell}: missing ${t}"; return 0; }
      require="${require} --require-reads ${t}"
    done
  fi
  for f in "${table}" "${alleles}" "${truth_tags}"; do
    [[ -s "${f}" ]] || { warn "${cell}: missing ${f}"; return 0; }
  done

  local jobs; jobs="$(mktemp -p "${TMP_DIR}")"
  local phaser name tags vcf out
  for phaser in ${PHASER_LIST}; do
    phaser_applies "${phaser}" "${platform}" || continue
    name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
    tags="${METRICS_DIR}/${name}.haplotag.tsv"
    vcf="${PHASE_DIR}/${name}.vcf.gz"
    out="${METRICS_READS_DIR}/${name}.recall_strata${RECALL_SFX}.tsv"
    [[ -s "${tags}" && -s "${vcf}" && ! -d "${METRICS_DIR}/${name}.row.tsv.lock" ]] || continue
    if newer_than "${out}" "${tags}" && newer_than "${out}" "${vcf}" && newer_than "${out}" "${table}" \
        && newer_than "${out}" "${truth_tags}" && newer_than "${out}" "${alleles}" \
        && newer_than "${out}" "${ref_tags}" && ! grep -q $'\treads$' "${out}"; then
      continue
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${REFERENCE_NAME}" "${SAMPLE}" \
        "${caller}" "${platform}" "${depth}" "${phaser}" "${tags}" "${vcf}" "${out}" >>"${jobs}"
  done
  if [[ ! -s "${jobs}" ]]; then
    rm -f "${jobs}"
    log "skip ${cell} ${caller}"
    return 0
  fi
  log "${cell} ${caller}: read recall per region for $(wc -l <"${jobs}") runs"
  # shellcheck disable=SC2086
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/read_recall_strata.py" \
      --reads-table "${table}" --truth-alleles "${alleles}" --truth-tags "${truth_tags}" \
      --truth-source "${READ_TRUTH}" --truth-vcf "${TRUTH_PHASED}" \
      ${STRATA_ARGS} ${require} --min-snps "${READACC_MIN_SNPS}" \
      --min-read-snps "${READ_RECALL_MIN_SNPS}" --jobs "${jobs}" \
    2>>"${LOGS_DIR}/reads/${cell}_${caller}.recall_strata${RECALL_SFX}.log" \
    || { rm -f "${jobs}"; die "read recall failed for ${cell} ${caller}, see ${LOGS_DIR}/reads/${cell}_${caller}.recall_strata${RECALL_SFX}.log"; }
  rm -f "${jobs}"
}

PHASER_LIST="${PHASERS[*]}"
export -f recall_one bam_key vcf_name newer_than phaser_applies stratum_bed log warn die
export METRICS_DIR METRICS_READS_DIR PHASE_DIR READS_DIR LOGS_DIR TMP_DIR REFERENCE_NAME \
       REGION_NAME SAMPLE TRUTH_PHASED STRATA_ARGS READACC_MIN_SNPS READ_RECALL_MIN_SNPS READ_TRUTH REGION_REF RECALL_SFX PREPARED_DIR PLOTS_PYTHON PROJECT_ROOT PHASER_LIST

log "${SAMPLE} ${REFERENCE_NAME}: read recall per region, ${READ_TRUTH} labels, ${REGION_REF} reads and regions, for ${CALLERS[*]} x ${PLATFORMS[*]} x ${DEPTHS[*]}"
"${PARALLEL}" -j "${RECALL_JOBS}" recall_one {1} {2} {3} \
    ::: "${CALLERS[@]}" ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" || true

summary="${METRICS_READS_DIR}/summary_recall_strata${RECALL_SFX}.tsv"
{
  files=("${METRICS_READS_DIR}"/*.recall_strata"${RECALL_SFX}".tsv)
  [[ -e "${files[0]}" ]] && head -1 "${files[0]}"
  for f in "${files[@]}"; do
    [[ -e "${f}" ]] && tail -n +2 "${f}"
  done
  true
} >"${summary}.$$"
mv -f "${summary}.$$" "${summary}"
log "read recall per region: ${summary}"
