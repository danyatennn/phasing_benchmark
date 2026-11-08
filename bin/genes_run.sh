#!/usr/bin/env bash
# Phasing metrics inside a set of genes (CMRG by default) on GRCh38 and T2T,
# from the runs already on disk: nothing is called or phased again, and the
# grid's scripts are run as they are.
#
# Genes: the names of the CMRG gene BED -> NCBI GeneID in the RefSeq annotation
# of GRCh38 -> the gene with that GeneID in the RefSeq annotation of each
# reference (NCBI release RS_2025_08 of GRCh38.p14 and of T2T-CHM13v2.0),
# chr1-22, only genes found on both (plots/gene_regions.py).
#
# Two region modes, GENE_MODE:
#   conf  each gene cut to the confident BED of its reference
#   body  the whole gene. The callers call inside the confident BED only, so a
#         VCF has no variant outside it yet: for now body differs from conf in
#         which reads count (read recall) and in genes fully phased only
#
#   bash bin/genes_run.sh build                    gene BEDs, genes.tsv (% of each gene in the confident BED)
#   bash bin/genes_run.sh prepare                  the truth cut to each gene, both modes and references
#   GENE_MODE=conf REFERENCE_NAME=GRCh38 bash bin/genes_run.sh run
#   GENE_MODE=body REFERENCE_NAME=T2T    bash bin/genes_run.sh run
#
# run, per mode and reference (GENES_STEPS, default conf: regions recall, body: recall fully):
#   regions  regions.sh: callers and phasing per gene and over the whole set
#            (switch, Hamming, phased fraction, ...)
#   recall   bin/read_recall_strata.sh, consensus labels, every labelled read
#            (READ_RECALL_MIN_SNPS=0) that overlaps the region by one base or
#            more (READ_RECALL_FRACTION=0, not half of the read as for the
#            genome strata), the same reads on both references and then each
#            reference's own (GENES_RECALL_COMMON="1 0")
#   fully    plots/genes_fully_phased.py on the whole genes (body only)
#
# Everything is written under
#   regions/genes/<set>/             genes.tsv, <ref>_genes.bed, strata.tsv, <mode>/prepared, <mode>/truth
#   metrics_genes/<set>/<mode>/      strata/summary_{phase,calls}_strata.tsv,
#                                    reads/summary_recall_strata_consensus[_own]_min0_tier*.tsv,
#                                    genes_fully_phased_<set>_<ref>.tsv, logs/
# through a source() hook as in bin/supported_grid.sh: right after a script
# sources config.sh, its strata list, strata BEDs and truth, and the places it
# writes metrics and logs to are moved there. The inputs (phased/, vcf/,
# metrics/ haplotag lists, reads/) are read where they are.
#
#   GENE_SET (CMRG), GENE_NAMES (the CMRG gene BED), PHASING_DEPTHS (30),
#   PHASING_CALLERS (truth, the callers, the assemblies), PHASING_PLATFORMS,
#   PHASING_PHASERS, GENES_STRATA_TIER (regions: empty = the set and every
#   gene, 1 = the set only), GENES_RECALL_TIER (recall: empty = the set and
#   every gene, the default; 1 = the set only)
#
# A trial on a few genes first: GENES_ONLY names them (the whole set is added),
# and everything goes to metrics_genes/<set>_test/<mode>/ instead
#   GENES_ONLY="SMN1 KCNE1 CBS" PHASING_PLATFORMS="hifi ont5k" PHASING_CALLERS="truth clair3" \
#     GENE_MODE=conf REFERENCE_NAME=GRCh38 bash bin/genes_run.sh run

set -euo pipefail

PROJECT_ROOT="/autofs/bal36/dten/phasing_benchmark"
DATA_CONFIG="/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
export PHASING_DEPTHS="${PHASING_DEPTHS:-30}"
export PHASING_CALLERS="${PHASING_CALLERS:-truth clair3 deepvariant longcalld hifiasm hifiasmtrio verkkotrio}"

source "${DATA_CONFIG}"
source "${PROJECT_ROOT}/lib.sh"
source "${PROJECT_ROOT}/config.sh"

[[ "${SAMPLE}" == HG002 ]] || die "genes_run.sh is for HG002 only"

GENE_SET="${GENE_SET:-CMRG}"
GENE_MODE="${GENE_MODE:-conf}"
ANNOT_DIR="${PROJECT_ROOT}/data/annotation"
GENE_NAMES="${GENE_NAMES:-${ANNOT_DIR}/cmrg/GRCh38_CMRG_benchmark_gene_coordinates.bed}"
GENES_DIR="${PROJECT_ROOT}/regions/genes/${GENE_SET}"
GENES_ONLY="${GENES_ONLY:-}"
if [[ -z "${GENES_ONLY}" ]]; then
  GENES_OUT="${PROJECT_ROOT}/metrics_genes/${GENE_SET}/${GENE_MODE}"
else
  GENES_OUT="${PROJECT_ROOT}/metrics_genes/${GENE_SET}_test/${GENE_MODE}"
  # strata_names (lib.sh) keeps only these, in regions.sh and read_recall_strata.sh
  export STRATA_ONLY="${GENE_SET} ${GENES_ONLY}"
fi
[[ "${GENE_MODE}" == conf || "${GENE_MODE}" == body ]] || die "GENE_MODE must be conf or body, not ${GENE_MODE}"

# ref_value <reference> <variable> - a value of config.sh for that reference, in
# a fresh shell (config.sh keeps CONFIDENT_BED and TRUTH_PHASED if already set)
ref_value() {
  env -u CONFIDENT_BED -u TRUTH_PHASED -u TRUTH_SRC_VCF REFERENCE_NAME="$1" bash -c '
    source "$2" >/dev/null 2>&1; source "$3/lib.sh" >/dev/null 2>&1; source "$3/config.sh" >/dev/null 2>&1
    printf "%s\n" "${!1}"' _ "$2" "${DATA_CONFIG}" "${PROJECT_ROOT}"
}

build() {
  local gff_grch38="${ANNOT_DIR}/GRCh38/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
  local rep_grch38="${ANNOT_DIR}/GRCh38/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
  local gff_t2t="${ANNOT_DIR}/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz"
  local rep_t2t="${ANNOT_DIR}/T2T/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt"
  local conf_grch38 conf_t2t f
  conf_grch38="$(ref_value GRCh38 CONFIDENT_BED)"
  conf_t2t="$(ref_value T2T CONFIDENT_BED)"
  for f in "${GENE_NAMES}" "${gff_grch38}" "${rep_grch38}" "${gff_t2t}" "${rep_t2t}" "${conf_grch38}" "${conf_t2t}"; do
    [[ -s "${f}" ]] || die "missing ${f}"
  done
  log "${GENE_SET}: gene regions from $(basename "${GENE_NAMES}")"
  log "  confident BEDs: ${conf_grch38}, ${conf_t2t}"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/gene_regions.py" --set "${GENE_SET}" --names "${GENE_NAMES}" \
      --gff GRCh38="${gff_grch38}" --gff T2T="${gff_t2t}" \
      --assembly-report GRCh38="${rep_grch38}" --assembly-report T2T="${rep_t2t}" \
      --confident GRCh38="${conf_grch38}" --confident T2T="${conf_t2t}" \
      --out "${GENES_DIR}"
  log "${GENES_DIR}/genes.tsv"
}

prepare() {
  [[ -s "${GENES_DIR}/strata.tsv" ]] || die "run build first: no ${GENES_DIR}/strata.tsv"
  local ref truth mode stratum bed out n=0
  for ref in GRCh38 T2T; do
    truth="$(ref_value "${ref}" TRUTH_PHASED)"
    [[ -s "${truth}" ]] || die "missing ${truth}"
    for mode in body conf; do
      mkdir -p "${GENES_DIR}/${mode}/truth"
      for stratum in $(awk -F'\t' 'NR > 1 { print $3 }' "${GENES_DIR}/strata.tsv"); do
        bed="${GENES_DIR}/${mode}/prepared/${ref}_${stratum}.bed"
        out="${GENES_DIR}/${mode}/truth/${ref}_${stratum}.phased.vcf.gz"
        [[ -s "${bed}" ]] || continue
        newer_than "${out}.tbi" "${bed}" && continue
        "${BCFTOOLS}" view -R "${bed}" -Oz -o "${out}" "${truth}"
        "${TABIX}" -f -p vcf "${out}"
        n=$(( n + 1 ))
      done
      log "${ref} ${mode}: truth of $(ls "${GENES_DIR}/${mode}/truth/${ref}_"*.phased.vcf.gz 2>/dev/null | wc -l) regions"
    done
  done
  log "prepare: ${n} truth files written"
}

# after config.sh: the gene set's strata, and every output under GENES_OUT
genes_paths() {
  STRATA_LIST="${GENES_DIR}/strata.tsv"
  PREPARED_DIR="${GENES_DIR}/${GENE_MODE}/prepared"
  STRATA_TRUTH_DIR="${GENES_DIR}/${GENE_MODE}/truth"
  METRICS_STRATA_DIR="${GENES_OUT}/strata"
  METRICS_READS_DIR="${GENES_OUT}/reads"
  LOGS_DIR="${GENES_OUT}/logs"
  mkdir -p "${METRICS_STRATA_DIR}" "${METRICS_READS_DIR}" "${LOGS_DIR}"
}

source() {
  builtin source "$@" || return
  if [[ "$1" == "${PROJECT_ROOT}/config.sh" ]]; then
    genes_paths
  fi
}

run() {
  [[ -d "${GENES_DIR}/${GENE_MODE}/truth" ]] || die "run build and prepare first"
  local g
  for g in ${GENES_ONLY}; do
    awk -F'\t' -v g="${g}" 'NR > 1 && $3 == g { f = 1 } END { exit !f }' "${GENES_DIR}/strata.tsv" \
      || die "GENES_ONLY: ${g} is not in ${GENES_DIR}/strata.tsv (see genes.tsv)"
  done
  local steps="${GENES_STEPS:-}"
  if [[ -z "${steps}" ]]; then
    [[ "${GENE_MODE}" == conf ]] && steps="regions recall" || steps="recall fully"
  fi
  export -f source genes_paths
  mkdir -p "${GENES_OUT}"
  export GENES_DIR GENE_MODE GENES_OUT PROJECT_ROOT
  log "${GENE_SET} ${GENE_MODE} ${REFERENCE_NAME}: ${steps}; ${PHASING_CALLERS} x ${PLATFORMS[*]} x ${PHASING_DEPTHS}x -> ${GENES_OUT}"
  local step common
  for step in ${steps}; do
    case "${step}" in
      regions)
        STRATA_TIER="${GENES_STRATA_TIER:-}" bash "${PROJECT_ROOT}/regions.sh"
        ;;
      recall)
        for common in ${GENES_RECALL_COMMON:-1 0}; do
          READ_TRUTH=consensus READ_RECALL_MIN_SNPS=0 READ_RECALL_COMMON="${common}" READ_RECALL_FRACTION=0 \
            READ_RECALL_TIER="${GENES_RECALL_TIER-}" bash "${PROJECT_ROOT}/bin/read_recall_strata.sh"
        done
        ;;
      fully)
        [[ "${GENE_MODE}" == body ]] || { warn "fully: whole genes only, GENE_MODE=body"; continue; }
        local genes_bed="${GENES_DIR}/${REFERENCE_NAME}_genes.bed"
        if [[ -n "${GENES_ONLY}" ]]; then
          genes_bed="${GENES_OUT}/${REFERENCE_NAME}_genes_only.bed"
          awk -F'\t' -v only="${GENES_ONLY}" 'BEGIN { n = split(only, w, / +/); for (i = 1; i <= n; i++) k[w[i]] = 1 }
              $4 in k' "${GENES_DIR}/${REFERENCE_NAME}_genes.bed" >"${genes_bed}"
        fi
        "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/genes_fully_phased.py" \
            --genes "${genes_bed}" --set "${GENE_SET}_${REFERENCE_NAME}" \
            --refs "${REFERENCE_NAME}" --callers "${CALLERS[*]}" --platforms "${PLATFORMS[*]}" \
            --depths "${DEPTHS[*]}" --phasers "${PHASERS[*]}" --out-dir "${GENES_OUT}" -j 8
        ;;
      *) die "unknown step ${step}: regions, recall or fully" ;;
    esac
  done
  log "${GENE_SET} ${GENE_MODE} ${REFERENCE_NAME}: done, ${GENES_OUT}"
}

case "${1:-}" in
  build) build ;;
  prepare) prepare ;;
  run) run ;;
  *) die "usage: bash bin/genes_run.sh build | prepare | run" ;;
esac
