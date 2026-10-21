#!/usr/bin/env bash
# The tables of the paper (plots/make_tables.py), from what is already on disk:
# nothing is computed again except the read statistics of the 30x cells and
# the variant classes of the truth VCFs. Read-only for every pipeline output.
#
#   bash bin/make_tables.sh
#   TABLE_PHASERS="whatshap longphase hapcut2" bash bin/make_tables.sh
#   TABLES_DIR=/somewhere bash bin/make_tables.sh
#
# Writes <TABLES_DIR>/*.tsv, tables.xlsx and README.tsv (default tables/), plus
# the two files it hands to the python side:
#   inputs.tsv     sample, reference, item, value: truth VCF, BED, phased truth,
#                  source BAM, basecaller model and aligner of every platform
#   software.tsv   tool, version, where it was read from

set -euo pipefail

cd /autofs/bal36/dten/phasing_benchmark

OUT="${TABLES_DIR:-/autofs/bal36/dten/phasing_benchmark/tables}"
TABLE_PHASERS="${TABLE_PHASERS:-whatshap longphase margin hapcut2}"
mkdir -p "${OUT}"

# config.sh keeps CONFIDENT_BED and TRUTH_SRC_VCF if they are already set, so
# every sample and reference reads it in a shell of its own
sample_inputs() {   # <sample> <reference>
  (
    export SAMPLE="$1" REFERENCE_NAME="$2"
    source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
    source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
    source "/autofs/bal36/dten/phasing_benchmark/config.sh"
    set +e
    [[ -s "${TRUTH_SRC_VCF:-}" ]] || exit 0
    printf '%s\t%s\t%s\t%s\n' "${SAMPLE}" "${REFERENCE_NAME}" truth_vcf "${TRUTH_SRC_VCF}"
    printf '%s\t%s\t%s\t%s\n' "${SAMPLE}" "${REFERENCE_NAME}" bed "${CONFIDENT_BED}"
    printf '%s\t%s\t%s\t%s\n' "${SAMPLE}" "${REFERENCE_NAME}" truth_phased "${TRUTH_PHASED}"
    for p in "${PLATFORMS[@]}"; do
      bam="$(source_bam "${p}")"
      [[ -s "${bam}" ]] || continue
      printf '%s\t%s\t%s\t%s\n' "${SAMPLE}" "${REFERENCE_NAME}" "bam_${p}" "${bam}"
      header="$("${SAMTOOLS}" view -H "${bam}" 2>/dev/null)"
      model="$(grep -oE 'basecall_model=[^[:space:]]+|dna_r[0-9][^[:space:],]*' <<<"${header}" | sed 's/^basecall_model=//' | sort -u | head -1)"
      aligner="$(grep '^@PG' <<<"${header}" | grep -oiE 'ID:(minimap2|pbmm2|winnowmap)[^[:space:]]*|VN:[^[:space:]]+' | head -2 | tr '\n' ' ' | sed 's/ $//')"
      printf '%s\t%s\t%s\t%s\n' "${SAMPLE}" "${REFERENCE_NAME}" "model_${p}" "${model}"
      printf '%s\t%s\t%s\t%s\n' "${SAMPLE}" "${REFERENCE_NAME}" "aligner_${p}" "${aligner}"
    done
  )
}

samples=(HG002)
while read -r s; do
  [[ -d "samples/${s}" ]] && samples+=("${s}")
done < <(awk -F'\t' '!/^#/ && NF >= 3 && $1 != "HG002" { print $1 }' samples.tsv | sort -u)

{
  printf 'sample\treference\titem\tvalue\n'
  for s in "${samples[@]}"; do
    for r in GRCh38 T2T; do
      sample_inputs "${s}" "${r}"
    done
  done
} >"${OUT}/inputs.tsv.tmp"
mv -f "${OUT}/inputs.tsv.tmp" "${OUT}/inputs.tsv"

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"
set +e

# versions: conda packages from conda-meta (name-version-build.json), the rest
# from file names or --version; whatever cannot be read is '?'
conda_versions() {   # <env bin dir> <package>...
  local meta="$1/../conda-meta" pkg f
  shift
  for pkg in "$@"; do
    f="$(ls "${meta}" 2>/dev/null | grep -E "^${pkg}-[0-9]" | head -1)"
    [[ -n "${f}" ]] || continue
    f="${f%.json}"; f="${f#${pkg}-}"
    printf '%s\t%s\tconda %s\n' "${pkg}" "${f%%-*}" "$(basename "$(dirname "$1")")"
  done
}
first_version() {   # <command...>: the first x.y[.z] it prints
  timeout 60 "$@" 2>&1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1
}
{
  printf 'tool\tversion\twhere\n'
  conda_versions "${PHASING_BIN}" whatshap bcftools htslib pysam python mosdepth parallel
  conda_versions "${HAPCUT2_BIN}" hapcut2
  conda_versions "${TRIO_BIN}" yak meryl merqury
  conda_versions "${VERKKO_BIN}" verkko
  conda_versions "${DIPCALL_BIN}" dipcall
  printf 'samtools\t%s\t%s\n' "$(first_version "${SAMTOOLS}" --version)" "${SAMTOOLS}"
  printf 'minimap2\t%s\t%s\n' "$(first_version "${MINIMAP2}" --version)" "${MINIMAP2}"
  printf 'LongPhase\t%s\t%s\n' "$(first_version "${LONGPHASE}" --version)" "${LONGPHASE}"
  printf 'longcallD\t%s\t%s\n' "$(grep -oE 'v[0-9]+\.[0-9]+\.[0-9]+' <<<"${LONGCALLD}" | head -1 | tr -d v)" "${LONGCALLD}"
  printf 'DeepVariant\t%s\t%s\n' "$(grep -oE '[0-9]+\.[0-9]+\.[0-9]+' <<<"$(basename "${DEEPVARIANT_SIF}")" | head -1)" "${DEEPVARIANT_SIF}"
  printf 'Margin\t%s\t%s\n' "$(git -C "$(dirname "$(dirname "${MARGIN}")")" describe --tags 2>/dev/null || echo '?')" "${MARGIN}"
  printf 'hifiasm\t%s\t%s\n' "$(first_version "${HIFIASM}" --version)" "${HIFIASM}"
  printf 'Clair3\t?\t%s (models: %s, %s, %s, %s)\n' "${CLAIR3_SIF}" \
    "${CLAIR3_MODEL_ont4k}" "${CLAIR3_MODEL_ont5k}" "${CLAIR3_MODEL_hifi}" "${CLAIR3_MODEL_revio}"
} | awk -F'\t' -v OFS='\t' '{ if ($2 == "") $2 = "?"; print }' >"${OUT}/software.tsv"

"${PLOTS_PYTHON}" plots/make_tables.py \
    --root "${PROJECT_ROOT}" --out "${OUT}" \
    --inputs "${OUT}/inputs.tsv" --software "${OUT}/software.tsv" \
    --phasers "${TABLE_PHASERS}" --threads "${TABLE_THREADS:-8}" \
    ${TABLE_ARGS:-}
