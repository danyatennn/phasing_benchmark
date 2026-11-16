#!/usr/bin/env bash
# Pairwise phasing precision and recall (plots/pairwise.py) of every phased VCF of a sample and
# reference in phased/: the grid, experiments A (truthindel, truthall), B
# (MethPhaser, LongPhase 5mC) and C (assemblies), at every depth. Against the
# truth the grid is scored with (TRUTH_PHASED). Read-only for phased/ and the
# truth; writes metrics_pairwise/ only. A row newer than its VCF and the truth
# is kept.
#
#   bash bin/pairwise_run.sh                              HG002 GRCh38
#   REFERENCE_NAME=T2T bash bin/pairwise_run.sh
#   SAMPLE=HG005r2 bash bin/pairwise_run.sh
#
# -> metrics_pairwise/<run>.pairwise.tsv, metrics_pairwise/summary_pairwise_<ref>.tsv

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

require_file "${TRUTH_PHASED}" "phased truth of ${SAMPLE} ${REFERENCE_NAME}"
PAIRWISE_DIR="${SAMPLE_ROOT}/metrics_pairwise"
mkdir -p "${PAIRWISE_DIR}"

pair_one() {
  set -euo pipefail
  local vcf="$1" name out
  name="$(basename "${vcf}" .vcf.gz)"
  out="${PAIRWISE_DIR}/${name}.pairwise.tsv"
  newer_than "${out}" "${vcf}" && newer_than "${out}" "${TRUTH_PHASED}" && return 0
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/pairwise.py" --truth "${TRUTH_PHASED}" --vcf "${vcf}" \
      --name "${name}" >"${out}.$$" && mv -f "${out}.$$" "${out}"
}
export -f pair_one newer_than
export PAIRWISE_DIR TRUTH_PHASED PLOTS_PYTHON PROJECT_ROOT

vcfs=("${PHASE_DIR}"/*_"${REFERENCE_NAME}"_"${SAMPLE}"_*.vcf.gz)
[[ -e "${vcfs[0]}" ]] || die "no phased VCFs of ${SAMPLE} ${REFERENCE_NAME} in ${PHASE_DIR}"
log "${SAMPLE} ${REFERENCE_NAME}: pairwise phasing of ${#vcfs[@]} VCFs, ${PAIRWISE_JOBS:-16} at once"
printf '%s\n' "${vcfs[@]}" | "${PARALLEL}" -j "${PAIRWISE_JOBS:-16}" pair_one {} || true

summary="${PAIRWISE_DIR}/summary_pairwise_${REFERENCE_NAME}.tsv"
files=("${PAIRWISE_DIR}"/*_"${REFERENCE_NAME}"_"${SAMPLE}"_*.pairwise.tsv)
{ head -1 "${files[0]}"; for f in "${files[@]}"; do tail -n +2 "${f}"; done; } >"${summary}.$$"
mv -f "${summary}.$$" "${summary}"
log "${#files[@]} of ${#vcfs[@]} runs -> ${summary}"
