#!/usr/bin/env bash
# Where the reads of each stratum land on the other reference, in both
# directions (T2T-CHM13 -> GRCh38, GRCh38 -> T2T-CHM13), 30x, every platform:
# plots/strata_crossref.py, for the upset plot of plots/figures_upset.py.
# Read-only for every pipeline output; writes one table.
#
#   bash bin/strata_crossref.sh
#
# -> metrics_reads/strata_crossref_sets_30x.tsv (an older one is renamed, not removed)

set -euo pipefail

cd /autofs/bal36/dten/phasing_benchmark
PY=/autofs/bal36/dten/miniforge3/envs/plots/bin/python
OUT=metrics_reads/strata_crossref_sets_30x.tsv
[[ -e "${OUT}" ]] && mv "${OUT}" "${OUT}.old.$(date +%s)"
STRATA="notinalldifficultregions alldifficultregions segdups lowmappabilityall AllTandemRepeatsandHomopolymers_slop5 satellites_slop5 MHC"

# config.sh of one reference, in a shell of its own: confident BED, the
# strata cut to it (prepared/), the folder and name of the GIAB strata
ref_paths() {
  env -u CONFIDENT_BED -u PREPARED_DIR -u STRATA_DIR REFERENCE_NAME="$1" bash -c '
    source /autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh >/dev/null 2>&1
    source lib.sh >/dev/null 2>&1; source config.sh >/dev/null 2>&1
    echo "${CONFIDENT_BED} ${PREPARED_DIR} ${STRATA_DIR} ${STRATA_REF}"'
}
read -r G_CONF G_PREP G_SDIR G_SREF < <(ref_paths GRCh38)
read -r T_CONF T_PREP T_SDIR T_SREF < <(ref_paths T2T)
echo "GRCh38: ${G_CONF} ${G_PREP} ${G_SREF}" >&2
echo "T2T:    ${T_CONF} ${T_PREP} ${T_SREF}" >&2

# one direction: <source> <other> <source prepared dir> <other confident> <other strata dir> <other strata name>
direction() {
  local src=$1 dst=$2 sprep=$3 dconf=$4 dsdir=$5 dsref=$6 p s args sc dc
  for p in ont4k ont5k hifi revio; do
    args=()
    for s in ${STRATA}; do
      args+=(--stratum "${s}=${sprep}/${src}_${s}.bed=$(find "${dsdir}/${dsref}@all" -name "${dsref}_${s}.bed.gz" | head -1)")
    done
    args+=(--stratum "CMRG=regions/genes/CMRG/conf/prepared/${src}_CMRG.bed=regions/genes/CMRG/body/prepared/${dst}_CMRG.bed")
    sc=${p}_${src}_HG002_30x_wgs
    dc=${p}_${dst}_HG002_30x_wgs
    "${PY}" plots/strata_crossref.py --labels "reads/${sc}.consensus.tsv.gz" \
        --src-table "reads/${sc}.reads.tsv.gz" --dst-table "reads/${dc}.reads.tsv.gz" \
        --src-alleles "reads/${sc}_truth.reads.tsv.gz" --dst-conf "${dconf}" "${args[@]}" \
        --direction "${src}_to_${dst}" --platform "${p}" --out "${OUT}"
  done
}

direction T2T GRCh38 "${T_PREP}" "${G_CONF}" "${G_SDIR}" "${G_SREF}"
direction GRCh38 T2T "${G_PREP}" "${T_CONF}" "${T_SDIR}" "${T_SREF}"
echo "-> ${OUT}" >&2
