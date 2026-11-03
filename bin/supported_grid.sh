#!/usr/bin/env bash
# The grid again with the truth cut to the het SNPs that reads support: a truth
# site counts only if at least one primary read of the cell (platform x
# reference x depth) has a base there - any MAPQ, any base quality, a deletion
# over the site does not count (samtools depth without -J, secondary,
# supplementary, QC-failed and duplicate reads left out). The support is per
# cell, so at 10x more sites drop out than at 50x.
#
# Nothing of the grid is changed. Its scripts are run as they are; right after
# each of them sources config.sh, every place it writes to is moved under
# SUPPORTED_ROOT (the source() hook below), and the truth files point to the
# supported ones of the cell. Inputs are read where they are: the downsampled
# BAMs (data/), the callers' VCFs (vcf/), the read tables (reads/), the strata
# (regions/strata). A run whose output path would fall outside SUPPORTED_ROOT
# stops before it writes anything.
#
#   supported/derived/<ref>_<platform>_<depth>x/  the cell's truth files:
#       supported_sites.tsv          chrom, pos of every supported truth record
#       source.supported.vcf.gz      the source truth VCF on those positions
#       truth*.{phased,unphased}.vcf.gz   built from it by make_truth.sh and
#                                    make_truth_variants.sh, as for the grid
#   supported/phased/      truth, truthindel, truthall inputs phased again (the
#                          phasers get the supported sites only); the callers'
#                          VCFs copied from phased/ - their input is unchanged
#   supported/metrics/, metrics_strata_auto/, metrics_switchbed/, meth/
#                          as in the grid, computed against the cell's truth
#   supported/logs/<cell>/, tmp/<cell>/, regions/<cell>/
#
# Stages, per cell in order (SUPPORTED_STAGES picks some):
#   truth     supported sites, the cell's truth files
#   phase     phase_run.sh on truth (and truthindel, truthall at 30x)
#   reuse     copies the callers' phased VCFs and their haplotag lists, and the
#             cell's truth haplotag list, from the grid: the reads and the VCFs
#             are the same, so eval_run.sh does not haplotag them again
#   eval      eval_run.sh, every input
#   analyses  autosomes.sh, readacc_run.sh, calls_fdr.sh, prepare_strata.sh,
#             read_recall_strata.sh with READ_TRUTH=matpat (the read recall of
#             the tables), regions.sh; at 30x also variants_eval.sh and
#             switch_bed.sh
#   meth      experiment B at 30x (ONT, Revio): meth_run.sh on the supported
#             truth-input VCFs, then eval_run.sh and autosomes.sh of its runs.
#             Long; not in the default stages
# and once per reference after the cells:
#   ng50      plots/ng50_autosomes.py over supported/
# Not run: steps that do not depend on the truth sites (read tables, tagged
# bases, read classes, the parental read labels) - a site no read covers
# changes nothing a read carries.
#
# Several servers: each cell takes a lock (supported/locks/<cell>.lock) and
# leaves a mark per finished stage (supported/done/<cell>.<stage>); a cell
# locked by another server is skipped, a finished stage is not run again. A
# lock left by a killed run is removed by hand (its owner file says where).
#
#   REFERENCE_NAME=T2T bash bin/supported_grid.sh                 every cell of T2T
#   SUPPORTED_JOBS=2 PHASING_PLATFORMS=hifi bash bin/supported_grid.sh
#   SUPPORTED_STAGES=paths bash bin/supported_grid.sh             prints where one cell writes, runs nothing
#   SUPPORTED_STAGES=meth REFERENCE_NAME=T2T bash bin/supported_grid.sh
#   python plots/make_plots.py --root /autofs/bal36/dten/phasing_benchmark/supported \
#       --out /autofs/bal36/dten/phasing_benchmark/supported/plots/pdf

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

[[ "${SAMPLE}" == HG002 ]] || die "supported_grid.sh is for HG002 only"

SUPPORTED_ROOT="${SUPPORTED_ROOT:-${PROJECT_ROOT}/supported}"
SUPPORTED_JOBS="${SUPPORTED_JOBS:-4}"
SUPPORTED_STAGES="${SUPPORTED_STAGES:-truth phase reuse eval analyses ng50}"
EXPERIMENT_DEPTH=30
METH_PLATFORMS="ont4k ont5k revio"
METH_PHASERS="lpmod mpwhatshap mplongphase mpmargin mphapcut2 mphiphase"

# the grid's own places, read only
GRID_PHASE_DIR="${PHASE_DIR}"
GRID_METRICS_DIR="${METRICS_DIR}"
GRID_PREPARED_DIR="${PREPARED_DIR}"

case "${SUPPORTED_ROOT}" in
  "${PROJECT_ROOT}"/?*) ;;
  *) die "SUPPORTED_ROOT must be a folder inside ${PROJECT_ROOT}: ${SUPPORTED_ROOT}" ;;
esac
{ "${SAMTOOLS}" depth 2>&1 || true; } | grep -q -- '--excl-flags' \
  || die "${SAMTOOLS} depth has no --excl-flags (samtools 1.13 or newer needed)"
mkdir -p "${SUPPORTED_ROOT}"/{locks,done,logs,data}

# ------------------------------------------------------------ the path hook

# supported_paths - after config.sh: every output place under SUPPORTED_ROOT,
# the truth files those of the cell SUPPORTED_CELL (<platform>_<depth>x)
supported_paths() {
  [[ -n "${SUPPORTED_CELL:-}" && -n "${SUPPORTED_ROOT:-}" ]] \
    || { echo "supported_paths: SUPPORTED_CELL or SUPPORTED_ROOT unset" >&2; exit 1; }
  local cell="${REFERENCE_NAME}_${SUPPORTED_CELL}"
  SAMPLE_ROOT="${SUPPORTED_ROOT}"
  LOGS_DIR="${SUPPORTED_ROOT}/logs/${cell}"
  TMP_DIR="${SUPPORTED_ROOT}/tmp/${cell}"
  SCRATCH_DIR="${TMP_DIR}"
  PHASE_DIR="${SUPPORTED_ROOT}/phased"
  METRICS_DIR="${SUPPORTED_ROOT}/metrics"
  METRICS_STRATA_DIR="${SUPPORTED_ROOT}/$(basename "${METRICS_STRATA_DIR}")"
  METRICS_READS_DIR="${SUPPORTED_ROOT}/metrics_reads"
  SWITCHBED_DIR="${SUPPORTED_ROOT}/metrics_switchbed"
  METH_DIR="${SUPPORTED_ROOT}/meth"
  PREPARED_DIR="${SUPPORTED_ROOT}/regions/${cell}/$(basename "${PREPARED_DIR}")"
  STRATA_TRUTH_DIR="${SUPPORTED_ROOT}/regions/${cell}/$(basename "${STRATA_TRUTH_DIR}")"
  DERIVED_DIR="${SUPPORTED_ROOT}/derived/${cell}"
  TRUTH_SRC_VCF="${DERIVED_DIR}/source.supported.vcf.gz"
  TRUTH_PHASED="${DERIVED_DIR}/truth${_ref_sfx}.het_snps.phased.vcf.gz"
  TRUTH_UNPHASED="${DERIVED_DIR}/truth${_ref_sfx}.het_snps.unphased.vcf.gz"
  TRUTH_INDEL_PHASED="${DERIVED_DIR}/truth${_ref_sfx}.indel.phased.vcf.gz"
  TRUTH_INDEL_UNPHASED="${DERIVED_DIR}/truth${_ref_sfx}.indel.unphased.vcf.gz"
  TRUTH_ALL_PHASED="${DERIVED_DIR}/truth${_ref_sfx}.all.phased.vcf.gz"
  TRUTH_ALL_UNPHASED="${DERIVED_DIR}/truth${_ref_sfx}.all.unphased.vcf.gz"
  local v
  for v in SAMPLE_ROOT LOGS_DIR TMP_DIR SCRATCH_DIR PHASE_DIR METRICS_DIR METRICS_STRATA_DIR \
           METRICS_READS_DIR SWITCHBED_DIR METH_DIR PREPARED_DIR STRATA_TRUTH_DIR DERIVED_DIR \
           TRUTH_SRC_VCF TRUTH_PHASED TRUTH_UNPHASED TRUTH_INDEL_PHASED TRUTH_INDEL_UNPHASED \
           TRUTH_ALL_PHASED TRUTH_ALL_UNPHASED; do
    [[ "${!v}" == "${SUPPORTED_ROOT}" || "${!v}" == "${SUPPORTED_ROOT}"/* ]] \
      || { echo "supported_paths: ${v}=${!v} is outside ${SUPPORTED_ROOT}" >&2; exit 1; }
  done
  mkdir -p "${LOGS_DIR}" "${TMP_DIR}"
}

# every script started from here sources config.sh through this
source() {
  builtin source "$@" || return
  if [[ "$1" == "${PROJECT_ROOT}/config.sh" ]]; then
    supported_paths
  fi
}

export -f source supported_paths
export SUPPORTED_ROOT PROJECT_ROOT

# ------------------------------------------------------------------ stages

# a script of the grid on one cell; PHASING_CALLERS / PHASING_PHASERS as given
grid() {
  PHASING_PLATFORMS="${PLATFORM}" PHASING_DEPTHS="${DEPTH}" bash "$@"
}

stage_paths() {
  bash -c 'source "${PROJECT_ROOT}/config.sh"
           for v in LOGS_DIR TMP_DIR PHASE_DIR METRICS_DIR METRICS_STRATA_DIR SWITCHBED_DIR METH_DIR \
                    PREPARED_DIR STRATA_TRUTH_DIR DERIVED_DIR TRUTH_SRC_VCF TRUTH_PHASED TRUTH_UNPHASED \
                    DATA_DIR VCF_DIR READS_DIR STRATA_DIR; do printf "%-20s %s\n" "${v}" "${!v}"; done'
}

stage_truth() {
  local derived="${SUPPORTED_ROOT}/derived/${CELL}"
  local src="${derived}/source.supported.vcf.gz" sites="${derived}/supported_sites.tsv"
  local bam; bam="$(sub_bam "${PLATFORM}" "${REFERENCE_NAME}" "${DEPTH}")"
  require_file "${bam}" "BAM of ${CELL}"
  mkdir -p "${derived}"
  if ! is_done "${src}"; then
    local tmp; tmp="$(mktemp -d -p "${SUPPORTED_ROOT}/tmp")"
    # every record of the grid's source truth, as one-base intervals
    "${BCFTOOLS}" query -f '%CHROM\t%POS\n' "${TRUTH_SRC_VCF}" \
      | awk -v OFS='\t' '{ print $1, $2 - 1, $2 }' >"${tmp}/sites.bed"
    # positions with a base of a primary read: without -a only those are printed
    "${SAMTOOLS}" depth -b "${tmp}/sites.bed" --excl-flags 0xF04 "${bam}" \
      | cut -f1,2 >"${tmp}/sites.tsv"
    "${BCFTOOLS}" view -T "${tmp}/sites.tsv" -Oz -o "${tmp}/src.vcf.gz" "${TRUTH_SRC_VCF}"
    "${TABIX}" -f -p vcf "${tmp}/src.vcf.gz"
    mv -f "${tmp}/sites.tsv" "${sites}"
    mv -f "${tmp}/src.vcf.gz.tbi" "${src}.tbi"
    mv -f "${tmp}/src.vcf.gz" "${src}"
    rm -rf "${tmp}"
    log "${CELL}: $("${BCFTOOLS}" index -n "${src}") of $("${BCFTOOLS}" index -n "${TRUTH_SRC_VCF}") truth records supported"
  fi
  # TRUTH_SRC_VCF of the scripts is the supported source of the cell (hook)
  grid bin/make_truth.sh "${REFERENCE_NAME}"
  [[ "${DEPTH}" != "${EXPERIMENT_DEPTH}" ]] || grid bin/make_truth_variants.sh
}

stage_phase() {
  local callers="truth"
  [[ "${DEPTH}" != "${EXPERIMENT_DEPTH}" ]] || callers="truth truthindel truthall"
  PHASING_CALLERS="${callers}" grid phase_run.sh
}

stage_reuse() {
  local caller phaser name dst="${SUPPORTED_ROOT}/phased"
  mkdir -p "${dst}" "${SUPPORTED_ROOT}/metrics"
  for caller in "${VARIANT_CALLERS[@]}"; do
    for phaser in "${PHASERS[@]}"; do
      phaser_applies "${phaser}" "${PLATFORM}" || continue
      name="$(vcf_name "${PLATFORM}" "${DEPTH}" "${caller}" "${phaser}")"
      if ! is_done "${GRID_PHASE_DIR}/${name}.vcf.gz"; then
        warn "${name}: no phased VCF in the grid, not reused"
        continue
      fi
      # the VCF keeps its old time, the haplotag list is copied after it, so
      # eval_run.sh sees the list as made from this VCF
      is_done "${dst}/${name}.vcf.gz" \
        || { cp -p "${GRID_PHASE_DIR}/${name}.vcf.gz" "${dst}/"; cp -p "${GRID_PHASE_DIR}/${name}.vcf.gz.tbi" "${dst}/"; }
      [[ -s "${GRID_METRICS_DIR}/${name}.haplotag.tsv" && ! -s "${SUPPORTED_ROOT}/metrics/${name}.haplotag.tsv" ]] \
        && cp "${GRID_METRICS_DIR}/${name}.haplotag.tsv" "${SUPPORTED_ROOT}/metrics/"
    done
  done
  # the truth haplotag list: a site no read covers tags no read, so the
  # supported truth tags the reads as the full one
  name="$(vcf_name "${PLATFORM}" "${DEPTH}" truthtags).tsv"
  [[ -s "${GRID_METRICS_DIR}/${name}" && ! -s "${SUPPORTED_ROOT}/metrics/${name}" ]] \
    && cp "${GRID_METRICS_DIR}/${name}" "${SUPPORTED_ROOT}/metrics/"
  return 0
}

stage_eval() {
  local callers="truth ${VARIANT_CALLERS[*]}"
  [[ "${DEPTH}" != "${EXPERIMENT_DEPTH}" ]] || callers="${callers} truthindel truthall"
  PHASING_CALLERS="${callers}" grid eval_run.sh
}

# read haplotagging recall against the parental labels (matpat), as in the
# grid's tables: the reads, their labels, allele counts and strata are those of
# T2T (READ_RECALL_COMMON), read where the grid keeps them; only the blocks are
# oriented by the cell's supported truth. A read covers only sites a read
# supports, so the reads counted are the grid's. The T2T strata BEDs are copied
# into the cell's prepared folder, where read_recall_strata.sh looks for them.
recall_matpat() {
  local prepared b
  prepared="${SUPPORTED_ROOT}/regions/${CELL}/$(basename "${GRID_PREPARED_DIR}")"
  mkdir -p "${prepared}"
  for b in "${GRID_PREPARED_DIR}"/T2T_*.bed; do
    if [[ -e "${b}" && ! -s "${prepared}/$(basename "${b}")" ]]; then
      cp -p "${b}" "${prepared}/"
    fi
  done
  READ_TRUTH=matpat grid bin/read_recall_strata.sh
}

stage_analyses() {
  local callers="truth ${VARIANT_CALLERS[*]}"
  [[ "${DEPTH}" != "${EXPERIMENT_DEPTH}" ]] || callers="${callers} truthindel truthall"
  PHASING_CALLERS="${callers}" grid bin/autosomes.sh
  PHASING_CALLERS="${callers}" grid readacc_run.sh
  grid bin/calls_fdr.sh
  grid bin/prepare_strata.sh
  PHASING_CALLERS="${callers}" recall_matpat
  PHASING_CALLERS="truth ${VARIANT_CALLERS[*]}" grid regions.sh
  if [[ "${DEPTH}" == "${EXPERIMENT_DEPTH}" ]]; then
    grid bin/variants_eval.sh
    grid bin/switch_bed.sh
  fi
}

stage_meth() {
  [[ "${DEPTH}" == "${EXPERIMENT_DEPTH}" && " ${METH_PLATFORMS} " == *" ${PLATFORM} "* ]] || return 0
  grid meth_run.sh
  PHASING_CALLERS=truth PHASING_PHASERS="${METH_PHASERS}" grid eval_run.sh
  PHASING_CALLERS=truth PHASING_PHASERS="${METH_PHASERS}" grid bin/autosomes.sh
  PHASING_CALLERS=truth PHASING_PHASERS="${METH_PHASERS}" recall_matpat
}

# cell_run <platform> <depth> - the stages of one cell, under its lock
cell_run() {
  set -euo pipefail
  PLATFORM="$1" DEPTH="$2"
  export SUPPORTED_CELL="${PLATFORM}_${DEPTH}x"
  CELL="${REFERENCE_NAME}_${SUPPORTED_CELL}"
  if [[ " ${SUPPORTED_STAGES} " == *" paths "* ]]; then
    stage_paths
    return 0
  fi
  local lock="${SUPPORTED_ROOT}/locks/${CELL}.lock"
  mkdir "${lock}" 2>/dev/null || { log "busy elsewhere ${CELL} ($(cat "${lock}/owner" 2>/dev/null))"; return 0; }
  printf '%s %s %s\n' "$(hostname -s)" "$$" "$(date '+%F %H:%M')" >"${lock}/owner"
  trap "rm -f '${lock}/owner'; rmdir '${lock}' 2>/dev/null" EXIT INT TERM HUP
  mkdir -p "${SUPPORTED_ROOT}/logs/${CELL}" "${SUPPORTED_ROOT}/tmp"
  local stage mark status logf="${SUPPORTED_ROOT}/logs/${CELL}/supported_grid.log"
  for stage in truth phase reuse eval analyses meth; do
    [[ " ${SUPPORTED_STAGES} " == *" ${stage} "* ]] || continue
    mark="${SUPPORTED_ROOT}/done/${CELL}.${stage}"
    [[ -e "${mark}" ]] && { log "${CELL} ${stage}: done before"; continue; }
    log "${CELL} ${stage}"
    # in a subshell of its own, so set -e stops the stage at its first error
    set +e
    ( set -euo pipefail; "stage_${stage}" ) >>"${logf}" 2>&1
    status=$?
    set -e
    [[ "${status}" -eq 0 ]] || die "${CELL} ${stage} failed, see ${logf}"
    date '+%F %H:%M' >"${mark}"
  done
  log "${CELL}: ${SUPPORTED_STAGES}"
}

export -f grid cell_run stage_paths stage_truth stage_phase stage_reuse stage_eval stage_analyses \
          stage_meth recall_matpat sub_bam bam_key vcf_name phaser_applies is_done require_file log warn die
export REFERENCE_NAME SAMPLE REGION_NAME SUPPORTED_STAGES EXPERIMENT_DEPTH METH_PLATFORMS METH_PHASERS \
       GRID_PHASE_DIR GRID_METRICS_DIR GRID_PREPARED_DIR TRUTH_SRC_VCF DATA_DIR BCFTOOLS TABIX SAMTOOLS
export VARIANT_CALLERS_LIST="${VARIANT_CALLERS[*]}" PHASERS_LIST="${PHASERS[*]}"

# arrays do not cross into the parallel jobs: rebuilt there from the lists
cell_job() {
  VARIANT_CALLERS=(${VARIANT_CALLERS_LIST})
  PHASERS=(${PHASERS_LIST})
  cell_run "$@"
}
export -f cell_job

# ------------------------------------------------------------------- run

if [[ " ${SUPPORTED_STAGES} " == *" paths "* ]]; then
  cell_job "${PLATFORMS[0]}" "${DEPTHS[0]}"
  exit 0
fi

log "${SAMPLE} ${REFERENCE_NAME}, supported truth sites: ${PLATFORMS[*]} x ${DEPTHS[*]}x, stages ${SUPPORTED_STAGES}, ${SUPPORTED_JOBS} cells at a time -> ${SUPPORTED_ROOT}"

"${PARALLEL}" -j "${SUPPORTED_JOBS}" --line-buffer cell_job {1} {2} \
    ::: "${PLATFORMS[@]}" ::: "${DEPTHS[@]}" \
  || warn "some cells failed, see ${SUPPORTED_ROOT}/logs/<cell>/supported_grid.log"

# each cell rewrites the read recall summary with what is there at that moment;
# once more after all of them, on whatever every server has finished
recall_files=("${SUPPORTED_ROOT}"/metrics_reads/*.recall_strata_matpat.tsv)
if [[ -e "${recall_files[0]}" ]]; then
  {
    head -1 "${recall_files[0]}"
    for f in "${recall_files[@]}"; do tail -n +2 "${f}"; done
  } >"${SUPPORTED_ROOT}/metrics_reads/summary_recall_strata_matpat.tsv.$$"
  mv -f "${SUPPORTED_ROOT}/metrics_reads/summary_recall_strata_matpat.tsv.$$" \
        "${SUPPORTED_ROOT}/metrics_reads/summary_recall_strata_matpat.tsv"
fi

if [[ " ${SUPPORTED_STAGES} " == *" ng50 "* ]]; then
  # ng50_autosomes.py reads <root>/data/chrom_lengths_<ref>.tsv
  cp -f "${CHROM_LENGTHS}" "${SUPPORTED_ROOT}/data/"
  log "${REFERENCE_NAME}: NG50 over supported/"
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/ng50_autosomes.py" --root "${SUPPORTED_ROOT}" \
      --refs "${REFERENCE_NAME}" -j 16
  "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/ng50_autosomes.py" --root "${SUPPORTED_ROOT}" \
      --refs "${REFERENCE_NAME}" --callers "truthindel truthall" --depths "${EXPERIMENT_DEPTH}" -j 16
  if [[ " ${SUPPORTED_STAGES} " == *" meth "* ]]; then
    "${PLOTS_PYTHON}" "${PROJECT_ROOT}/plots/ng50_autosomes.py" --root "${SUPPORTED_ROOT}" \
        --refs "${REFERENCE_NAME}" --callers truth --phasers "${METH_PHASERS}" \
        --depths "${EXPERIMENT_DEPTH}" -j 16
  fi
fi

log "supported grid ${REFERENCE_NAME}: ${SUPPORTED_STAGES} done"
