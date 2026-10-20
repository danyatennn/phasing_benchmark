#!/usr/bin/env bash
# Moves the outputs of chosen runs out of the way, so the pipeline computes
# them again. Nothing is deleted: everything goes to set_aside/<label>/ under
# the sample's root, with the same relative paths. A dry run by default.
#
#   bash bin/set_aside.sh <label> run|row
#
#   run   the phased VCF and its index, every metrics/<run>.* file (row,
#         stats, compare, haplotag list, read accuracy, NG50, ...), and the
#         per-stratum rows phase_<run>__* - the run is phased and scored anew.
#         An assembly's dipcall work directory (<run>.dipcall) stays, so
#         dipcall_run.sh rebuilds the VCF from it without running dipcall.
#   row   metrics/<run>.row.tsv only - the run is scored again from the files
#         it already has (eval_run.sh reuses stats and haplotag lists that are
#         newer than the VCF).
#
# Which runs: the grid of config.sh (PHASING_* variables, REFERENCE_NAME,
# SAMPLE), and SET_ASIDE_REFS for several references at once.
#
#   PHASING_PHASERS=margin SET_ASIDE_REFS="GRCh38 T2T" bash bin/set_aside.sh margin_p30 run
#   ... CONFIRM=1 bash bin/set_aside.sh margin_p30 run      to actually move
#
# A run with a lock (being written by some host) stops the whole move.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

label="${1:-}"
what="${2:-}"
[[ "${label}" =~ ^[A-Za-z0-9._-]+$ ]] || die "usage: bash bin/set_aside.sh <label> run|row"
[[ "${what}" == run || "${what}" == row ]] || die "usage: bash bin/set_aside.sh <label> run|row"

REFS=(${SET_ASIDE_REFS:-${REFERENCE_NAME}})
dest_root="${SAMPLE_ROOT}/set_aside/${label}"

names=()
for ref in "${REFS[@]}"; do
  for caller in "${CALLERS[@]}"; do
    for platform in "${PLATFORMS[@]}"; do
      for depth in "${DEPTHS[@]}"; do
        for phaser in "${PHASERS[@]}"; do
          names+=("${platform}_${ref}_${SAMPLE}_${depth}x_${REGION_NAME}_${caller}_${phaser}")
        done
      done
    done
  done
done

# every file of a run, relative to SAMPLE_ROOT
files_of() {
  local name="$1" f
  if [[ "${what}" == row ]]; then
    [[ -e "${METRICS_DIR}/${name}.row.tsv" ]] && printf '%s\n' "metrics/${name}.row.tsv"
    return 0
  fi
  for f in "${PHASE_DIR}/${name}.vcf.gz" "${PHASE_DIR}/${name}.vcf.gz.tbi" "${PHASE_DIR}/${name}.dip.bed"; do
    [[ -e "${f}" ]] && printf '%s\n' "${f#"${SAMPLE_ROOT}"/}"
  done
  for f in "${METRICS_DIR}/${name}".* \
           "${SAMPLE_ROOT}"/metrics_strata*/phase_"${name}"__* \
           "${METRICS_READS_DIR}/${name}".* "${SWITCHBED_DIR}/${name}".*; do
    [[ -f "${f}" ]] && printf '%s\n' "${f#"${SAMPLE_ROOT}"/}"
  done
  return 0
}

locked=()
moves=()
for name in "${names[@]}"; do
  for lock in "${PHASE_DIR}/${name}.vcf.gz.lock" "${METRICS_DIR}/${name}.row.tsv.lock"; do
    [[ -d "${lock}" ]] && locked+=("${lock} ($(cat "${lock}/owner" 2>/dev/null))")
  done
  while IFS= read -r rel; do
    [[ -n "${rel}" ]] && moves+=("${rel}")
  done < <(files_of "${name}")
done

[[ ${#locked[@]} -eq 0 ]] || { printf '  %s\n' "${locked[@]}" >&2; die "${#locked[@]} runs are being written right now - nothing moved"; }

log "${SAMPLE}: ${#names[@]} runs (${REFS[*]} x ${CALLERS[*]} x ${PLATFORMS[*]} x ${DEPTHS[*]} x ${PHASERS[*]}), ${what}: ${#moves[@]} files -> ${dest_root}"
[[ ${#moves[@]} -gt 0 ]] || { log "nothing to move"; exit 0; }

if [[ -z "${CONFIRM:-}" ]]; then
  printf '%s\n' "${moves[@]}" | head -n 20 | sed 's/^/  /'
  [[ ${#moves[@]} -gt 20 ]] && printf '  ... and %d more\n' "$(( ${#moves[@]} - 20 ))"
  log "dry run, nothing moved: CONFIRM=1 to move"
  exit 0
fi

for rel in "${moves[@]}"; do
  [[ -e "${dest_root}/${rel}" ]] && die "${dest_root}/${rel} exists already - pick another label"
done
for rel in "${moves[@]}"; do
  mkdir -p "$(dirname "${dest_root}/${rel}")"
  mv "${SAMPLE_ROOT}/${rel}" "${dest_root}/${rel}"
done
log "moved ${#moves[@]} files to ${dest_root}"
