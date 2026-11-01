#!/usr/bin/env bash
# Wall-clock time and peak memory of every phaser at 30x, every caller input
# (the full grid of one reference), measured with GNU time.
#
# The phasing commands are not copied: phase_one is read from phase_run.sh as
# it stands, so the runs here are the runs of the grid. Everything this script
# writes goes to runtime/ (RUNTIME_DIR) - its own phased VCFs, logs and table;
# phased/, metrics/ and logs/ of the grid are only read, never written.
#
# One run at a time, so the runs do not share the machine with each other. Run
# it on a node with nothing else heavy on it; the host and its load average
# before each run are written next to the numbers.
#
#   runtime/runtime_30x_<ref>.tsv, one line per run:
#     reference platform depth caller phaser threads wall_s user_s sys_s
#     max_rss_mb exit host load1 date
#   wall_s       elapsed time of the whole phase_one (phaser, plus reading the
#                input VCF and bgzip/tabix of the output: seconds)
#   max_rss_mb   peak resident memory of the largest process (GNU time)
#   threads      PHASE_THREADS, given to LongPhase, Margin and HiPhase;
#                WhatsHap and HapCUT2 run on one thread
# A run already in the table with exit 0 is skipped; a failed one is run again.
#
#   bash bin/runtime_run.sh                             GRCh38
#   REFERENCE_NAME=T2T bash bin/runtime_run.sh
#   PHASING_CALLERS=truth PHASING_PHASERS="whatshap longphase" bash bin/runtime_run.sh
#
# Heatmaps: python plots/figures_runtime.py

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

RUNTIME_DEPTH=30
RUNTIME_DIR="${RUNTIME_DIR:-${SAMPLE_ROOT}/runtime}"
GNU_TIME="${GNU_TIME:-/usr/bin/time}"
[[ -x "${GNU_TIME}" ]] && "${GNU_TIME}" -v true >/dev/null 2>&1 \
  || die "no GNU time at ${GNU_TIME} (set GNU_TIME)"

# phase_one writes its VCF to PHASE_DIR and its log to LOGS_DIR/phase: both
# point into runtime/ before the function is read
PHASE_DIR="${RUNTIME_DIR}/phased"
LOGS_DIR="${RUNTIME_DIR}/logs"
mkdir -p "${PHASE_DIR}" "${LOGS_DIR}/phase" "${LOGS_DIR}/time" "${TMP_DIR}"

eval "$(sed -n '/^phase_one() {$/,/^}$/p' "${PROJECT_ROOT}/phase_run.sh")"
declare -F phase_one >/dev/null || die "phase_one not found in ${PROJECT_ROOT}/phase_run.sh"

# as phase_run.sh exports them, for the bash that GNU time starts
export -f phase_one vcf_name sub_bam bam_key longphase_flag margin_params \
          hapcut2_flag phaser_applies is_done log warn die
export PHASE_DIR VCF_DIR LOGS_DIR DATA_DIR REF SAMPLE REFERENCE_NAME REGION_NAME \
       TRUTH_UNPHASED TRUTH_INDEL_UNPHASED TRUTH_ALL_UNPHASED WHATSHAP LONGPHASE MARGIN MARGIN_PARAMS_DIR HAPCUT2_BIN HIPHASE \
       BCFTOOLS BGZIP TABIX PHASE_THREADS TMP_DIR WHATSHAP_DOWNSAMPLING MARGIN_DEPTH

table="${RUNTIME_DIR}/runtime_30x_${REFERENCE_NAME}.tsv"
[[ -s "${table}" ]] || printf 'reference\tplatform\tdepth\tcaller\tphaser\tthreads\twall_s\tuser_s\tsys_s\tmax_rss_mb\texit\thost\tload1\tdate\n' >"${table}"

# time_value <GNU time -v output> <label> - the value after "label: "
time_value() { awk -F': ' -v k="$2" '$1 ~ k { print $NF; exit }' "$1"; }

# h:mm:ss or m:ss.ss -> seconds
to_seconds() { awk -F: '{ s = 0; for (i = 1; i <= NF; i++) s = s * 60 + $i; printf "%.1f", s }' <<<"$1"; }

log "${SAMPLE} ${REFERENCE_NAME} ${RUNTIME_DEPTH}x: ${CALLERS[*]} x ${PLATFORMS[*]} x ${PHASERS[*]}, one run at a time -> ${table}"

for caller in "${CALLERS[@]}"; do
  for platform in "${PLATFORMS[@]}"; do
    for phaser in "${PHASERS[@]}"; do
      phaser_applies "${phaser}" "${platform}" || continue
      if awk -F'\t' -v p="${platform}" -v c="${caller}" -v f="${phaser}" -v d="${RUNTIME_DEPTH}" \
          '$2 == p && $3 == d && $4 == c && $5 == f && $11 == 0 { found = 1 } END { exit !found }' "${table}"; then
        log "skip ${platform} ${caller} ${phaser}: measured"
        continue
      fi
      [[ -s "$(sub_bam "${platform}" "${REFERENCE_NAME}" "${RUNTIME_DEPTH}")" ]] \
        || { warn "no ${RUNTIME_DEPTH}x BAM for ${platform}, skipped"; continue; }

      name="$(vcf_name "${platform}" "${RUNTIME_DEPTH}" "${caller}" "${phaser}")"
      # phase_one skips a VCF that exists: measuring means running it again
      rm -f "${PHASE_DIR}/${name}.vcf.gz" "${PHASE_DIR}/${name}.vcf.gz.tbi"
      rm -rf "${PHASE_DIR}/${name}.vcf.gz.lock"
      tfile="${LOGS_DIR}/time/${name}.txt"
      load1="$(cut -d' ' -f1 /proc/loadavg 2>/dev/null || echo NA)"
      log "${name}"
      set +e
      "${GNU_TIME}" -v -o "${tfile}" bash -c 'phase_one "$@"' _ \
          "${caller}" "${platform}" "${RUNTIME_DEPTH}" "${phaser}"
      status=$?
      set -e

      wall="$(to_seconds "$(time_value "${tfile}" 'Elapsed \\(wall clock\\)')")"
      user="$(time_value "${tfile}" 'User time')"
      sys="$(time_value "${tfile}" 'System time')"
      rss="$(time_value "${tfile}" 'Maximum resident set size')"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "${REFERENCE_NAME}" "${platform}" "${RUNTIME_DEPTH}" "${caller}" "${phaser}" "${PHASE_THREADS}" \
          "${wall}" "${user}" "${sys}" "$(awk -v k="${rss}" 'BEGIN { printf "%.0f", k / 1024 }')" \
          "${status}" "$(hostname -s)" "${load1}" "$(date '+%F %H:%M')" >>"${table}"
      [[ "${status}" -eq 0 ]] || warn "${name} failed (exit ${status}), see ${LOGS_DIR}/phase/${name}.log"
    done
  done
done

log "runtime done: ${table}"
