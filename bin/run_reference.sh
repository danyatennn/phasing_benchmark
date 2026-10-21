#!/usr/bin/env bash
# One sample on one reference, from its sources to every analysis - the same
# steps on GRCh38 and T2T. Each step skips what is done, so a stopped run is
# resumed by starting it again.
#
#   bash bin/run_reference.sh                          HG002, GRCh38
#   REFERENCE_NAME=T2T bash bin/run_reference.sh
#   SAMPLE=HG005 bash bin/run_reference.sh             (its platforms, 30x)
#   RUN_ANALYSES=0 bash bin/run_reference.sh           the grid only
#
#   1. sources       GRCh38: the providers' BAMs; T2T: bin/align_to_t2t.sh
#   2. depth         bin/measure_sources.sh
#   3. downsampling  downsample.sh
#   4. truth         bin/make_truth.sh
#   5. truth input   phase_run.sh and eval_run.sh with the truth as input
#                    (no caller and no GPU needed, so it comes first)
#   6. callers       call_run.sh, then phase_run.sh and eval_run.sh
#   7. analyses      bin/run_analyses.sh
# The grid (platforms, depths, callers, phasers) is that of config.sh; the
# PHASING_* variables narrow it, as for every script.

set -euo pipefail

cd /autofs/bal36/dten/phasing_benchmark

export SAMPLE="${SAMPLE:-HG002}"
export REFERENCE_NAME="${REFERENCE_NAME:-GRCh38}"
RUN_ANALYSES="${RUN_ANALYSES:-1}"

step() {
  printf '[%s] ===== %s %s: %s\n' "$(date '+%F %H:%M:%S')" "${SAMPLE}" "${REFERENCE_NAME}" "$*" >&2
  "$@"
}

[[ "${REFERENCE_NAME}" == GRCh38 ]] || step bash bin/align_to_t2t.sh
step bash bin/measure_sources.sh
step bash downsample.sh
step bash bin/make_truth.sh
step env PHASING_CALLERS=truth bash phase_run.sh
step env PHASING_CALLERS=truth bash eval_run.sh
step bash call_run.sh
step bash phase_run.sh
step bash eval_run.sh
[[ "${RUN_ANALYSES}" == 0 ]] || step bash bin/run_analyses.sh
printf '[%s] ===== %s %s: done\n' "$(date '+%F %H:%M:%S')" "${SAMPLE}" "${REFERENCE_NAME}" >&2
