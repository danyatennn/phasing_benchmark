#!/usr/bin/env bash

set -euo pipefail

cd /autofs/bal36/dten/phasing_benchmark
mkdir -p logs

stage() {
  log_line "=== $* ==="
  "$@" 2>&1 | tee -a logs/run_all.log
}

log_line() { printf '[%s] %s\n' "$(date '+%F %H:%M:%S')" "$*" | tee -a logs/run_all.log; }

log_line "start"

stage bash call_run.sh
stage bash phase_run.sh
stage bash eval_run.sh

stage bash assembly_run.sh
stage bash dipcall_run.sh
PHASING_CALLERS=hifiasm PHASING_PHASERS=dipcall stage bash eval_run.sh

log_line "all done"
