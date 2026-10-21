#!/usr/bin/env bash
# HG002 on every reference: bin/run_reference.sh for each of RUN_REFERENCES.
#
#   bash bin/run_wgs.sh
#   RUN_REFERENCES=T2T bash bin/run_wgs.sh
#   RUN_ANALYSES=0 bash bin/run_wgs.sh       the grid only

set -euo pipefail

cd /autofs/bal36/dten/phasing_benchmark

RUN_REFERENCES="${RUN_REFERENCES:-GRCh38 T2T}"
export PHASE_JOBS="${PHASE_JOBS:-12}"
export EVAL_JOBS="${EVAL_JOBS:-12}"

for reference in ${RUN_REFERENCES}; do
  REFERENCE_NAME="${reference}" bash bin/run_reference.sh
done
echo "wgs grid done: ${RUN_REFERENCES}"
