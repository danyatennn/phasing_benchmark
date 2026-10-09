#!/usr/bin/env bash
set -euo pipefail

cd /autofs/bal36/dten/phasing_benchmark

RUN_REFERENCES="${RUN_REFERENCES:-GRCh38 T2T}"

export PHASE_JOBS="${PHASE_JOBS:-12}"
export EVAL_JOBS="${EVAL_JOBS:-12}"

for reference in ${RUN_REFERENCES}; do
  export REFERENCE_NAME="${reference}"

  d="logs/joblog_${reference}_$(date +%F_%H%M)"
  mkdir -p "${d}" logs
  mv logs/parallel_*.log "${d}"/ 2>/dev/null || true

  echo "######## ${reference} ########"
  bash bin/make_truth.sh
  bash call_run.sh
  bash phase_run.sh
  bash eval_run.sh
done

echo "wgs grid done: ${RUN_REFERENCES}"
