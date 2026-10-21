#!/usr/bin/env bash
# A GIAB sample other than HG002, from its source BAMs to every analysis: the
# same steps as HG002 (bin/run_reference.sh), on the platforms it has in
# samples.tsv, at 30x. GRCh38 by default; a sample has a T2T run only if
# samples.tsv gives it a truth_T2T and a bed_T2T.
#
#   SAMPLE=HG005 bash bin/run_sample.sh
#   SAMPLE=HG005 GPU_DEVICE=1 bash bin/run_sample.sh
#
# Check the sample first: bash bin/check_samples.sh HG005

set -euo pipefail

cd /autofs/bal36/dten/phasing_benchmark

[[ -n "${SAMPLE:-}" && "${SAMPLE}" != HG002 ]] \
  || { echo "set SAMPLE to a sample of samples.tsv other than HG002" >&2; exit 1; }
export SAMPLE
export REFERENCE_NAME="${REFERENCE_NAME:-GRCh38}"

exec bash bin/run_reference.sh
