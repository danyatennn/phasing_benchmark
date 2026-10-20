#!/usr/bin/env bash
# Lists every metrics/<run>.row.tsv that eval_run.sh would not accept any more
# (row_is_valid in lib.sh) - and so would recompute. Read-only.
#
#   bash bin/check_rows.sh                 HG002
#   SAMPLE=HG005 bash bin/check_rows.sh

set -uo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"
set +e

total=0
bad=0
for row in "${METRICS_DIR}"/*.row.tsv; do
  [[ -e "${row}" ]] || continue
  total=$((total + 1))
  if ! row_is_valid "${row}"; then
    bad=$((bad + 1))
    printf '%-70s %s fields\n' "$(basename "${row}")" "$(awk -F'\t' 'NR == 1 { print NF }' "${row}")"
  fi
done
printf '%s: %d rows, %d would be recomputed\n' "${METRICS_DIR}" "${total}" "${bad}"
