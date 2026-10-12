#!/usr/bin/env bash

#   watch -n 300 'bash bin/progress.sh'

set -uo pipefail

RUN_REFERENCES="${RUN_REFERENCES:-GRCh38 T2T}"

printf '%s\n' "$(date '+%F %H:%M')"

for reference in ${RUN_REFERENCES}; do
(
  export REFERENCE_NAME="${reference}"
  source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
  source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
  source "/autofs/bal36/dten/phasing_benchmark/config.sh"
  cd "${PROJECT_ROOT}"

  cells=$(( ${#PLATFORMS[@]} * ${#DEPTHS[@]} ))
  want_calls=$(( (${#CALLERS[@]} - 1) * cells ))
  want_runs=$(( ${#CALLERS[@]} * cells * ${#PHASERS[@]} ))

  count() { ls $1 2>/dev/null | wc -l | tr -d ' '; }

  printf '%-7s вызовы %3s/%-3s  фазирования %3s/%-3s  строки %3s/%-3s\n' \
      "${REFERENCE_NAME}" \
      "$(count "vcf/*_${REFERENCE_NAME}_*_${REGION_NAME}_*.het_snps.vcf.gz")" "${want_calls}" \
      "$(count "phased/*_${REFERENCE_NAME}_*_${REGION_NAME}_*.vcf.gz")"       "${want_runs}" \
      "$(count "metrics/*_${REFERENCE_NAME}_*_${REGION_NAME}_*.row.tsv")"     "${want_runs}"
)
done

cd /autofs/bal36/dten/phasing_benchmark

printf '\nсейчас работает:\n'
for b in clair3 run_deepvariant longcallD whatshap longphase margin HAPCUT2 extractHAIRS; do
  n="$(pgrep -c -f "${b}" 2>/dev/null)" || n=0
  [[ "${n:-0}" -gt 0 ]] && printf '  %-16s %s\n' "${b}" "${n}"
done

printf '\nсбои в joblog (7-я колонка - код выхода):\n'
for f in logs/parallel_*.log; do
  [[ -s "${f}" ]] || continue
  awk -v f="${f}" 'NR > 1 { n++; if ($7 != 0) bad++ }
    END { if (n) printf "  %-32s задач=%-4d сбоев=%d\n", f, n, bad + 0 }' "${f}"
done

printf '\nпамять: '
free -g | awk 'NR == 2 { printf "занято %s ГБ, свободно %s ГБ\n", $3, $7 }'

printf 'GPU:    '
nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader
nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader \
  | sed 's|^|        |'
