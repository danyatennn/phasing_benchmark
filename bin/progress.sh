#!/usr/bin/env bash

#   watch -n 300 'bash bin/progress.sh'

set -uo pipefail

RUN_REFERENCES="${RUN_REFERENCES:-GRCh38 T2T}"

printf '%s   %s\n' "$(hostname)" "$(date '+%F %H:%M')"

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

printf '\nсейчас работает (clair3 внутри себя тоже зовёт whatshap и longphase):\n'
for b in clair3 run_deepvariant longcallD whatshap longphase margin HAPCUT2 extractHAIRS; do
  n="$(pgrep -c -f "${b}" 2>/dev/null)" || n=0
  [[ "${n:-0}" -gt 0 ]] && printf '  %-16s %s\n' "${b}" "${n}"
done

# logs/ is shared by every host, so without a time window this lists the dead
# slices forever. Backups made by keep_joblog are skipped too - same work twice.
printf '\njoblog за последние 6 часов (7-я колонка - код выхода):\n'
shown=0
while IFS= read -r f; do
  short="$(basename "${f}" | sed 's/^parallel_//; s/_10-20-30-40-50//; s/\.log$//')"
  awk -v f="${short}" 'NR > 1 { n++; if ($7 != 0) bad++ }
    END { if (n) printf "  %-52s задач=%-4d сбоев=%d\n", f, n, bad + 0 }' "${f}"
  shown=1
done < <(find logs -maxdepth 1 -name 'parallel_*.log' \
    ! -name '*_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9][0-9][0-9].log' \
    -mmin -360 2>/dev/null | sort)
[[ "${shown}" -eq 1 ]] || printf '  нет свежих\n'

printf 'всего сбоев за всю историю: %s\n' \
    "$(find logs -maxdepth 1 -name 'parallel_*.log' ! -name '*_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9][0-9][0-9].log' -exec awk 'FNR > 1 && $7 != 0 { n++ } END { print n + 0 }' {} + 2>/dev/null | awk '{ s += $1 } END { print s + 0 }')"

# a lock whose owner process is gone makes that cell skip forever, and only the
# host named in the lock can tell whether the pid is still alive
printf '\nзамки с умершим процессом этого хоста:\n'
stale=0
anon=0
for l in vcf/*.lock phased/*.lock metrics/*.lock; do
  if [[ ! -f "${l}/owner" ]]; then
    [[ -d "${l}" ]] && anon=$((anon + 1))
    continue
  fi
  read -r lhost lpid _ <"${l}/owner"
  [[ "${lhost}" == "$(hostname)" ]] || continue
  kill -0 "${lpid}" 2>/dev/null && continue
  printf '  %s  (%s)\n' "${l}" "$(cat "${l}/owner")"
  stale=$((stale + 1))
done
[[ "${stale}" -eq 0 ]] && printf '  нет\n'
[[ "${anon}" -eq 0 ]] || printf '  плюс %s замков без подписи - от версии до правки, судить по времени\n' "${anon}"

printf '\nпоследние ошибки из логов всех хостов:\n'
grep -hE 'ERROR|WARN' logs/run_*.log 2>/dev/null | tail -5 | sed 's|^|  |'

printf '\nпамять: '
free -g | awk 'NR == 2 { printf "занято %s ГБ, свободно %s ГБ\n", $3, $7 }'

printf 'место:  проект %s свободно' "$(df -h . | awk 'NR == 2 { print $4 }')"
if [[ -n "${SCRATCH_DIR:-}" && "${SCRATCH_DIR}" != "${TMP_DIR:-}" ]]; then
  printf ', scratch %s свободно' "$(df -h "${SCRATCH_DIR}" 2>/dev/null | awk 'NR == 2 { print $4 }')"
fi
printf '\n'

if command -v nvidia-smi >/dev/null 2>&1; then
  printf 'GPU:    '
  nvidia-smi --query-gpu=index,memory.used,memory.total --format=csv,noheader | paste -sd' | ' -
  nvidia-smi --query-compute-apps=gpu_uuid,pid,used_memory --format=csv,noheader \
    | sed 's|^|        |'
else
  printf 'GPU:    нет\n'
fi
