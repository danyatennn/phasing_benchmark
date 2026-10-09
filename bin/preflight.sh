#!/usr/bin/env bash
# Would anything left over from an earlier run get in the way? Read-only:
# prints a report and changes nothing.
#
#   bash bin/preflight.sh      checks every reference run_wgs.sh will do

# no -e, so one failing check does not cut the report short
set -uo pipefail

RUN_REFERENCES="${RUN_REFERENCES:-GRCh38 T2T}"

for reference in ${RUN_REFERENCES}; do
(
  export REFERENCE_NAME="${reference}"
  source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
  source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
  source "/autofs/bal36/dten/phasing_benchmark/config.sh"
  cd "${PROJECT_ROOT}"

  runs=$(( ${#CALLERS[@]} * ${#PLATFORMS[@]} * ${#DEPTHS[@]} * ${#PHASERS[@]} ))
  calls=$(( (${#CALLERS[@]} - 1) * ${#PLATFORMS[@]} * ${#DEPTHS[@]} ))

  printf '\n################ %s / %s ################\n' "${REFERENCE_NAME}" "${REGION_NAME}"
  printf 'платформы %s\nглубины   %s\nколлеры   %s\nфазеры    %s\nждём      %s вызовов, %s фазирований, %s строк\n' \
      "${PLATFORMS[*]}" "${DEPTHS[*]}" "${CALLERS[*]}" "${PHASERS[*]}" "${calls}" "${runs}" "${runs}"

  printf '\n-- референс, BED, истина --\n'
  for v in REF CONFIDENT_BED TRUTH_SRC_VCF TRUTH_PHASED; do
    printf '  %-15s %s\n' "${v}" "${!v}"
    [[ -s "${!v}" ]] || printf '  %-15s ^^^ ФАЙЛА НЕТ\n' ''
  done
  [[ -s "${REF}.fai" ]] || printf '  нет %s.fai - eval не построит chrom_lengths\n' "${REF}"
  if [[ -s "${TRUTH_PHASED}" ]]; then
    printf '  истина: %s сайтов, %s хромосом\n' \
        "$("${BCFTOOLS}" index -n "${TRUTH_PHASED}")" \
        "$("${BCFTOOLS}" index -s "${TRUTH_PHASED}" | wc -l)"
  else
    printf '  истина соберётся первым шагом run_wgs.sh\n'
  fi

  printf '\n-- готовые вызовы, которые запуск пропустит --\n'
  found=0
  for caller in "${CALLERS[@]}"; do
    [[ "${caller}" == truth ]] && continue
    for p in "${PLATFORMS[@]}"; do
      for dep in "${DEPTHS[@]}"; do
        name="$(vcf_name "${p}" "${dep}" "${caller}")"
        v="${VCF_DIR}/${name}.vcf.gz"
        [[ -s "${v}" ]] || continue
        found=$((found + 1))
        cnt="$("${BCFTOOLS}" index -n "${v}" 2>/dev/null)" || cnt=СЛОМАН
        tbi=НЕТ; [[ -s "${v}.tbi" ]] && tbi=ok
        het=НЕТ; [[ -s "${VCF_DIR}/${name}.het_snps.vcf.gz" ]] && het=ok
        printf '  %-56s вариантов=%-9s tbi=%-3s het_snps=%s\n' "${name}" "${cnt}" "${tbi}" "${het}"
      done
    done
  done
  [[ "${found}" -eq 0 ]] && printf '  нет, всё посчитается с нуля\n'

  printf '\n-- готовые фазирования, которые запуск пропустит --\n'
  found=0
  for caller in "${CALLERS[@]}"; do
    for p in "${PLATFORMS[@]}"; do
      for dep in "${DEPTHS[@]}"; do
        for ph in "${PHASERS[@]}"; do
          name="$(vcf_name "${p}" "${dep}" "${caller}" "${ph}")"
          v="${PHASE_DIR}/${name}.vcf.gz"
          [[ -s "${v}" ]] || continue
          found=$((found + 1))
          cnt="$("${BCFTOOLS}" index -n "${v}" 2>/dev/null)" || cnt=СЛОМАН
          tbi=НЕТ; [[ -s "${v}.tbi" ]] && tbi=ok
          row=НЕТ; [[ -s "${METRICS_DIR}/${name}.row.tsv" ]] && row=ok
          printf '  %-60s вариантов=%-9s tbi=%-3s row=%s\n' "${name}" "${cnt}" "${tbi}" "${row}"
        done
      done
    done
  done
  [[ "${found}" -eq 0 ]] && printf '  нет, всё посчитается с нуля\n'
)
done

cd /autofs/bal36/dten/phasing_benchmark

printf '\n################ общее ################\n'

printf '\n-- что лежит в vcf/ phased/ metrics/, по регионам --\n'
for d in vcf phased metrics; do
  printf '  %s/\n' "${d}"
  ls "${d}" 2>/dev/null \
    | awk -F_ 'NF >= 5 { c[$5]++ } END { for (r in c) printf "    region=%-8s %d файлов\n", r, c[r] }'
done

printf '\n-- строки метрик: 20 полей значит строка от старой версии eval_run.sh --\n'
awk -F'\t' 'FNR == 1 { c[NF]++ } END { for (n in c) printf "  %s полей: %s файлов\n", n, c[n] }' \
    metrics/*.row.tsv 2>/dev/null || printf '  строк пока нет\n'

printf '\n-- место --\n'
df -h . | tail -1
du -sh data derived vcf phased metrics asm tmp logs 2>/dev/null
printf '  каталогов в tmp/: %s\n' "$(ls -d tmp/*/ 2>/dev/null | wc -l)"
