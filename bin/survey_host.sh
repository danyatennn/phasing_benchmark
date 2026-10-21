#!/usr/bin/env bash

printf 'host      %s\n' "$(hostname)"
printf 'ядер      %s потоков, %s физических\n' "$(nproc)" \
    "$(lscpu 2>/dev/null | awk -F: '/^Core\(s\) per socket/ {c=$2} /^Socket\(s\)/ {s=$2} END {print (c && s) ? c * s : "?"}' | tr -d ' ')"
printf 'память    %s ГБ всего, %s ГБ свободно\n' \
    "$(free -g | awk 'NR == 2 { print $2 }')" "$(free -g | awk 'NR == 2 { print $7 }')"

if command -v nvidia-smi >/dev/null 2>&1; then
  printf 'GPU       %s\n' "$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | paste -sd'; ' -)"
  busy="$(nvidia-smi --query-compute-apps=process_name,used_memory --format=csv,noheader | paste -sd'; ' -)"
  printf 'GPU занят %s\n' "${busy:-свободна}"
else
  printf 'GPU       нет nvidia-smi\n'
fi

printf '\n-- куда распаковывать SIF (нужен локальный диск, ~30 ГБ) --\n'
for p in /tmp /var/tmp /scratch /local /data "${HOME}"; do
  [[ -d "${p}" ]] || continue
  line="$(df -hT "${p}" 2>/dev/null | awk 'NR == 2 { print $2, $3, $5 }')" || continue
  set -- ${line}
  case "$1" in
    ext2|ext3|ext4|xfs|btrfs|tmpfs|overlay) kind="локальный" ;;
    nfs|nfs4|autofs|fuse*)                  kind="СЕТЕВОЙ" ;;
    *)                                      kind="$1 ?" ;;
  esac
  printf '  %-12s %-10s %-9s всего %-7s свободно %s\n' "${p}" "$1" "${kind}" "$2" "$3"
done

for p in \
    /autofs/bal36/dten/phasing_benchmark \
    /autofs/bal36/dten/phasing_benchmark/data \
    /autofs/bal19/zxzheng/testData/ont/data \
    /autofs/bal33/zxzheng/data/vcf/v5.0q \
    /autofs/bal19/lchen/data/chm13/ref/chm13v2.0.fa \
    /autofs/bal33/zxzheng/env/miniconda2/envs/singularity-env/bin/singularity \
    /autofs/bal36/dten/miniforge3/envs/phasing/bin/bcftools \
    /autofs/bal36/dten/miniforge3/envs/hapcut2/bin/HAPCUT2 \
    /autofs/bal36/dten/miniforge3/envs/asmeval/bin/run-dipcall ; do
  if [[ -e "${p}" ]]; then printf 'есть      %s\n' "${p}"
  else                     printf 'НЕТ       %s\n' "${p}"; fi
done

printf '\n-- запускаются ли бинарники на этом процессоре и этой системе --\n'
smoke() {
  local label="$1"; shift
  local all rc first
  if all="$("$@" 2>&1)"; then rc=0; else rc=$?; fi
  # these tools print their usage to stderr and exit non-zero when called with
  # no real arguments, and some start it with a blank line - take the first line
  # that actually says something, and show rc so a banner is not mistaken for a
  # failure
  first="$(printf '%s\n' "${all}" | awk 'NF { print; exit }')"
  case "${rc}" in
    126|127) printf '  НЕ ЗАПУСКАЕТСЯ %-14s rc=%-3s %s\n' "${label}" "${rc}" "${first}"; return ;;
    132)     printf '  SIGILL         %-14s собран под более новый процессор\n' "${label}"; return ;;
  esac
  case "${all}" in
    *"error while loading shared libraries"*|*GLIBC*|*"Illegal instruction"*|*"No such file or directory"*)
             printf '  БИБЛИОТЕКИ     %-14s rc=%-3s %s\n' "${label}" "${rc}" "${first}" ;;
    *)       printf '  ok             %-14s rc=%-3s %s\n' "${label}" "${rc}" "${first}" ;;
  esac
}

if [[ -r /autofs/bal36/dten/phasing_benchmark/config.sh ]]; then
  source /autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh
  source /autofs/bal36/dten/phasing_benchmark/lib.sh
  source /autofs/bal36/dten/phasing_benchmark/config.sh

  smoke bcftools   "${BCFTOOLS}" --version
  smoke samtools   "${SAMTOOLS}" --version
  smoke whatshap   "${WHATSHAP}" --version
  smoke parallel   "${PARALLEL}" --version
  smoke mosdepth   "${MOSDEPTH}" --version
  smoke minimap2   "${MINIMAP2}" --version
  smoke singularity "${SING}" --version
  smoke longphase  "${LONGPHASE}" --help
  smoke longcallD  "${LONGCALLD}"
  smoke margin     "${MARGIN}" --help
  smoke HAPCUT2    "${HAPCUT2_BIN}/HAPCUT2" --help
  smoke extractHAIRS "${HAPCUT2_BIN}/extractHAIRS"
  smoke hifiasm    "${HIFIASM}" --version
  smoke dipcall    "${DIPCALL_BIN}/run-dipcall"
else
  printf '  проект не виден, проверить нечем\n'
fi
