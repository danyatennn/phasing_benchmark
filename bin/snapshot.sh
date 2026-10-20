#!/usr/bin/env bash
# A copy of every result table: the per-run rows and read accuracy, the summaries, the stratum tables
# bash bin/snapshot.sh -> snapshots/<date>.tar.gz, plus a git tag

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

cd "${PROJECT_ROOT}"

stamp="$(date +%Y-%m-%d_%H%M)"
out="snapshots/${stamp}.tar.gz"
mkdir -p snapshots

list="$(mktemp)"
trap 'rm -f "${list}"' EXIT

for root in . samples/*; do
  [[ -d "${root}" ]] || continue
  for d in metrics metrics_reads; do
    [[ -d "${root}/${d}" ]] || continue
    find "${root}/${d}" -maxdepth 1 -type f \( \
        -name '*.row.tsv' -o -name '*.readacc*.tsv' -o -name '*.ng50.tsv' \
        -o -name '*.tagbases*.tsv' -o -name '*.auto.tsv' -o -name 'summary*.tsv' \
        -o -name 'calls_*.tsv' -o -name 'calibration.tsv' -o -name '*.classes.tsv' \
        -o -name '*.anchor_*.tsv' \) \
      -print
  done
  for d in metrics_strata metrics_strata_auto; do
    if [[ -d "${root}/${d}" ]]; then
      find "${root}/${d}" -maxdepth 1 -type f -name '*.tsv' -print
    fi
  done
done | sed 's|^\./||' | sort >"${list}"

n="$(wc -l <"${list}" | tr -d ' ')"
[[ "${n}" -gt 0 ]] || die "nothing to snapshot under ${PROJECT_ROOT}"

tar czf "${out}.tmp" -T "${list}"
mv -f "${out}.tmp" "${out}"
log "snapshot: ${n} files -> ${PROJECT_ROOT}/${out} ($(du -h "${out}" | cut -f1))"

if git -C "${PROJECT_ROOT}" rev-parse --git-dir >/dev/null 2>&1; then
  tag="snapshot-${stamp}"
  git -C "${PROJECT_ROOT}" tag "${tag}" \
    && log "git tag ${tag} on $(git -C "${PROJECT_ROOT}" rev-parse --short HEAD)" \
    || warn "could not create git tag ${tag}"
fi
