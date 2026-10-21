#!/usr/bin/env bash
# Checks the inputs listed in samples.tsv before a sample is run. Read-only:
# prints a report and changes nothing.
#
#   bash bin/check_samples.sh               every sample in samples.tsv
#   bash bin/check_samples.sh HG005 HG006   these only
#
# For each sample:
#   BAMs   present, indexed as <bam>.bai (downsample.sh wants exactly that),
#          the same chr1-22 names and lengths as REF_GRCH38, haplotype tags
#          (HP) on them (downsample.sh strips them), the aligner in @PG
#   truth  present, its sample name, het SNPs on chr1-22, and how many of them
#          are phased - only phased sites enter the truth (make_truth.sh -p)
#   BED    present, its size on chr1-22
#   parents' reads (pat_reads, mat_reads) present

set -uo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"
set +e

[[ -s "${SAMPLES_TSV}" ]] || { echo "no ${SAMPLES_TSV}"; exit 1; }
FAI="${REF_GRCH38}.fai"
[[ -s "${FAI}" ]] || { echo "no ${FAI}"; exit 1; }

AUTO='^chr([1-9]|1[0-9]|2[0-2])$'

samples=("$@")
[[ ${#samples[@]} -gt 0 ]] || mapfile -t samples < <(awk -F'\t' '!/^#/ && NF >= 3 { print $1 }' "${SAMPLES_TSV}" | awk '!seen[$0]++')

flag() { printf '  %-4s %s\n' "$1" "$2"; }

check_bam() {
  local bam="$1"
  [[ -s "${bam}" ]] || { flag NO "file missing: ${bam}"; return; }
  if [[ -s "${bam}.bai" ]]; then
    flag ok "index ${bam}.bai"
  elif [[ -s "${bam%.bam}.bai" || -s "${bam}.csi" ]]; then
    flag NO "index is not <bam>.bai (downsample.sh needs that name): $(ls "${bam%.bam}.bai" "${bam}.csi" 2>/dev/null | head -1)"
  else
    flag NO "no index"
  fi
  local diff
  diff="$(comm -3 \
      <("${SAMTOOLS}" view -H "${bam}" | awk -F'\t' -v re="${AUTO}" '$1 == "@SQ" {
          for (i = 2; i <= NF; i++) { if ($i ~ /^SN:/) n = substr($i, 4); if ($i ~ /^LN:/) l = substr($i, 4) }
          if (n ~ re) print n "\t" l }' | sort) \
      <(awk -F'\t' -v re="${AUTO}" '$1 ~ re { print $1 "\t" $2 }' "${FAI}" | sort) | wc -l)"
  [[ "${diff}" -eq 0 ]] && flag ok "chr1-22 match the reference" \
                        || flag NO "chr1-22 differ from $(basename "${REF_GRCH38}") in ${diff} lines: another reference build"
  local hp
  hp="$("${SAMTOOLS}" view "${bam}" chr20:30000000-30200000 2>/dev/null | head -200 | grep -c $'\tHP:i:')"
  [[ "${hp}" -gt 0 ]] && flag note "HP tags on ${hp}/200 reads (stripped by downsample.sh)" \
                      || flag ok "no HP tags"
  "${SAMTOOLS}" view -H "${bam}" | grep '^@PG' | grep -m2 -iE 'minimap2|pbmm2|winnowmap|ngmlr|dorado|ccs' \
    | cut -c1-140 | sed 's/^/       /'
  # ONT: the basecaller model says the chemistry (r1041_e82_400bps = R10.4.1,
  # 4 kHz unless _5khz / v5 models; r941 = R9.4.1)
  "${SAMTOOLS}" view -H "${bam}" | grep -oE 'basecall_model=[^[:space:]]+|(dna_r[0-9][^[:space:],]*)' \
    | sort -u | head -3 | sed 's/^/       model: /'
}

check_truth() {
  local vcf="$1"
  [[ -s "${vcf}" ]] || { flag NO "file missing: ${vcf}"; return; }
  flag ok "sample name: $("${BCFTOOLS}" query -l "${vcf}" | tr '\n' ' ')"
  "${BCFTOOLS}" view -H -v snps -m2 -M2 -g het "${vcf}" 2>/dev/null \
    | awk -F'\t' -v re="${AUTO}" '$1 ~ re { n++; split($10, a, ":"); if (a[1] ~ /\|/) p++ }
        END { printf "  %-4s het biallelic SNPs on chr1-22: %d, phased %d (%.1f%%)\n", "ok", n, p, (n ? 100 * p / n : 0) }'
}

check_bed() {
  local bed="$1"
  [[ -s "${bed}" ]] || { flag NO "file missing: ${bed} - candidates: $(ls "$(dirname "${bed}")" 2>/dev/null | grep -i "$(basename "${bed}" | cut -d_ -f1).*bed" | tr '\n' ' ')"; return; }
  awk -F'\t' -v re="${AUTO}" '$1 ~ re { n++; s += $3 - $2 } END { printf "  %-4s chr1-22: %d intervals, %.1f Mb\n", "ok", n, s / 1e6 }' "${bed}"
}

for s in "${samples[@]}"; do
  printf '\n######## %s ########\n' "${s}"
  # a value $NAME of samples.tsv is a data_config variable: an unset one would
  # drop its platform without a word, so it is reported here
  while IFS=$'\t' read -r item raw; do
    name="${raw#\$}"; name="${name#\{}"; name="${name%\}}"
    if [[ -n "${!name:-}" ]]; then
      printf '  %-14s %s\n' "${item}" "${raw}"
    else
      flag NO "${item}: ${raw} is not set in data_config_germline.sh - this input is skipped"
    fi
  done < <(awk -F'\t' -v s="${s}" '!/^#/ && $1 == s && $3 ~ /^\$/ { print $2 "\t" $3 }' "${SAMPLES_TSV}")
  for p in ont4k ont5k hifi revio; do
    bam="$(sample_path "${s}" "bam_${p}")"
    [[ -n "${bam}" ]] || continue
    printf '%s  %s\n' "${p}" "${bam}"
    check_bam "${bam}"
  done
  truth="$(sample_path "${s}" truth_GRCh38)"
  [[ -n "${truth}" ]] && { printf 'truth  %s\n' "${truth}"; check_truth "${truth}"; }
  bed="$(sample_path "${s}" bed_GRCh38)"
  [[ -n "${bed}" ]] && { printf 'bed    %s\n' "${bed}"; check_bed "${bed}"; }
  for role in pat_reads mat_reads; do
    for f in $(sample_path "${s}" "${role}"); do
      [[ -s "${f}" ]] && flag ok "${role}: ${f}" || flag NO "${role} missing: ${f}"
    done
  done
done

# the ONT 30x BAMs of data_config whose chemistry is not known yet: candidates
# for bam_ont4k of HG001 and HG005-HG007 (commented out in samples.tsv)
printf '\n######## ONT 30x of unknown chemistry ########\n'
for v in HG001_ONT_30x_GRCH38_BAM HG005_ONT_30x_GRCH38_BAM HG006_ONT_30x_GRCH38_BAM HG007_ONT_30x_GRCH38_BAM; do
  [[ -n "${!v:-}" ]] || continue
  printf '%s  %s\n' "${v}" "${!v}"
  check_bam "${!v}"
done

