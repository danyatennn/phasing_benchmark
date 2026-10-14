#!/usr/bin/env bash
# What each stratum actually is, so the per-stratum metrics can be read.
# A rate over 300 variants and over 300000 look the same without this table.
#
# metrics_strata/regions_statistics_<ref>.tsv   one row per stratum
# metrics_strata/jaccard_<ref>.tsv              pairwise overlap, long format
#
# Two sources of truth on purpose:
#   truth_het_snps   from derived/truth*.het_snps.phased.vcf.gz - the biallelic
#                    het SNPs, i.e. the denominator of every metric we report
#   src_*            from the raw GIAB VCF - the real composition of the stratum,
#                    including 1/2 and indels that our SNP-only benchmark drops.
# The gap between them says how much of the variation the benchmark is blind to.
#
# The overlap matrix covers JACCARD_TIERS (default the headline strata); strata
# are heavily nested, and without it they read as independent observations.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

JACCARD_TIERS="${JACCARD_TIERS:-1 2}"

mkdir -p "${METRICS_STRATA_DIR}" "${TMP_DIR}"

GENOME_BP="$(awk '{ s += $2 } END { print s }' "${CHROM_LENGTHS}")"
DIFFICULT_BED="$(stratum_bed alldifficultregions)"

stats_one() {
  set -euo pipefail

  local stratum="$1"
  local row="${METRICS_STRATA_DIR}/stats_${REFERENCE_NAME}_${stratum}.row.tsv"
  [[ -s "${row}" ]] && return 0

  local bed; bed="$(stratum_bed "${stratum}")"
  local truth; truth="$(stratum_truth "${stratum}")"
  [[ -s "${truth}.tbi" ]] || return 0

  local meta; meta="$(awk -F'\t' -v s="${stratum}" '$3 == s { print $1 "\t" $2 }' "${STRATA_LIST}")"

  # the stratum before it met the confident BED
  local src bp_raw intervals_raw
  src="$(find "${STRATA_DIR}/${STRATA_REF}@all" -name "${STRATA_REF}_${stratum}.bed.gz" | head -1)"
  read -r intervals_raw bp_raw < <(gzip -dc "${src}" | awk '{ n++; s += $3 - $2 } END { print n + 0, s + 0 }')

  # and after, which is what every metric is computed over
  local intervals bp median_interval n50_interval
  read -r intervals bp < <(awk '{ n++; s += $3 - $2 } END { print n + 0, s + 0 }' "${bed}")
  read -r median_interval n50_interval < <(awk '{ print $3 - $2 }' "${bed}" | sort -n | awk -v bp="${bp}" '
      { a[NR] = $1 }
      END {
        if (NR == 0) { print 0, 0; exit }
        med = (NR % 2) ? a[(NR + 1) / 2] : int((a[NR / 2] + a[NR / 2 + 1]) / 2)
        half = bp / 2
        for (i = NR; i >= 1; i--) { acc += a[i]; if (acc >= half) { n50 = a[i]; break } }
        print med, n50
      }')

  local truth_het_snps; truth_het_snps="$("${BCFTOOLS}" index -n "${truth}")"

  # composition from the raw truth: 0/1 against 1/2, SNP against indel
  local src_variants src_01 src_12 src_11 src_hemi src_other src_snp src_indel src_multi
  read -r src_variants src_01 src_12 src_11 src_hemi src_other src_snp src_indel src_multi < <(
    "${BCFTOOLS}" view -T "${bed}" "${TRUTH_SRC_VCF}" \
      | "${BCFTOOLS}" query -f '%N_ALT\t%TYPE\t[%GT]\n' \
      | awk '
          { n++
            if ($1 > 1) multi++
            if ($2 ~ /SNP/)   snp++
            if ($2 ~ /INDEL/) indel++
            gt = $3; gsub(/\|/, "/", gt); k = split(gt, g, "/")
            if (k < 2)                               hemi++   # chrX and chrY outside PAR
            else if (g[1] == "." || g[2] == ".")     other++   # missing
            else if ((g[1] == "0") != (g[2] == "0")) a01++     # 0/1, 1/0, 0|2
            else if (g[1] != g[2])                   a12++     # two different alts
            else if (g[1] != "0")                    a11++     # hom alt
            else                                     other++   # hom ref
          }
          END { print n + 0, a01 + 0, a12 + 0, a11 + 0, hemi + 0, other + 0, snp + 0, indel + 0, multi + 0 }')

  # distance to the next het SNP: phase blocks break where they are sparse
  local gaps="${TMP_DIR}/gaps_${REFERENCE_NAME}_${stratum}.$$"
  "${BCFTOOLS}" query -f '%CHROM\t%POS\n' "${truth}" \
    | awk 'NR > 1 && $1 == c { print $2 - p } { c = $1; p = $2 }' >"${gaps}"
  local median_gap frac_near
  median_gap="$(sort -n "${gaps}" | awk '{ a[NR] = $1 } END { print (NR ? ((NR % 2) ? a[(NR + 1) / 2] : int((a[NR / 2] + a[NR / 2 + 1]) / 2)) : 0) }')"
  frac_near="$(awk '{ n++; if ($1 < 50) k++ } END { printf "%.6f", (n ? k / n : 0) }' "${gaps}")"
  rm -f "${gaps}"

  local frac_difficult=0
  [[ -s "${DIFFICULT_BED}" && "${bp}" -gt 0 ]] && frac_difficult="$(
    bed_intersect "${DIFFICULT_BED}" "${bed}" \
      | awk -v bp="${bp}" '{ s += $3 - $2 } END { printf "%.6f", s / bp }')"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${REFERENCE_NAME}" "${stratum}" "${meta}" \
      "${intervals_raw}" "${bp_raw}" \
      "${intervals}" "${bp}" \
      "$(awk -v a="${bp}" -v b="${GENOME_BP}" 'BEGIN{printf "%.6f", a/b}')" \
      "${median_interval}" "${n50_interval}" \
      "${truth_het_snps}" \
      "$(awk -v n="${truth_het_snps}" -v bp="${bp}" 'BEGIN{printf "%.2f", (bp ? n / (bp / 1000000) : 0)}')" \
      "${src_variants}" "${src_01}" "${src_12}" "${src_11}" "${src_hemi}" "${src_other}" \
      "${src_snp}" "${src_indel}" "${src_multi}" \
      "${median_gap}" "${frac_near}" "${frac_difficult}" \
    >"${row}.$$"
  mv -f "${row}.$$" "${row}"
}

export -f stats_one stratum_bed stratum_truth bed_intersect log warn die
export METRICS_STRATA_DIR PREPARED_DIR STRATA_TRUTH_DIR STRATA_DIR STRATA_REF \
       STRATA_LIST REFERENCE_NAME TRUTH_SRC_VCF TMP_DIR GENOME_BP DIFFICULT_BED \
       BCFTOOLS

log "${REFERENCE_NAME}: statistics for $(strata_names | wc -l) strata"
strata_names | "${PARALLEL}" --joblog "${LOGS_DIR}/parallel_regions_stats_${REFERENCE_NAME}.log" \
    -j "${REGIONS_JOBS}" stats_one {} || true

out="${METRICS_STRATA_DIR}/regions_statistics_${REFERENCE_NAME}.tsv"
{
  printf 'reference\tstratum\ttier\tcategory'
  printf '\tintervals_raw\tbp_raw\tintervals\tbp\tfrac_genome\tmedian_interval\tn50_interval'
  printf '\ttruth_het_snps\thet_snps_per_mb'
  printf '\tsrc_variants\tsrc_0_1\tsrc_1_2\tsrc_1_1\tsrc_hemi\tsrc_other\tsrc_snp\tsrc_indel\tsrc_multiallelic'
  printf '\tmedian_inter_snp_dist\tfrac_snp_within_50bp\tfrac_in_alldifficult\n'
  cat "${METRICS_STRATA_DIR}/stats_${REFERENCE_NAME}"_*.row.tsv 2>/dev/null
} >"${out}.$$"
mv -f "${out}.$$" "${out}"
log "  ${out}"

# pairwise overlap, long format: stratum_a, stratum_b, shared bp, jaccard
jaccard_pair() {
  set -euo pipefail
  local a="$1" b="$2"
  local ba; ba="$(stratum_bed "${a}")"
  local bb; bb="$(stratum_bed "${b}")"
  [[ -s "${ba}" && -s "${bb}" ]] || return 0
  bed_intersect "${ba}" "${bb}" | awk -v a="${a}" -v b="${b}" -v OFS='\t' \
      -v bpa="$(bed_bp "${ba}")" -v bpb="$(bed_bp "${bb}")" '
      { s += $3 - $2 }
      END { u = bpa + bpb - s; print a, b, s + 0, (u ? sprintf("%.6f", s / u) : 0) }'
}
export -f jaccard_pair bed_bp

jac="${METRICS_STRATA_DIR}/jaccard_${REFERENCE_NAME}.tsv"
log "overlap matrix for tiers: ${JACCARD_TIERS}"
names="$(for t in ${JACCARD_TIERS}; do STRATA_TIER="${t}" strata_names; done)"
{
  printf 'stratum_a\tstratum_b\tshared_bp\tjaccard\n'
  for a in ${names}; do for b in ${names}; do
    [[ "${a}" < "${b}" ]] && printf '%s\t%s\n' "${a}" "${b}"
  done; done \
    | "${PARALLEL}" --colsep '\t' -j "${REGIONS_JOBS}" jaccard_pair {1} {2}
} >"${jac}.$$"
mv -f "${jac}.$$" "${jac}"
log "  ${jac}"
