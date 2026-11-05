#!/usr/bin/env bash
# Phased truth for a GIAB sample from its trio assemblies: the het SNPs of the
# sample's GIAB v4.2.1 benchmark, phased by the paternal and maternal haplotype
# assemblies (hifiasm-trio, HPRC release 2), keeping only the variants the
# assemblies phase with confidence.
#
#   1. dipcall of the two haplotypes against the reference, hap1 = paternal,
#      hap2 = maternal, so a call's GT is paternal|maternal
#   2. regions: the v4.2.1 benchmark BED x dip.bed (each haplotype has exactly
#      one long alignment there), chr1-22
#   3. truth: the v4.2.1 biallelic het SNPs in those regions that dipcall calls
#      as the same het SNP (same position, REF and ALT), with dipcall's phase.
#      A v4.2.1 het SNP where the assemblies carry the same allele, or another
#      one, is not phasable from them and is left out
#   4. check: the counts at each step, and whatshap compare of the new phase
#      against the phase of the v4.2.1 VCF where it has one
#
#   <ASM_TRUTH_DIR>/<sample>_<ref>_v4.2.1_asmphased.vcf.gz   the truth (GT pat|mat)
#   <ASM_TRUTH_DIR>/<sample>_<ref>_v4.2.1_asmphased.bed      its regions
#   <ASM_TRUTH_DIR>/check.txt                                counts and the compare
#   <dir of ASM_PAT>/<sample>.dip.fa.gz   both assemblies in one FASTA, contigs
#                     <name>_PATERNAL / <name>_MATERNAL: the diploid_fasta of
#                     samples.tsv, for bin/matpat_truth.sh and consensus_truth.sh
#
# Both are then given to a sample of their own in samples.tsv (truth_<ref>,
# bed_<ref>), so that the runs against the v4.2.1 phase stay as they are.
#
#   SAMPLE=HG005 ASM_PAT=data/hg005_hprc_r2/HG005_pat_hprc_r2_v1.0.1.fa.gz \
#     ASM_MAT=data/hg005_hprc_r2/HG005_mat_hprc_r2_v1.0.1.fa.gz bash bin/asm_truth.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

[[ "${SAMPLE}" != HG002 ]] || die "HG002 has its Q100 truth; this is for the samples with a v4.2.1 one"
[[ -n "${ASM_PAT:-}" && -n "${ASM_MAT:-}" ]] || die "set ASM_PAT and ASM_MAT, the paternal and maternal assemblies"
[[ "${ASM_PAT}" == /* ]] || ASM_PAT="${PROJECT_ROOT}/${ASM_PAT}"
[[ "${ASM_MAT}" == /* ]] || ASM_MAT="${PROJECT_ROOT}/${ASM_MAT}"
require_file "${ASM_PAT}" "paternal assembly"
require_file "${ASM_MAT}" "maternal assembly"
require_file "${TRUTH_SRC_VCF}" "v4.2.1 truth of ${SAMPLE} (samples.tsv)"
require_file "${CONFIDENT_BED}" "v4.2.1 benchmark BED of ${SAMPLE} (samples.tsv)"

ASM_TRUTH_DIR="${ASM_TRUTH_DIR:-$(dirname "${ASM_PAT}")/truth_${REFERENCE_NAME}}"
base="${ASM_TRUTH_DIR}/${SAMPLE}_${REFERENCE_NAME}_v4.2.1_asmphased"
work="${ASM_TRUTH_DIR}/work"
pre="${work}/dip"
check="${ASM_TRUTH_DIR}/check.txt"
logf="${ASM_TRUTH_DIR}/dipcall.log"
mkdir -p "${work}"

# the diploid FASTA for the read labels, named as the Q100 copies are
dip_fa="$(dirname "${ASM_PAT}")/${SAMPLE}.dip.fa.gz"
if ! newer_than "${dip_fa}" "${ASM_PAT}" || ! newer_than "${dip_fa}" "${ASM_MAT}"; then
  log "${SAMPLE}: $(basename "${dip_fa}"), paternal and maternal contigs in one FASTA"
  { gzip -cdf "${ASM_PAT}" | sed 's/^>\([^[:space:]]*\).*/>\1_PATERNAL/'
    gzip -cdf "${ASM_MAT}" | sed 's/^>\([^[:space:]]*\).*/>\1_MATERNAL/'; } \
    | gzip -c >"${dip_fa}.tmp"
  mv -f "${dip_fa}.tmp" "${dip_fa}"
  log "${SAMPLE}: $(gzip -dc "${dip_fa}" | grep -c '^>') contigs in ${dip_fa}"
fi

if newer_than "${base}.vcf.gz" "${ASM_PAT}" && newer_than "${base}.vcf.gz" "${ASM_MAT}" \
    && newer_than "${base}.vcf.gz" "${TRUTH_SRC_VCF}" && [[ -s "${base}.bed" ]]; then
  log "${SAMPLE} ${REFERENCE_NAME}: assembly-phased truth kept: ${base}.vcf.gz"
  cat "${check}" >&2
  exit 0
fi

# ---------------------------------------------------------------- 1. dipcall
if [[ ! -s "${pre}.dip.vcf.gz" || ! -s "${pre}.dip.bed" ]]; then
  for h in pat mat; do
    src="ASM_${h^^}"
    [[ -s "${work}/${h}.fa" ]] || { gzip -cdf "${!src}" >"${work}/${h}.fa.tmp" && mv -f "${work}/${h}.fa.tmp" "${work}/${h}.fa"; }
  done
  log "${SAMPLE}: dipcall, hap1 = paternal, hap2 = maternal, against $(basename "${REF}")"
  "${DIPCALL_BIN}/run-dipcall" -t "${DIPCALL_THREADS}" "${pre}" "${REF}" "${work}/pat.fa" "${work}/mat.fa" \
      >"${work}/dip.mak" 2>"${logf}" \
    || die "run-dipcall failed, see ${logf}"
  make -j2 -f "${work}/dip.mak" >>"${logf}" 2>&1 || die "dipcall make failed, see ${logf}"
fi

# ---------------------------------------------------------------- 2. regions
grep -E '^chr([1-9]|1[0-9]|2[0-2])[[:space:]]' "${pre}.dip.bed" | sort -k1,1 -k2,2n >"${work}/dip.auto.bed"
grep -E '^chr([1-9]|1[0-9]|2[0-2])[[:space:]]' "${CONFIDENT_BED}" | sort -k1,1 -k2,2n >"${work}/bench.auto.bed"
bed_intersect "${work}/bench.auto.bed" "${work}/dip.auto.bed" >"${work}/region.bed"
[[ -s "${work}/region.bed" ]] || die "the v4.2.1 BED and dip.bed do not overlap"

# ------------------------------------------------------------------ 3. truth
# v4.2.1 het SNPs, in the benchmark BED and in the regions
"${BCFTOOLS}" view -T "${work}/bench.auto.bed" -m2 -M2 -v snps -g het "${TRUTH_SRC_VCF}" -Oz -o "${work}/v421.bench.vcf.gz"
"${BCFTOOLS}" view -T "${work}/region.bed" -m2 -M2 -v snps -g het "${TRUTH_SRC_VCF}" \
  | "${BCFTOOLS}" annotate -x 'INFO,^FORMAT/GT' -Oz -o "${work}/v421.region.vcf.gz"
# the assemblies' het SNPs there, phased hap1|hap2, alleles in upper case as
# in dipcall_run.sh
"${BCFTOOLS}" view -T "${work}/region.bed" -m2 -M2 -v snps -g het "${pre}.dip.vcf.gz" \
  | "${BCFTOOLS}" annotate -x 'INFO,^FORMAT/GT' \
  | awk -F'\t' -v OFS='\t' '/^#/ { print; next } { $4 = toupper($4); $5 = toupper($5); print }' \
  | "${BGZIP}" -c >"${work}/asm.region.vcf.gz"
for f in v421.bench v421.region asm.region; do "${TABIX}" -f -p vcf "${work}/${f}.vcf.gz"; done

# the same het SNP in both, the record and phase of the assemblies, under
# the sample's name (dipcall writes its own), as is the v4.2.1 copy compared
printf '%s\n' "${SAMPLE}" >"${work}/sample.txt"
"${BCFTOOLS}" isec -c none -n=2 -w2 -Oz -o "${work}/kept.vcf.gz" \
  "${work}/v421.region.vcf.gz" "${work}/asm.region.vcf.gz"
"${BCFTOOLS}" reheader -s "${work}/sample.txt" -o "${base}.tmp.vcf.gz" "${work}/kept.vcf.gz"
"${TABIX}" -f -p vcf "${base}.tmp.vcf.gz"
"${BCFTOOLS}" reheader -s "${work}/sample.txt" -o "${work}/v421.cmp.vcf.gz" "${work}/v421.region.vcf.gz"
"${TABIX}" -f -p vcf "${work}/v421.cmp.vcf.gz"
unphased="$("${BCFTOOLS}" query -f '[%GT]\n' "${base}.tmp.vcf.gz" | awk '!/\|/ { n++ } END { print n + 0 }')"
[[ "${unphased}" -eq 0 ]] || die "${unphased} kept records are not phased"

# ------------------------------------------------------------------ 4. check
n_bench="$("${BCFTOOLS}" index -n "${work}/v421.bench.vcf.gz")"
n_region="$("${BCFTOOLS}" index -n "${work}/v421.region.vcf.gz")"
n_asm="$("${BCFTOOLS}" index -n "${work}/asm.region.vcf.gz")"
n_kept="$("${BCFTOOLS}" index -n "${base}.tmp.vcf.gz")"
bp() { awk '{ s += $3 - $2 } END { printf "%.1f Mb", s / 1e6 }' "$1"; }
{
  printf '%s %s, phase from %s and %s\n' "${SAMPLE}" "${REFERENCE_NAME}" "$(basename "${ASM_PAT}")" "$(basename "${ASM_MAT}")"
  printf 'v4.2.1 benchmark BED, chr1-22:   %s, %d het SNPs\n' "$(bp "${work}/bench.auto.bed")" "${n_bench}"
  printf 'x dip.bed:                       %s, %d het SNPs (%.2f%%)\n' "$(bp "${work}/region.bed")" "${n_region}" \
    "$(awk -v a="${n_region}" -v b="${n_bench}" 'BEGIN { print b ? 100 * a / b : 0 }')"
  printf 'assemblies, het SNPs there:      %d\n' "${n_asm}"
  printf 'kept, phased by the assemblies:  %d (%.2f%% of the region, %.2f%% of the benchmark)\n' "${n_kept}" \
    "$(awk -v a="${n_kept}" -v b="${n_region}" 'BEGIN { print b ? 100 * a / b : 0 }')" \
    "$(awk -v a="${n_kept}" -v b="${n_bench}" 'BEGIN { print b ? 100 * a / b : 0 }')"
  printf 'v4.2.1 het SNPs left out (not het, or other alleles, in the assemblies): %d\n' "$((n_region - n_kept))"
} >"${check}.tmp"
# the new phase against v4.2.1's, where v4.2.1 has one
"${WHATSHAP}" compare --names asm,v4.2.1 --only-snvs --tsv-pairwise "${work}/compare.tsv" \
    "${base}.tmp.vcf.gz" "${work}/v421.cmp.vcf.gz" >"${work}/compare.txt" 2>&1 \
  || warn "whatshap compare failed, see ${work}/compare.txt"
if [[ -s "${work}/compare.tsv" ]]; then
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    { s += $h["all_switches"]; a += $h["all_assessed_pairs"]; f += $h["blockwise_hamming"]; c += $h["covered_variants"] }
    END { printf "against the v4.2.1 phase: %d switches in %d pairs (%.3f%%), Hamming %d of %d (%.3f%%)\n",
          s, a, a ? 100 * s / a : 0, f, c, c ? 100 * f / c : 0 }' "${work}/compare.tsv" >>"${check}.tmp"
fi
mv -f "${check}.tmp" "${check}"

cp "${work}/region.bed" "${base}.bed"
mv -f "${base}.tmp.vcf.gz.tbi" "${base}.vcf.gz.tbi"
mv -f "${base}.tmp.vcf.gz" "${base}.vcf.gz"
cat "${check}" >&2
log "${SAMPLE} ${REFERENCE_NAME}: ${base}.vcf.gz"
