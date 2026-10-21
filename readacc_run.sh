#!/usr/bin/env bash
# Read-level haplotag accuracy, with every phase block oriented by its variants.
#
# HP labels mean something only within a phase block: a phaser fixes the phase
# of variants within a block, but which haplotype it calls 1 is arbitrary and
# can differ from block to block. Each block of the scored VCF is oriented
# against the GIAB truth from its own variants - phased as in the truth or the
# other way round, majority wins - and the labels of its reads are flipped with
# it. Then every read is compared with its label from haplotagging against the
# truth. haplotag_accuracy in eval_run.sh instead flips each block to fit its
# reads, which can never score a block below 50% and always scores a one-read
# block as right; it also keys blocks by phaseset alone, so two blocks with the
# same start on different chromosomes merge. This script replaces it.
#
# No whatshap is run: the inputs are what eval_run.sh left on disk, and a run
# is scored only once its row.tsv is there.
#
# metrics/<name>.readacc.tsv, one line:
#   reference region caller platform depth phaser
#   reads_truth_tagged  reads with H1/H2 against the truth
#   reads_tagged        reads with H1/H2 against this VCF
#   reads_compared      tagged by both, in a block that shares a SNP with the truth
#   reads_correct       of those, on the truth's haplotype once the block is oriented
#   reads_unoriented    tagged by both, in a block that shares no SNP with the truth
#   read_precision      reads_correct / reads_compared
#   read_recall         reads_correct / reads_truth_tagged
#   blocks_oriented     blocks that share a SNP with the truth
#   blocks_flipped      of those, the ones phased the other way round
# metrics/summary_readacc.tsv - every line, with a header
#
#   bash readacc_run.sh                    the grid of config.sh, REFERENCE_NAME
#   READACC_FORCE=1 bash readacc_run.sh    recompute lines that already exist
#
# READACC_MIN_SNPS and READACC_AUTOSOMES (config.sh) choose the rules, and the
# file: with 1 and 0 the original <name>.readacc.tsv, otherwise
# <name>.readacc[_auto][_min<N>].tsv and summary_readacc[...].tsv beside it.
#   READACC_MIN_SNPS=N   blocks of fewer than N phased het SNPs orient nothing
#                        and their reads count as untagged
#   READACC_AUTOSOMES=1  reads and blocks on chr1-22 only
# blocks_oriented and blocks_flipped then count the blocks that pass.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

READACC_JOBS="${READACC_JOBS:-8}"
READACC_FORCE="${READACC_FORCE:-}"

mkdir -p "${TMP_DIR}"

READACC_SUFFIX="$(readacc_suffix)"
log "read accuracy rules: min ${READACC_MIN_SNPS} SNPs per block, autosomes only ${READACC_AUTOSOMES} -> *.${READACC_SUFFIX}.tsv"

# the truth genotypes, read once and shared by every run of this reference
TRUTH_GT="${TMP_DIR}/truth_gt_${REFERENCE_NAME}.tsv"
if [[ ! -s "${TRUTH_GT}" || "${TRUTH_PHASED}" -nt "${TRUTH_GT}" ]]; then
  log "${REFERENCE_NAME}: truth genotypes"
  "${BCFTOOLS}" query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "${TRUTH_PHASED}" >"${TRUTH_GT}.$$"
  mv -f "${TRUTH_GT}.$$" "${TRUTH_GT}"
fi

readacc_one() {
  set -euo pipefail

  local caller="$1" platform="$2" depth="$3" phaser="$4"
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
  local out="${METRICS_DIR}/${name}.${READACC_SUFFIX}.tsv"
  local row="${METRICS_DIR}/${name}.row.tsv"

  [[ -s "${out}" && -z "${READACC_FORCE}" ]] && return 0
  # scored only once eval_run.sh has finished this run
  [[ -s "${row}" && ! -d "${row}.lock" ]] || return 0

  local vcf="${PHASE_DIR}/${name}.vcf.gz"
  local tags="${METRICS_DIR}/${name}.haplotag.tsv"
  local truth_tags="${METRICS_DIR}/$(vcf_name "${platform}" "${depth}" truthtags).tsv"
  require_file "${vcf}" "phased VCF"
  require_file "${tags}" "haplotag list"
  require_file "${truth_tags}" "truth haplotag list"

  local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
  "${BCFTOOLS}" query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' "${vcf}" >"${tmp}/test_gt.tsv"
  block_orientation "${TRUTH_GT}" "${tmp}/test_gt.tsv" >"${tmp}/orient.tsv"

  if [[ "${READACC_SUFFIX}" != readacc ]]; then
    block_sizes "${tmp}/test_gt.tsv" >"${tmp}/blocks.tsv"
    awk -F'\t' -v OFS='\t' -v min="${READACC_MIN_SNPS}" -v auto="${READACC_AUTOSOMES}" '
      NR == FNR { size[$1 OFS $2] = $3; next }
      auto && $1 !~ /^chr([1-9]|1[0-9]|2[0-2])$/ { next }
      size[$1 OFS $2] + 0 >= min' "${tmp}/blocks.tsv" "${tmp}/orient.tsv" >"${tmp}/orient.kept.tsv"
    mv -f "${tmp}/orient.kept.tsv" "${tmp}/orient.tsv"
  fi

  local truth_tagged tagged compared correct unoriented
  if [[ "${READACC_SUFFIX}" == readacc ]]; then
    read -r truth_tagged tagged compared correct unoriented \
      < <(read_accuracy "${tmp}/orient.tsv" "${truth_tags}" "${tags}")
  else
    read -r truth_tagged tagged compared correct unoriented \
      < <(read_accuracy_filtered "${tmp}/orient.tsv" "${tmp}/blocks.tsv" "${truth_tags}" "${tags}" \
            "${READACC_MIN_SNPS}" "${READACC_AUTOSOMES}")
  fi
  [[ "${truth_tagged}" =~ ^[0-9]+$ && "${compared}" =~ ^[0-9]+$ ]] \
    || die "${name}: read accuracy came out empty"

  local blocks flipped
  blocks="$(wc -l <"${tmp}/orient.tsv" | tr -d ' ')"
  flipped="$(awk -F'\t' '$5 == 1' "${tmp}/orient.tsv" | wc -l | tr -d ' ')"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${REFERENCE_NAME}" "${REGION_NAME}" "${caller}" "${platform}" "${depth}" "${phaser}" \
      "${truth_tagged}" "${tagged}" "${compared}" "${correct}" "${unoriented}" \
      "$(awk -v a="${correct}" -v b="${compared}" 'BEGIN { printf "%.6f", (b ? a / b : 0) }')" \
      "$(awk -v a="${correct}" -v b="${truth_tagged}" 'BEGIN { printf "%.6f", (b ? a / b : 0) }')" \
      "${blocks}" "${flipped}" \
    >"${out}.$$"
  mv -f "${out}.$$" "${out}"
  rm -rf "${tmp}"

  log "${name}: precision $(cut -f12 "${out}"), recall $(cut -f13 "${out}")"
}

export -f readacc_one block_orientation read_accuracy vcf_name require_file log warn die \
          block_sizes read_accuracy_filtered
export METRICS_DIR PHASE_DIR TMP_DIR TRUTH_GT REFERENCE_NAME REGION_NAME SAMPLE \
       BCFTOOLS READACC_FORCE READACC_SUFFIX READACC_MIN_SNPS READACC_AUTOSOMES

log "read accuracy: ${#CALLERS[@]} callers x ${#PLATFORMS[@]} platforms x ${#DEPTHS[@]} depths x ${#PHASERS[@]} phasers on ${REFERENCE_NAME}"

joblog="${LOGS_DIR}/parallel_readacc_${RUN_TAG}.log"
keep_joblog "${joblog}"

"${PARALLEL}" \
    --joblog "${joblog}" \
    -j "${READACC_JOBS}" \
    readacc_one {1} {2} {3} {4} \
    ::: "${CALLERS[@]}" \
    ::: "${PLATFORMS[@]}" \
    ::: "${DEPTHS[@]}" \
    ::: "${PHASERS[@]}" \
  || true

failed="$(awk 'NR > 1 && $7 != 0 { n++ } END { print n + 0 }' "${joblog}")"
[[ "${failed}" -eq 0 ]] || warn "${failed} runs failed, see ${joblog}"

{
  printf 'reference\tregion\tcaller\tplatform\tdepth\tphaser'
  printf '\treads_truth_tagged\treads_tagged\treads_compared\treads_correct\treads_unoriented'
  printf '\tread_precision\tread_recall\tblocks_oriented\tblocks_flipped\n'
  cat "${METRICS_DIR}"/*."${READACC_SUFFIX}".tsv 2>/dev/null || true
} >"${METRICS_DIR}/summary_${READACC_SUFFIX}.tsv.$$"
mv -f "${METRICS_DIR}/summary_${READACC_SUFFIX}.tsv.$$" "${METRICS_DIR}/summary_${READACC_SUFFIX}.tsv"

log "read accuracy: ${METRICS_DIR}/summary_${READACC_SUFFIX}.tsv"
