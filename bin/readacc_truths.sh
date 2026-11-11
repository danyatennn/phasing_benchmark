#!/usr/bin/env bash
# Read-level accuracy against each read truth, every block oriented by its SNPs
# as in readacc_run.sh (lib.sh: block_orientation, the same rules: chr1-22,
# blocks of READACC_MIN_SNPS phased het SNPs or more). Only the label each read
# is compared with changes:
#   truthtags  whatshap haplotag against the phased truth (what readacc_run.sh uses)
#   matpat     bin/matpat_truth.sh: alignment to the two parental copies of Q100
#   yak        bin/yak_truth.sh: parental k-mers
# So a difference between truths is a difference of the labels, not of how the
# blocks were oriented (readacc_anchor.sh orients each block by its longest
# labelled read, and the parental truths err most on the longest reads).
#
# A parental label is turned into the truth's haplotypes once per cell: on each
# chromosome, PAT is H1 or H2 as most reads labelled by both PAT/MAT and the
# truth VCF say (the truth VCF keeps one parent on one side of a chromosome).
# These tables are reads/<cell>.<truth>.htags.tsv, whatshap haplotag-list layout.
#
# metrics/<run>.readacc_<truth>_<rules>.tsv, one line:
#   reference region caller platform depth phaser truth_source
#   reads_truth_tagged reads_tagged reads_compared reads_correct reads_unoriented
#   read_precision read_recall          as readacc_run.sh, against this truth
#   reads_truth_common                  reads labelled by this truth and the truth VCF
#   reads_compared_common reads_correct_common
#   read_precision_common read_recall_common   the same over those reads only:
#                                       the numbers to compare between truths
# metrics/summary_readacc_<truth>_<rules>.tsv - every line, with a header
#
#   bash bin/readacc_truths.sh                                 30x, HG002, GRCh38
#   REFERENCE_NAME=T2T bash bin/readacc_truths.sh
#   PHASING_CALLERS="truthindel truthall" bash bin/readacc_truths.sh
#   READ_TRUTHS=truthtags SAMPLE=HG005 bash bin/readacc_truths.sh   other samples: truthtags only
#
# A line newer than the run's VCF, its haplotag list and the labels is kept.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

READ_TRUTHS=(${READ_TRUTHS:-truthtags matpat yak})
TRUTH_DEPTHS=(${PHASING_DEPTHS:-30})
TRUTHS_JOBS="${TRUTHS_JOBS:-6}"
RULES="$(readacc_suffix)"
RULES="${RULES#readacc}"            # _auto_min2
mkdir -p "${TMP_DIR}" "${LOGS_DIR}/readacc"

# the truth genotypes, as readacc_run.sh reads them
TRUTH_GT="${TMP_DIR}/truth_gt_${REFERENCE_NAME}.tsv"
if [[ ! -s "${TRUTH_GT}" || "${TRUTH_PHASED}" -nt "${TRUTH_GT}" ]]; then
  log "${REFERENCE_NAME}: truth genotypes"
  "${BCFTOOLS}" query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "${TRUTH_PHASED}" >"${TRUTH_GT}.$$"
  mv -f "${TRUTH_GT}.$$" "${TRUTH_GT}"
fi

# labels_of <cell> <truth> - the truth's labels as a haplotag list
labels_of() {
  if [[ "$2" == truthtags ]]; then
    echo "${METRICS_DIR}/$1_truthtags.tsv"
  else
    echo "${READS_DIR}/$1.$2.htags.tsv"
  fi
}

# make_labels <cell> <truth> - PAT/MAT -> H1/H2 of the truth, per chromosome
make_labels() {
  local cell="$1" truth="$2"
  local src="${READS_DIR}/${cell}.${truth}.tsv.gz"
  local table="${READS_DIR}/${cell}.reads.tsv.gz"
  local tags="${METRICS_DIR}/${cell}_truthtags.tsv"
  local out; out="$(labels_of "${cell}" "${truth}")"
  [[ -s "${src}" && -s "${table}" && -s "${tags}" ]] || { warn "${cell}: no ${truth} labels, read table or truth tags"; return 0; }
  newer_than "${out}" "${src}" && newer_than "${out}" "${tags}" && return 0
  awk -F'\t' -v OFS='\t' -v log_out="${LOGS_DIR}/readacc/${cell}.${truth}.htags.log" '
    FILENAME == ARGV[1] { if (FNR > 1) chrom[$1] = $2; next }
    FILENAME == ARGV[2] { if (!($1 in tt) && ($2 == "H1" || $2 == "H2")) tt[$1] = $2; next }
    ($2 == "PAT" || $2 == "MAT") && ($1 in chrom) && !($1 in lab) {
      lab[$1] = $2; c = chrom[$1]
      if ($1 in tt) { n[c]++; if (($2 == "PAT") == (tt[$1] == "H1")) same[c]++ }
    }
    END {
      for (c in n) {
        pat_h1[c] = (same[c] >= n[c] - same[c])
        printf("%s\t%s\t%d reads labelled by both, %.2f%% agree\n", c, pat_h1[c] ? "PAT=H1" : "PAT=H2",
               n[c], 100 * (pat_h1[c] ? same[c] : n[c] - same[c]) / n[c]) > log_out
      }
      print "#readname", "haplotype", "phaseset", "chromosome"
      for (r in lab) {
        c = chrom[r]
        if (!(c in n)) continue
        print r, ((lab[r] == "PAT") == pat_h1[c]) ? "H1" : "H2", ".", c
      }
    }' <(zcat "${table}") "${tags}" <(zcat "${src}") >"${out}.$$"
  mv -f "${out}.$$" "${out}"
  log "${cell}: ${truth} labels -> $(basename "${out}")"
}

# accuracy <orient> <blocks> <labels> <truth tags> <test tags> <min> <autosomes>
# read_accuracy_filtered of lib.sh against <labels>, plus the same counts over the
# reads that also carry a label of <truth tags>
accuracy() {
  awk -F'\t' -v min="$6" -v auto="$7" '
    function autosomal(c) { return !auto || c ~ /^chr([1-9]|1[0-9]|2[0-2])$/ }
    FILENAME == ARGV[1] { flip[$1 FS $2] = $5; next }
    FILENAME == ARGV[2] { size[$1 FS $2] = $3; next }
    /^#/ { next }
    FILENAME == ARGV[3] {
      if ($1 in in_truth) next
      in_truth[$1] = 1
      if (!autosomal($4)) next
      if ($2 == "H1" || $2 == "H2") { truth[$1] = $2; truth_tagged++ }
      next
    }
    FILENAME == ARGV[4] {
      if ($1 in in_vcf) next
      in_vcf[$1] = 1
      if (autosomal($4) && ($2 == "H1" || $2 == "H2") && ($1 in truth)) { common[$1] = 1; truth_common++ }
      next
    }
    {
      if ($1 in in_test) next
      in_test[$1] = 1
      if ($2 != "H1" && $2 != "H2") next
      if (!autosomal($4)) next
      b = $4 FS $3
      if (size[b] + 0 < min) next
      tagged++
      if (!($1 in truth)) next
      if (!(b in flip)) { unoriented++; next }
      label = $2
      if (flip[b] == 1) label = (label == "H1") ? "H2" : "H1"
      compared++
      ok = (label == truth[$1])
      correct += ok
      if ($1 in common) { compared_common++; correct_common += ok }
    }
    END { print truth_tagged + 0, tagged + 0, compared + 0, correct + 0, unoriented + 0,
                truth_common + 0, compared_common + 0, correct_common + 0 }
  ' "$1" "$2" "$3" "$4" "$5"
}

truth_one() {
  set -euo pipefail
  local truth="$1" caller="$2" platform="$3" depth="$4" phaser="$5"
  phaser_applies "${phaser}" "${platform}" || return 0
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
  local row="${METRICS_DIR}/${name}.row.tsv"
  local vcf="${PHASE_DIR}/${name}.vcf.gz"
  local tags="${METRICS_DIR}/${name}.haplotag.tsv"
  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local labels; labels="$(labels_of "${cell}" "${truth}")"
  local truth_tags; truth_tags="$(labels_of "${cell}" truthtags)"
  local out="${METRICS_DIR}/${name}.readacc_${truth}${RULES}.tsv"
  [[ -s "${row}" && ! -d "${row}.lock" && -s "${vcf}" && -s "${tags}" && -s "${labels}" && -s "${truth_tags}" ]] || return 0
  newer_than "${out}" "${vcf}" && newer_than "${out}" "${tags}" && newer_than "${out}" "${labels}" && return 0

  local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
  "${BCFTOOLS}" query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' "${vcf}" >"${tmp}/test_gt.tsv"
  block_orientation "${TRUTH_GT}" "${tmp}/test_gt.tsv" >"${tmp}/orient.tsv"
  block_sizes "${tmp}/test_gt.tsv" >"${tmp}/blocks.tsv"
  local n
  # the truth tags through a pipe: for truthtags they are also the labels, and
  # awk tells its inputs apart by name
  read -r -a n < <(accuracy "${tmp}/orient.tsv" "${tmp}/blocks.tsv" "${labels}" <(cat "${truth_tags}") "${tags}" \
                            "${READACC_MIN_SNPS}" "${READACC_AUTOSOMES}")
  rm -rf "${tmp}"
  [[ ${#n[@]} -eq 8 ]] || die "${name}: ${truth} accuracy came out empty"
  r() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.6f", (b ? a / b : 0) }'; }
  printf '%s\t' "${REFERENCE_NAME}" "${REGION_NAME}" "${caller}" "${platform}" "${depth}" "${phaser}" "${truth}" \
      "${n[0]}" "${n[1]}" "${n[2]}" "${n[3]}" "${n[4]}" "$(r "${n[3]}" "${n[2]}")" "$(r "${n[3]}" "${n[0]}")" \
      "${n[5]}" "${n[6]}" "${n[7]}" "$(r "${n[7]}" "${n[6]}")" >"${out}.$$"
  printf '%s\n' "$(r "${n[7]}" "${n[5]}")" >>"${out}.$$"
  mv -f "${out}.$$" "${out}"
  log "${name} ${truth}: precision $(cut -f13 "${out}"), recall $(cut -f14 "${out}"), recall on common reads $(cut -f19 "${out}")"
}

export -f truth_one accuracy labels_of block_orientation block_sizes vcf_name bam_key phaser_applies \
          newer_than log warn die
export METRICS_DIR PHASE_DIR READS_DIR LOGS_DIR TMP_DIR TRUTH_GT REFERENCE_NAME REGION_NAME SAMPLE \
       BCFTOOLS RULES READACC_MIN_SNPS READACC_AUTOSOMES

log "${SAMPLE} ${REFERENCE_NAME}: read accuracy against ${READ_TRUTHS[*]}, blocks oriented by SNPs (${RULES#_}), ${CALLERS[*]} x ${PLATFORMS[*]} x ${TRUTH_DEPTHS[*]}x x ${PHASERS[*]}"
for platform in "${PLATFORMS[@]}"; do
  for depth in "${TRUTH_DEPTHS[@]}"; do
    cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
    for truth in "${READ_TRUTHS[@]}"; do
      [[ "${truth}" == truthtags ]] || make_labels "${cell}" "${truth}"
    done
  done
done

"${PARALLEL}" -j "${TRUTHS_JOBS}" truth_one {1} {2} {3} {4} {5} \
    ::: "${READ_TRUTHS[@]}" ::: "${CALLERS[@]}" ::: "${PLATFORMS[@]}" ::: "${TRUTH_DEPTHS[@]}" ::: "${PHASERS[@]}" || true

for truth in "${READ_TRUTHS[@]}"; do
  summary="${METRICS_DIR}/summary_readacc_${truth}${RULES}.tsv"
  {
    printf 'reference\tregion\tcaller\tplatform\tdepth\tphaser\ttruth_source'
    printf '\treads_truth_tagged\treads_tagged\treads_compared\treads_correct\treads_unoriented'
    printf '\tread_precision\tread_recall\treads_truth_common\treads_compared_common'
    printf '\treads_correct_common\tread_precision_common\tread_recall_common\n'
    files=("${METRICS_DIR}"/*.readacc_${truth}${RULES}.tsv)
    [[ -e "${files[0]}" ]] && cat "${files[@]}"
    true
  } >"${summary}.$$"
  mv -f "${summary}.$$" "${summary}"
  log "read accuracy against ${truth}: ${summary}"
done
