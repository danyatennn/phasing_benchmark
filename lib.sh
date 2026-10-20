#!/usr/bin/env bash
# Shared helpers

log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date '+%H:%M:%S')" "$*" >&2; exit 1; }
warn() { printf '[%s] WARN: %s\n'  "$(date '+%H:%M:%S')" "$*" >&2; }

measure_depth() {
    local bam="$1"
    local tmpdir
    local prefix
    local mean_depth
    tmpdir="$(mktemp -d "${TMP_DIR:-${TMPDIR:-/tmp}}/mosdepth.XXXXXX")"
    # shellcheck disable=SC2064
    trap "rm -rf '${tmpdir}'" RETURN
    prefix="${tmpdir}/depth"
    "${MOSDEPTH}" \
        -n \
        -t "${MOSDEPTH_THREADS:-${THREADS_PER_JOB:-1}}" \
        "${prefix}" \
        "${bam}"
    mean_depth="$(
        awk '$1 == "total" {print $4}' \
            "${prefix}.mosdepth.summary.txt"
    )"
    [[ -n "${mean_depth}" ]] \
        || die "mosdepth did not report total mean depth for ${bam}"
    printf '%s\n' "${mean_depth}"
}

# for_each_bam <fn>  - calls fn <platform> <depth> <seed>
for_each_bam() {
  local fn="$1" platform depth seed
  for platform in "${PLATFORMS[@]}"; do
    for depth in "${DEPTHS[@]}"; do
        "${fn}" "${platform}" "${depth}"
    done
  done
}

# for_each_run <fn> - calls fn <exp> <platform> <depth> <seed> <phaser>
for_each_run() {
  local fn="$1" exp platform depth seed phaser
  for exp in "${EXPERIMENTS[@]}"; do
    for platform in "${PLATFORMS[@]}"; do
      for depth in "${DEPTHS[@]}"; do
        for phaser in "${PHASERS[@]}"; do
          "${fn}" "${exp}" "${platform}" "${depth}" "${phaser}"
        done
      done
    done
  done
}

require_file() { [[ -f "$1" ]] || die "missing $2: $1"; }


# sample_path <sample> <item> - the value of <item> for <sample> in samples.tsv
# (sample, item, value; # starts a comment); empty when it is not there. A value
# written as $NAME is the variable NAME of data_config_germline.sh, so the paths
# are kept in one place; it is empty if that variable is not set.
sample_path() {
  [[ -s "${SAMPLES_TSV:-}" ]] || return 0
  local value name
  value="$(awk -F'\t' -v s="$1" -v i="$2" '
    /^#/ || NF < 3 { next }
    $1 == s && $2 == i { print $3; exit }' "${SAMPLES_TSV}")"
  if [[ "${value}" =~ ^\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?$ ]]; then
    name="${BASH_REMATCH[1]}"
    value="${!name:-}"
  fi
  [[ -n "${value}" ]] && printf '%s\n' "${value}"
  return 0
}

# sample_platforms <sample> - the platforms with a source BAM in samples.tsv,
# in the order the grid uses
sample_platforms() {
  local p out=()
  for p in ont4k ont5k hifi revio; do
    [[ -n "$(sample_path "$1" "bam_${p}")" ]] && out+=("${p}")
  done
  printf '%s\n' "${out[*]}"
}

# source_bam <platform> [reference]
# reference defaults to REFERENCE_NAME
source_bam() {
    local platform="$1"
    local reference="${2:-${REFERENCE_NAME}}"

    if [[ "${reference}" != "GRCh38" ]]; then
        printf '%s/%s.bam\n' "${T2T_SRC_DIR}" "${platform}"
        return 0
    fi

    # every sample but HG002 lists its BAMs in samples.tsv; config.sh resolves
    # them into SAMPLE_BAM_<platform>, exported, so a function run by GNU
    # parallel sees them without data_config_germline.sh
    if [[ "${SAMPLE}" != HG002 ]]; then
        local path var="SAMPLE_BAM_${platform}"
        path="${!var:-}"
        [[ -n "${path}" ]] || path="$(sample_path "${SAMPLE}" "bam_${platform}")"
        [[ -n "${path}" ]] || {
            printf 'no bam_%s for %s in %s\n' "${platform}" "${SAMPLE}" "${SAMPLES_TSV}" >&2
            return 1
        }
        printf '%s\n' "${path}"
        return 0
    fi

    case "${platform}" in
        ont5k) printf '%s\n' "${ONT_HG002_V14_5KHZ_NO_ALT_BAM}";;
        ont4k) printf '%s\n' "${ONT_HG002_dorado_NO_ALT_BAM}";;
        hifi)  printf '%s\n' "${PB_HG002_20K_NO_ALT_BAM}";;
        revio) printf '%s\n' "${HG002_REVIO_BAM}";;
        *)
            printf 'unknown platform: %s\n' "${platform}" >&2
            return 1
            ;;
    esac
}

source_depth() {
  awk -v p="$1" '$1==p {print $2}' "${DEPTH_TSV}"
}

# require_vars <name>...  - fail with a clear message if any is unset or empty
require_vars() {
  local missing=() v
  for v in "$@"; do
    [[ -n "${!v:-}" ]] || missing+=("${v}")
  done
  [[ ${#missing[@]} -eq 0 ]] \
    || die "unset/empty config variable(s): ${missing[*]} (check data_config_germline.sh)"
}

# require_cmds <cmd>...  - fail with a clear message if any is not executable
require_cmds() {
  local missing=() c
  for c in "$@"; do
    command -v "${c}" >/dev/null 2>&1 || missing+=("${c}")
  done
  [[ ${#missing[@]} -eq 0 ]] || die "command(s) not found: ${missing[*]}"
}
# bam_key <platform> <ref> <depth> - the downsample directory name
bam_key() {
  printf '%s_%s_%s_%sx_wgs\n' "$1" "$2" "${SAMPLE}" "$3"
}

# sub_bam <platform> <ref> <depth> - the downsampled BAM for that cell
sub_bam() {
  printf '%s/%s/input.bam\n' "${DATA_DIR}" "$(bam_key "$1" "$2" "$3")"
}

# clair3_platform <platform> - Clair3's --platform value
clair3_platform() {
  case "$1" in
    ont4k|ont5k) printf 'ont\n' ;;
    hifi|revio)  printf 'hifi\n' ;;
    *) die "no Clair3 platform mapping for '$1'" ;;
  esac
}

# clair3_model <platform> - model directory name inside /opt/models
clair3_model() {
  local var="CLAIR3_MODEL_$1"
  local model="${!var:-}"
  [[ -n "${model}" ]] || die "no Clair3 model configured for '$1' (set ${var})"
  printf '%s\n' "${model}"
}

# is_done <vcf.gz> - non-empty VCF with a non-empty index
keep_joblog() {
    [[ -s "$1" ]] && cp -f "$1" "${1%.log}_$(date +%F_%H%M).log"
    return 0
}

stratum_bed()   { printf '%s/%s_%s.bed\n' "${PREPARED_DIR}" "${REFERENCE_NAME}" "$1"; }
stratum_truth() { printf '%s/%s_%s.phased.vcf.gz\n' "${STRATA_TRUTH_DIR}" "${REFERENCE_NAME}" "$1"; }

# strata.tsv is tier<TAB>category<TAB>stratum, written by bin/prepare_strata.sh.
# STRATA_ONLY and STRATA_TIER narrow it for splitting work across hosts.
strata_names() {
    awk -F'\t' -v tier="${STRATA_TIER}" -v only="${STRATA_ONLY}" '
        NR == 1 { next }
        tier != "" && $1 != tier { next }
        { if (only == "") { print $3; next }
          n = split(only, want, /[ ,]+/)
          for (i = 1; i <= n; i++) if (want[i] == $3) print $3 }
    ' "${STRATA_LIST}"
}

# both files must be sorted by chrom then start and non-overlapping within a
# file, which is how GIAB ships them
bed_intersect() {
    awk -v OFS='\t' '
        NR == FNR { n[$1]++; s[$1, n[$1]] = $2; e[$1, n[$1]] = $3; next }
        {
          c = $1; a = $2; b = $3
          if (!(c in n)) next
          i = (c in ptr) ? ptr[c] : 1
          while (i <= n[c] && e[c, i] <= a) i++
          ptr[c] = i
          for (j = i; j <= n[c] && s[c, j] < b; j++) {
            lo = (s[c, j] > a) ? s[c, j] : a
            hi = (e[c, j] < b) ? e[c, j] : b
            if (hi > lo) print c, lo, hi
          }
        }' "$1" "$2"
}

bed_bp() { awk '{ s += $3 - $2 } END { print s + 0 }' "$1"; }

is_done() {
  [[ -s "$1" && -s "$1.tbi" ]]
}

# newer_than <file> <reference> - the file exists, is not empty, and was
# written after the reference: a result derived from it is still current
newer_than() {
  [[ -s "$1" && "$1" -nt "$2" ]]
}

# the columns of metrics/<run>.row.tsv (eval_run.sh, metrics.md)
ROW_FIELDS=32

# row_is_valid <row.tsv> - one line, ROW_FIELDS fields, the six labels set and every other field a number 
row_is_valid() {
  awk -F'\t' -v n="${ROW_FIELDS:-32}" '
    NR == 1 {
      if (NF != n) bad = 1
      for (i = 1; i <= 6 && !bad; i++) if ($i == "") bad = 1
      for (i = 7; i <= NF && !bad; i++)
        if (tolower($i) !~ /^[-+]?([0-9]+(\.[0-9]*)?([e][-+]?[0-9]+)?|nan|inf)$/) bad = 1
    }
    END { exit (bad || NR != 1) }' "$1"
}

# vcf_name <platform> <depth> <caller> [phaser]
# platform_reference_sample_depth_region_caller[_phaser]
vcf_name() {
  local n="$1_${REFERENCE_NAME}_${SAMPLE}_$2x_${REGION_NAME}_$3"
  [[ -n "${4:-}" ]] && n="${n}_$4"
  printf '%s\n' "${n}"
}

longphase_flag() {
  case "$1" in
    ont4k|ont5k) printf '%s\n' --ont ;;
    hifi|revio)  printf '%s\n' --pb ;;
  esac
}

# stats_col <whatshap stats tsv> <column> 
stats_col() {
  awk -F'\t' -v w="$2" '
    NR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h); if (h == w) c = i } next }
    NR == 2 { v = $c }
    $2 == "ALL" { v = $c }
    END { print v }' "$1"
}

# cmp_sum <whatshap compare tsv> <column>  - summed over chromosomes
cmp_sum() {
  awk -F'\t' -v w="$2" '
    NR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h); if (h == w) c = i } next }
    { s += $c }
    END { printf "%d\n", s }' "$1"
}

deepvariant_model() {
  case "$1" in
    ont4k|ont5k) printf 'ONT_R104\n' ;;
    hifi|revio)  printf 'PACBIO\n' ;;
  esac
}

longcalld_preset() {
  case "$1" in
    ont4k|ont5k) printf -- '--ont\n' ;;
    hifi|revio)  printf -- '--hifi\n' ;;
  esac
}

margin_params() {
  case "$1" in
    ont4k|ont5k) printf '%s/allParams.phase_vcf.ont.json\n' "${MARGIN_PARAMS_DIR}" ;;
    hifi|revio)  printf '%s/allParams.phase_vcf.pb-hifi.json\n' "${MARGIN_PARAMS_DIR}" ;;
  esac
}

# extract read-type flag
hapcut2_flag() {
  case "$1" in
    ont4k|ont5k) printf -- '--ont\n' ;;
    hifi|revio)  printf -- '--pacbio\n' ;;
  esac
}

# how many cells of this caller to run at once
call_jobs() {
  case "$1" in
    clair3)      printf '%s\n' "${CLAIR3_JOBS}" ;;
    deepvariant) printf '%s\n' "${DEEPVARIANT_JOBS}" ;;
    longcalld)   printf '%s\n' "${LONGCALLD_JOBS}" ;;
  esac
}

# haplotag_accuracy <truth tags> <test tags>
# Prints: correct <tab> compared <tab> accuracy
haplotag_accuracy() {
  awk -F'\t' '
    NR == FNR { if (FNR > 1 && $2 != "none") truth[$1] = $2; next }
    FNR > 1 && $2 != "none" && ($1 in truth) {
      block[$3] = 1
      if (truth[$1] == $2) agree[$3]++; else disagree[$3]++
    }
    END {
      for (b in block) {
        a = agree[b] + 0; d = disagree[b] + 0
        total += a + d
        correct += (a > d ? a : d)
      }
      printf "%d\t%d\t%.6f", correct + 0, total + 0, (total ? correct / total : 0)
    }' "$1" "$2"
}

# block_orientation <truth_gt.tsv> <test_gt.tsv>
#   truth_gt  chrom pos ref alt GT        - the phased truth
#   test_gt   chrom pos ref alt GT PS     - the VCF being scored
# Prints chrom, PS, same, opposite, flip for every phase block of the test VCF
# that shares at least one het SNP with the truth: how many of those SNPs are
# phased as in the truth and how many the other way round, and flip=1 when the
# block is the other way round on most of them.
block_orientation() {
  awk -F'\t' -v OFS='\t' '
    function swapped(gt,   a) { split(gt, a, "|"); return a[2] "|" a[1] }
    NR == FNR { truth[$1 FS $2 FS $3 FS $4] = $5; next }
    $6 == "." || $5 !~ /\|/ { next }
    {
      k = $1 FS $2 FS $3 FS $4
      if (!(k in truth)) next
      b = $1 OFS $6
      if ($5 == truth[k])               { same[b]++; seen[b] = 1 }
      else if ($5 == swapped(truth[k])) { opp[b]++;  seen[b] = 1 }
    }
    END {
      for (b in seen)
        print b, same[b] + 0, opp[b] + 0, (opp[b] > same[b]) ? 1 : 0
    }' "$1" "$2"
}

# read_accuracy <orientation.tsv> <truth tags> <test tags>
# The tags are whatshap haplotag lists: read, haplotype, phaseset, chromosome.
# A read takes the orientation of its block (chromosome + phaseset) and is
# compared with its truth label. A read listed twice counts once.
# Prints: truth_tagged tagged compared correct unoriented
#   truth_tagged  reads with H1/H2 against the truth
#   tagged        reads with H1/H2 against the test VCF
#   compared      tagged by both, in a block that shares a SNP with the truth
#   correct       of those, the truth's haplotype after orienting the block
#   unoriented    tagged by both, but their block shares no SNP with the truth
read_accuracy() {
  awk -F'\t' '
    FILENAME == ARGV[1] { flip[$1 FS $2] = $5; next }
    /^#/ { next }
    FILENAME == ARGV[2] {
      if ($1 in in_truth) next
      in_truth[$1] = 1
      if ($2 == "H1" || $2 == "H2") { truth[$1] = $2; truth_tagged++ }
      next
    }
    {
      if ($1 in in_test) next
      in_test[$1] = 1
      if ($2 != "H1" && $2 != "H2") next
      tagged++
      if (!($1 in truth)) next
      b = $4 FS $3
      if (!(b in flip)) { unoriented++; next }
      label = $2
      if (flip[b] == 1) label = (label == "H1") ? "H2" : "H1"
      compared++
      if (label == truth[$1]) correct++
    }
    END { print truth_tagged + 0, tagged + 0, compared + 0, correct + 0, unoriented + 0 }
  ' "$1" "$2" "$3"
}

# readacc_suffix - the file name readacc_run.sh writes for the current
# READACC_MIN_SNPS and READACC_AUTOSOMES: "readacc" for the original 1 and 0,
# otherwise readacc[_auto][_min<N>], so no setting overwrites another
readacc_suffix() {
  local min="${READACC_MIN_SNPS:-1}" auto="${READACC_AUTOSOMES:-0}" s="readacc"
  [[ "${auto}" == 1 ]] && s="${s}_auto"
  [[ "${min}" -gt 1 ]] && s="${s}_min${min}"
  printf '%s\n' "${s}"
}

# block_sizes <test_gt.tsv> - chrom, PS, number of phased het SNPs, for every
# phase block of the VCF being scored (test_gt as for block_orientation)
block_sizes() {
  awk -F'\t' -v OFS='\t' '
    $6 == "." || $5 !~ /\|/ { next }
    { split($5, a, "|"); if (a[1] == a[2]) next; n[$1 OFS $6]++ }
    END { for (b in n) print b, n[b] }' "$1"
}

# read_accuracy_filtered <orientation.tsv> <block sizes> <truth tags> <test tags> <min snps> <autosomes>
# read_accuracy with two rules on top:
#   a read tagged in a block of fewer than <min snps> phased het SNPs counts as
#   untagged - its label rests on single alleles, not on phase
#   with <autosomes> 1, reads and blocks off chr1-22 are left out entirely
# Prints the same five numbers as read_accuracy.
read_accuracy_filtered() {
  awk -F'\t' -v min="$5" -v auto="$6" '
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
      if (label == truth[$1]) correct++
    }
    END { print truth_tagged + 0, tagged + 0, compared + 0, correct + 0, unoriented + 0 }
  ' "$1" "$2" "$3" "$4"
}

# asm_dir <platform> <depth> <assembler> - where an assembly lives. It is built
# from the reads of the ASM_READS_REF cell whichever reference it is scored on.
asm_dir() {
  printf '%s/%s_%s\n' "${ASM_DIR}" "$(bam_key "$1" "${ASM_READS_REF}" "$2")" "$3"
}

# bam_to_fastq <bam> <out.fq.gz> - every read once: samtools fastq drops the
# secondary and supplementary records
bam_to_fastq() {
  "${SAMTOOLS}" view -u -@ "${SORT_THREADS}" "$1" ${ASM_REGION} \
    | "${SAMTOOLS}" fastq -@ "${SORT_THREADS}" -n - \
    | "${BGZIP}" -@ "${SORT_THREADS}" -c >"$2.tmp"
  mv -f "$2.tmp" "$2"
}

# assembly_phase_sets <hap1.bed> <hap2.bed> <snps.vcf> [chromosome]
#   hapN.bed  chrom start end contig - where each contig of haplotype N aligns
#   snps.vcf  the assembly's het SNPs, sorted
# Prints chrom, pos, phase set, the set named, as whatshap does, after the
# position of its first SNP. SNPs not on a contig in both haplotypes are dropped.
#   default     an assembly from reads alone vouches for the phase of two SNPs
#               only while both sit on the same hap1 contig and the same hap2
#               contig, so a set starts wherever either contig changes
#   chromosome  a trio assembly keeps hap1 paternal from contig to contig, so
#               the set is the whole chromosome
assembly_phase_sets() {
  awk -F'\t' -v OFS='\t' -v whole="${4:-}" '
    # the contig of haplotype h over chrom:pos; SNPs come in increasing
    # position, so the search resumes where the previous one stopped
    function contig(h, chrom, pos,   j) {
      if (!((h, chrom) in n)) return ""
      if (!((h, chrom) in at)) at[h, chrom] = 1
      while (at[h, chrom] <= n[h, chrom] && e[h, chrom, at[h, chrom]] < pos) at[h, chrom]++
      for (j = at[h, chrom]; j <= n[h, chrom] && s[h, chrom, j] < pos; j++)
        if (e[h, chrom, j] >= pos) return c[h, chrom, j]
      return ""
    }
    FILENAME == ARGV[1] || FILENAME == ARGV[2] {
      h = (FILENAME == ARGV[1]) ? 1 : 2
      k = ++n[h, $1]; s[h, $1, k] = $2; e[h, $1, k] = $3; c[h, $1, k] = $4
      next
    }
    /^#/ { next }
    {
      a = contig(1, $1, $2); b = contig(2, $1, $2)
      if (a == "" || b == "") next
      if ($1 != chrom || (whole == "" && (a != last_a || b != last_b))) ps = $2
      print $1, $2, ps
      chrom = $1; last_a = a; last_b = b
    }' "$1" "$2" "$3"
}

# add_phase_sets <phase sets tsv> <snps.vcf> <sample> - the VCF with GT:PS,
# only the SNPs that got a phase set, and the sample renamed
add_phase_sets() {
  awk -F'\t' -v OFS='\t' -v sample="$3" '
    NR == FNR { ps[$1, $2] = $3; next }
    /^##/ { print; next }
    /^#CHROM/ {
      print "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set: one pair of hap1/hap2 contigs\">"
      $10 = sample
      print
      next
    }
    ($1, $2) in ps { $9 = "GT:PS"; $10 = $10 ":" ps[$1, $2]; print }
  ' "$1" "$2"
}
