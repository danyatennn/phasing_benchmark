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
    tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/mosdepth.XXXXXX")"
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


# source_bam <platform> [reference]
# reference defaults to REFERENCE_NAME
source_bam() {
    local platform="$1"
    local reference="${2:-${REFERENCE_NAME}}"

    if [[ "${reference}" != "GRCh38" ]]; then
        printf '%s/%s.bam\n' "${T2T_SRC_DIR}" "${platform}"
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
is_done() {
  [[ -s "$1" && -s "$1.tbi" ]]
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

# stats_col <whatshap stats tsv> <column>  - the ALL row, or the only row
stats_col() {
  awk -F'\t' -v w="$2" '
    NR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h); if (h == w) c = i } next }
    NR == 2 { v = $c }
    $2 == "ALL" { v = $c }
    END { print v }' "$1"
}

# cmp_sum <whatshap compare tsv> <column>  - summed over chromosomes, because
# compare prints one row per chromosome and no aggregate
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

# extractHAIRS read-type flag
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
