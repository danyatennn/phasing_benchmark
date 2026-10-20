#!/usr/bin/env bash
# Parental k-mer databases for the trio assemblies, built once per child from
# the parents' reads (PAT_READS, MAT_READS: config.sh, by default pat_reads and
# mat_reads of the sample in samples.tsv):
#
#   <TRIO_DIR>/pat.yak, mat.yak                         hifiasm -1 / -2
#   <TRIO_DIR>/pat.k30.only.meryl, mat.k30.only.meryl   verkko --hap-kmers ... trio
#
# TRIO_DIR is data/trio of the sample (samples/<sample>/data/trio for any but
# HG002). Only reads are needed, no VCF and no parental assembly. They can be
# Illumina R1/R2 FASTQ or accurate long reads; a BAM or CRAM is turned into
# FASTQ once (primary records only, every read once) and kept in TRIO_DIR.
#
#   SAMPLE=HG005 bash trio_kmers.sh
#   TRIO_SUBSAMPLE=0.5 bash trio_kmers.sh   keep half of each parental BAM
#                                           (61x Revio -> ~30x, plenty for k-mers)

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

TRIO_SUBSAMPLE="${TRIO_SUBSAMPLE:-}"

[[ -n "${PAT_READS}" && -n "${MAT_READS}" ]] \
  || die "no parental reads for ${SAMPLE}: set pat_reads and mat_reads in ${SAMPLES_TSV}, or PAT_READS and MAT_READS"
for f in ${PAT_READS} ${MAT_READS}; do
  require_file "${f}" "parental reads"
done
if [[ -n "${TRIO_SUBSAMPLE}" ]]; then
  awk -v f="${TRIO_SUBSAMPLE}" 'BEGIN { exit !(f > 0 && f < 1) }' \
    || die "TRIO_SUBSAMPLE must be a fraction between 0 and 1, not '${TRIO_SUBSAMPLE}'"
fi

export PATH="${TRIO_BIN}:${PATH}"
require_cmds yak meryl "${SAMTOOLS}"

MERQURY="${MERQURY:-$(ls -d "${TRIO_BIN}"/../share/merqury* 2>/dev/null | head -1)}"
[[ -s "${MERQURY}/trio/hapmers.sh" ]] \
  || die "merqury not found in ${TRIO_BIN}/../share - set MERQURY to its directory"
export MERQURY

mkdir -p "${TRIO_DIR}"
cd "${TRIO_DIR}"
log "${SAMPLE}: trio k-mers in ${TRIO_DIR}"

# as_fastq <role> <file>... - prints the FASTQ files to count: FASTQ as given,
# a BAM or CRAM converted once to <role>.<name>.fq.gz in TRIO_DIR
as_fastq() {
  local role="$1"; shift
  local f fq sub=()
  [[ -n "${TRIO_SUBSAMPLE}" ]] && sub=(-s "${SEED}${TRIO_SUBSAMPLE#0}")
  for f in "$@"; do
    case "${f}" in
      *.bam|*.cram)
        fq="${TRIO_DIR}/${role}.$(basename "${f%.*}")${TRIO_SUBSAMPLE:+.sub${TRIO_SUBSAMPLE}}.fq.gz"
        if [[ ! -s "${fq}" ]]; then
          log "${role}: $(basename "${f}") -> $(basename "${fq}")${TRIO_SUBSAMPLE:+ (fraction ${TRIO_SUBSAMPLE})}"
          "${SAMTOOLS}" view -u -@ "${TRIO_THREADS}" -F 0x900 ${sub[@]+"${sub[@]}"} "${f}" \
            | "${SAMTOOLS}" fastq -@ "${TRIO_THREADS}" -n - 2>"${fq}.log" \
            | "${BGZIP}" -@ "${TRIO_THREADS}" -c >"${fq}.tmp"
          mv -f "${fq}.tmp" "${fq}"
        fi
        printf '%s\n' "${fq}"
        ;;
      *)
        printf '%s\n' "${f}"
        ;;
    esac
  done
}

# yak counts in two passes, so it reads the same input twice
yak_db() {   # <name> <reads...>
  local name="$1"; shift
  [[ -s "${name}.yak" ]] && { log "skip ${name}.yak"; return 0; }
  log "${name}.yak from $# files"
  yak count -k31 -b37 -t "${TRIO_THREADS}" -o "${name}.yak.tmp" <(cat "$@") <(cat "$@")
  mv -f "${name}.yak.tmp" "${name}.yak"
}

# verkko wants the k-mers counted in homopolymer-compressed space
meryl_db() {   # <name> <reads...>
  local name="$1"; shift
  [[ -d "${name}.k30.meryl" ]] && { log "skip ${name}.k30.meryl"; return 0; }
  log "${name}.k30.meryl from $# files"
  rm -rf "${name}.k30.tmp.meryl"
  meryl count compress k=30 threads="${TRIO_THREADS}" memory="${TRIO_MEMORY}" \
      "$@" output "${name}.k30.tmp.meryl"
  mv "${name}.k30.tmp.meryl" "${name}.k30.meryl"
}

# through files, not process substitution, so a failed conversion stops here
as_fastq pat ${PAT_READS} >pat.fastq.list
as_fastq mat ${MAT_READS} >mat.fastq.list
mapfile -t pat_fq <pat.fastq.list
mapfile -t mat_fq <mat.fastq.list
[[ ${#pat_fq[@]} -gt 0 && ${#mat_fq[@]} -gt 0 ]] || die "no parental FASTQ to count"

yak_db pat "${pat_fq[@]}"
yak_db mat "${mat_fq[@]}"
meryl_db pat "${pat_fq[@]}"
meryl_db mat "${mat_fq[@]}"

# k-mers found in one parent only; without a child database hapmers.sh names
# them <input>.only.meryl
if [[ -d pat.k30.only.meryl && -d mat.k30.only.meryl ]]; then
  log "skip hapmers"
else
  log "hapmers"
  bash "${MERQURY}/trio/hapmers.sh" mat.k30.meryl pat.k30.meryl >hapmers.log 2>&1 \
    || die "hapmers.sh failed, see ${TRIO_DIR}/hapmers.log"
  [[ -d pat.k30.only.meryl && -d mat.k30.only.meryl ]] \
    || die "hapmers.sh wrote no *.only.meryl, see ${TRIO_DIR}/hapmers.log"
fi

log "trio k-mers ready in ${TRIO_DIR}"
