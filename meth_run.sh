#!/usr/bin/env bash
# Methylation-aware phasing, truth input, on the cells whose reads carry 5mC
# (MM/ML tags; a cell without them is skipped). Each method is written as a
# phaser of its own, so eval_run.sh and the analyses score it like any run:
#
#   lpmod        LongPhase with 5mC: longphase modcall, then longphase phase
#                --mod-file, default settings, on the cell's BAM. The phased
#                5mC VCF is kept as phased/<run>.mod.vcf.gz.
#   mp<phaser>   MethPhaser on top of <phaser>'s VCF of the same cell, default
#                settings:
#                  reads: primary alignments on chr1-22, 5mC calls only
#                         (modkit adjust-mods drops 5hmC, 6mA, ...), tagged with
#                         the phaser's VCF (whatshap haplotag, one chromosome
#                         per job: it tags each chromosome on its own anyway)
#                  blocks: whatshap stats --gtf of that VCF
#                  meth_phaser_parallel, then meth_phaser_post_processing
#                The VCF post-processing writes is the run's VCF; the haplotype
#                it gave each read is kept as metrics/<run>.methtags.tsv.
#
#   phased/<platform>_<ref>_<sample>_30x_wgs_truth_{lpmod,mpwhatshap,...}.vcf.gz
#
# MethPhaser works on autosomes only (random X inactivation), so its VCFs hold
# chr1-22 - the metrics we report anyway. Every output is locked while being
# written and kept once done, so several hosts can share the work.
#
#   bash meth_run.sh                                 HG002, GRCh38, 30x
#   REFERENCE_NAME=T2T bash meth_run.sh
#   SAMPLE=HG001 bash meth_run.sh
#   METH_METHODS=lpmod bash meth_run.sh              one method
#   METH_BASE="whatshap longphase" bash meth_run.sh  MethPhaser on these only
# Then score the runs:
#   PHASING_CALLERS=truth PHASING_DEPTHS=30 \
#   PHASING_PHASERS="lpmod mpwhatshap mplongphase mpmargin mphapcut2 mphiphase" bash eval_run.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

METH_DEPTHS=(${METH_DEPTHS:-30})
METH_METHODS=(${METH_METHODS:-lpmod methphaser})
METH_BASE=(${METH_BASE:-${PHASERS[*]}})
# a cell's 5mC reads are dropped once MethPhaser ran on all of these, whatever
# METH_BASE this host was given: another host may still be on the others
MP_ALL=(whatshap longphase margin hapcut2 hiphase)
CHROMS=(chr{1..22})
# reads are sampled here to look at the MM tags: four chromosomes, away from
# centromeres on GRCh38 and CHM13 alike
SAMPLE_REGIONS=(chr1:100000000-100200000 chr6:40000000-40200000 chr11:80000000-80200000 chr20:40000000-40200000)

# meth_phaser_parallel calls methphasing by its bare name, and both are python
# scripts of the methphaser env
export PATH="${METHPHASER_DIR}:${METH_ENV_BIN}:${PATH}"
MODKIT="${METH_ENV_BIN}/modkit"
require_cmds meth_phaser_parallel meth_phaser_post_processing methphasing "${MODKIT}"
mkdir -p "${PHASE_DIR}" "${METRICS_DIR}" "${METH_DIR}" "${LOGS_DIR}/meth" "${TMP_DIR}"

# haplotag_part <vcf> <bam> <chrom> <out> - whatshap haplotag on one chromosome.
# whatshap reads, tags and writes one chromosome at a time and carries nothing
# over (no BX tags here, supplementary alignments are gone), so the parts put
# together in the order of the BAM header are the whole-genome output
haplotag_part() {
  "${WHATSHAP}" haplotag --reference "${REF}" --sample "${SAMPLE}" --ignore-read-groups \
      --regions "$3" -o "$4" "$1" "$2" >"$4.log" 2>&1
}
export -f haplotag_part
export WHATSHAP REF SAMPLE

# sample_reads <bam> - up to 1000 reads from each of SAMPLE_REGIONS
sample_reads() {
  local r
  for r in "${SAMPLE_REGIONS[@]}"; do
    "${SAMTOOLS}" view "$1" "${r}" 2>/dev/null | head -1000 || true
  done
}

# has_5mc <bam> - do the reads carry 5mC calls (C+m in MM)?
has_5mc() {
  local n
  n="$(sample_reads "$1" | grep -cE $'\tMM:Z:[^\t]*C\\+[a-z]*m' || true)"
  [[ "${n}" -gt 0 ]]
}

# lock <path> - take <path>.lock or return 1; released by the subshell's trap
lock() {
  mkdir "$1.lock" 2>/dev/null || { log "busy elsewhere $(basename "$1") ($(cat "$1.lock/owner" 2>/dev/null))"; return 1; }
  printf '%s %s %s\n' "$(hostname)" "$$" "$(date '+%F %H:%M')" >"$1.lock/owner"
  trap "rm -f '$1.lock/owner'; rmdir '$1.lock' 2>/dev/null" EXIT INT TERM HUP
}

# mp_vcf <platform> <depth> <phaser> - the VCF of a run (phaser or mp<phaser>)
mp_vcf() {
  echo "${PHASE_DIR}/$(vcf_name "$1" "$2" truth "$3").vcf.gz"
}

# lpmod_one <platform> <depth>
lpmod_one() (
  set -euo pipefail
  local platform="$1" depth="$2"
  local name; name="$(vcf_name "${platform}" "${depth}" truth lpmod)"
  local out="${PHASE_DIR}/${name}.vcf.gz"
  is_done "${out}" && { log "skip ${name}"; return 0; }
  lock "${out}" || return 0

  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local logf="${LOGS_DIR}/meth/${name}.log"
  local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
  "${BCFTOOLS}" view "${TRUTH_UNPHASED}" >"${tmp}/in.vcf"

  log "${name}: longphase modcall"
  "${LONGPHASE}" modcall -b "${bam}" -r "${REF}" -t "${METH_THREADS}" -o "${tmp}/modcall" \
    >"${logf}" 2>&1 || { rm -rf "${tmp}"; die "longphase modcall failed for ${name}, see ${logf}"; }
  local mod; mod="$(ls "${tmp}"/modcall*.vcf 2>/dev/null | head -1)"
  [[ -s "${mod}" ]] || { rm -rf "${tmp}"; die "modcall wrote no VCF for ${name}, see ${logf}"; }

  log "${name}: longphase phase --mod-file"
  "${LONGPHASE}" phase -s "${tmp}/in.vcf" -b "${bam}" -r "${REF}" --mod-file "${mod}" \
      -t "${METH_THREADS}" -o "${tmp}/out" "$(longphase_flag "${platform}")" \
    >>"${logf}" 2>&1 || { rm -rf "${tmp}"; die "longphase phase failed for ${name}, see ${logf}"; }
  [[ -s "${tmp}/out.vcf" ]] || { rm -rf "${tmp}"; die "longphase wrote no out.vcf for ${name}, see ${logf}"; }

  local phased_mod; phased_mod="$(ls "${tmp}"/out*mod*.vcf 2>/dev/null | head -1)"
  if [[ -s "${phased_mod}" ]]; then
    "${BGZIP}" -f -c "${phased_mod}" >"${PHASE_DIR}/${name}.mod.vcf.gz"
  else
    warn "${name}: no phased 5mC VCF among $(ls "${tmp}" | tr '\n' ' ')"
  fi
  "${BGZIP}" -f -c "${tmp}/out.vcf" >"${out}.tmp"
  "${TABIX}" -f -p vcf "${out}.tmp"
  mv -f "${out}.tmp.tbi" "${out}.tbi"
  mv -f "${out}.tmp" "${out}"
  rm -rf "${tmp}"
  log "${name}: done"
)

# mp_todo <platform> <depth> - is a MethPhaser run of the cell (METH_BASE)
# still to do, with its phaser's VCF there to start from?
mp_todo() {
  local platform="$1" depth="$2" phaser
  for phaser in "${METH_BASE[@]}"; do
    phaser_applies "${phaser}" "${platform}" || continue
    is_done "$(mp_vcf "${platform}" "${depth}" "mp${phaser}")" && continue
    is_done "$(mp_vcf "${platform}" "${depth}" "${phaser}")" && return 0
  done
  return 1
}

# prep_reads <platform> <depth> - primary alignments on CHROMS with 5mC calls
# only, the reads MethPhaser is given; once per cell
prep_reads() (
  set -euo pipefail
  local platform="$1" depth="$2"
  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local work="${METH_DIR}/$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local out="${work}/reads.5mc.bam"
  mkdir -p "${work}"
  [[ -s "${out}.bai" ]] && newer_than "${out}" "${bam}" && return 0
  lock "${out}" || return 1
  local logf="${LOGS_DIR}/meth/$(basename "${work}").reads.log"

  log "$(basename "${work}"): primary alignments on ${#CHROMS[@]} chromosomes"
  # 2308 = unmapped + secondary + supplementary: MethPhaser's pysam reading of
  # MM fails on the others
  "${SAMTOOLS}" view -@ "${METH_THREADS}" -b -F 2308 -o "${work}/primary.bam" "${bam}" "${CHROMS[@]}"
  "${SAMTOOLS}" index "${work}/primary.bam"

  # every modification code but 5mC (m) the sampled reads carry: MM is
  # "C+h?,...;C+m?,...;" or "C+hm?,..."; one-letter codes may be joined, a
  # ChEBI code is a number. MethPhaser reads only C+m, so a code missed here
  # only leaves its probability where the basecaller put it
  local codes
  codes="$(sample_reads "${work}/primary.bam" | grep -oE $'\tMM:Z:[^\t]*' \
    | cut -c7- | tr ';' '\n' | grep -oE '^[ACGTUN][-+][a-z0-9]+' | cut -c3- \
    | awk '{ if ($0 ~ /^[0-9]+$/) print; else for (i = 1; i <= length($0); i++) print substr($0, i, 1) }' \
    | sort -u | grep -vx m || true)"
  local cur="${work}/primary.bam" c
  for c in ${codes}; do
    log "$(basename "${work}"): modkit adjust-mods --ignore ${c}"
    "${MODKIT}" adjust-mods "${cur}" "${work}/no_${c}.bam" --ignore "${c}" \
      >>"${logf}" 2>&1 || die "modkit failed on $(basename "${work}"), see ${logf}"
    [[ "${cur}" == "${work}/primary.bam" ]] || rm -f "${cur}"
    cur="${work}/no_${c}.bam"
  done
  if [[ "${cur}" == "${work}/primary.bam" ]]; then
    mv -f "${work}/primary.bam" "${out}"
  else
    mv -f "${cur}" "${out}"
    rm -f "${work}/primary.bam"
  fi
  rm -f "${work}/primary.bam.bai"
  "${SAMTOOLS}" index "${out}"
)

# mp_complete <work> <logf> - did every methphasing job of meth_phaser_parallel
# finish? It starts them without looking at how they end (subprocess.run). A
# job writes methphaser/<chrom>/<first>_<last>.csv when it is done, and the
# jobs of a chromosome are cut from its n GTF lines as meth_phaser_parallel
# does: ceil(n / threads) blocks each, overlapping by one
mp_complete() {
  local work="$1" logf="$2" chrom want have ok=0 crashed
  # a traceback from multiprocessing's exit clean-up is not a crash: its
  # temporary directory cannot be removed on NFS while files in it are still
  # open ("Directory not empty"), which Python reports and ignores
  crashed="$(awk '/^Traceback/ { t = 1; next }
      t { if ($0 !~ /multiprocessing\/util\.py.*_run_finalizers/) n++; t = 0 }
      END { print n + 0 }' "${logf}")"
  if [[ "${crashed}" -gt 0 ]]; then
    warn "$(basename "${work}"): ${crashed} methphasing jobs crashed, see ${logf}"
    ok=1
  fi
  while read -r chrom want; do
    have="$(ls "${work}/methphaser/${chrom}"/*.csv 2>/dev/null | wc -l)"
    [[ "${have}" -eq "${want}" ]] && continue
    warn "$(basename "${work}") ${chrom}: ${have} of ${want} methphasing jobs wrote a result"
    ok=1
  done < <(awk -F'\t' -v t="${METH_THREADS}" '{ n[$1]++ }
      END { for (c in n) { bs = int((n[c] + t - 1) / t); print c, int((n[c] - 1 + bs - 1) / bs) } }' \
      "${work}/blocks.gtf")
  return "${ok}"
}

# phased_sites <vcf> - chrom:pos of the records with a phased genotype
phased_sites() {
  "${BCFTOOLS}" query -f '%CHROM:%POS\t[%GT]\n' "$1" | awk -F'\t' '$2 ~ /\|/ { print $1 }' | sort -u
}

# mp_one <platform> <depth> <phaser> - MethPhaser on top of <phaser>'s VCF
mp_one() (
  set -euo pipefail
  local platform="$1" depth="$2" phaser="$3"
  local base; base="$(mp_vcf "${platform}" "${depth}" "${phaser}")"
  local name; name="$(vcf_name "${platform}" "${depth}" truth "mp${phaser}")"
  local out="${PHASE_DIR}/${name}.vcf.gz"
  is_done "${out}" && { log "skip ${name}"; return 0; }
  is_done "${base}" || { warn "${name}: no ${phaser} VCF to start from"; return 0; }
  lock "${out}" || return 0

  local cell; cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local reads="${METH_DIR}/${cell}/reads.5mc.bam"
  local work="${METH_DIR}/${cell}/mp${phaser}"
  local logf="${LOGS_DIR}/meth/${name}.log"
  local tagged="${work}/tagged.bam"
  # the tagged reads of an earlier attempt are kept if they are complete (the
  # index is written last) and newer than the phaser's VCF and the reads
  local reuse=""
  [[ -s "${tagged}.bai" ]] && newer_than "${tagged}" "${base}" && newer_than "${tagged}" "${reads}" && reuse=1
  if [[ -n "${reuse}" ]]; then
    find "${work}" -mindepth 1 -maxdepth 1 ! -name tagged.bam ! -name tagged.bam.bai -exec rm -rf {} +
  else
    rm -rf "${work}"
    mkdir -p "${work}"
  fi
  : >"${logf}"

  "${BCFTOOLS}" view -r "$(IFS=,; echo "${CHROMS[*]}")" -Oz -o "${work}/base.vcf.gz" "${base}"
  "${TABIX}" -f -p vcf "${work}/base.vcf.gz"

  if [[ -n "${reuse}" ]]; then
    log "${name}: reads tagged with the ${phaser} VCF already there"
  else
    log "${name}: haplotag the reads with the ${phaser} VCF, ${METH_THREADS} chromosomes at a time"
    local order=() parts=() c
    order=($("${SAMTOOLS}" view -H "${reads}" \
      | awk -F'\t' '$1 == "@SQ" { for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) print substr($i, 4) }' \
      | grep -xE 'chr([1-9]|1[0-9]|2[0-2])' || true))
    [[ ${#order[@]} -eq ${#CHROMS[@]} ]] || die "${name}: ${#order[@]} of chr1-22 in the header of ${reads}"
    for c in "${order[@]}"; do parts+=("${work}/tagged.${c}.bam"); done
    "${PARALLEL}" -j "${METH_THREADS}" --halt now,fail=1 \
        haplotag_part "${work}/base.vcf.gz" "${reads}" {} "${work}/tagged.{}.bam" ::: "${order[@]}" \
      || { cat "${work}"/tagged.*.bam.log >>"${logf}" 2>/dev/null; die "whatshap haplotag failed for ${name}, see ${logf}"; }
    cat "${work}"/tagged.*.bam.log >>"${logf}"
    # the parts go in the order of the header, the order of the sorted input
    "${SAMTOOLS}" cat -o "${work}/tagged.tmp.bam" "${parts[@]}"
    "${SAMTOOLS}" index "${work}/tagged.tmp.bam"
    # whatshap writes every read it fetches, tagged or not: all of them are back
    local n_in n_out
    n_in="$("${SAMTOOLS}" idxstats "${reads}" | awk '{ n += $3 } END { print n + 0 }')"
    n_out="$("${SAMTOOLS}" idxstats "${work}/tagged.tmp.bam" | awk '{ n += $3 } END { print n + 0 }')"
    [[ "${n_in}" -eq "${n_out}" ]] || die "${name}: ${n_out} of ${n_in} reads in the tagged BAM"
    rm -f "${parts[@]}" "${work}"/tagged.*.bam.log
    mv -f "${work}/tagged.tmp.bam.bai" "${tagged}.bai"
    mv -f "${work}/tagged.tmp.bam" "${tagged}"
  fi
  "${WHATSHAP}" stats --sample "${SAMPLE}" --gtf "${work}/blocks.gtf" "${work}/base.vcf.gz" >>"${logf}" 2>&1 \
    || die "whatshap stats --gtf failed for ${name}, see ${logf}"

  log "${name}: meth_phaser_parallel"
  meth_phaser_parallel -b "${work}/tagged.bam" -r "${REF}" -g "${work}/blocks.gtf" \
      -vc "${work}/base.vcf.gz" -t "${METH_THREADS}" -o "${work}/methphaser" >>"${logf}" 2>&1 \
    || die "meth_phaser_parallel failed for ${name}, see ${logf}"
  mp_complete "${work}" "${logf}" || die "${name}: MethPhaser did not finish, see ${logf}"
  log "${name}: meth_phaser_post_processing"
  meth_phaser_post_processing -ib "${work}/tagged.bam" -if "${work}/methphaser" \
      -ov "${work}/out.vcf" -ob "${work}/out" -vc "${work}/base.vcf.gz" -t "${METH_THREADS}" >>"${logf}" 2>&1 \
    || die "meth_phaser_post_processing failed for ${name}, see ${logf}"
  [[ -s "${work}/out.vcf" ]] || die "${name}: no out.vcf, see ${logf}"

  # the VCF as MethPhaser wrote it, sorted for the index. Two records at one
  # position would make the sort order matter for whatshap, so it is reported
  "${BCFTOOLS}" sort -T "${work}" -Oz -o "${out}.tmp" "${work}/out.vcf" >>"${logf}" 2>&1
  "${TABIX}" -f -p vcf "${out}.tmp"
  local dup
  dup="$("${BCFTOOLS}" query -f '%CHROM\t%POS\n' "${out}.tmp" | uniq -d | wc -l)"
  [[ "${dup}" -eq 0 ]] || warn "${name}: ${dup} positions with more than one record"
  # every phased record of the base VCF comes back once (each GTF line is
  # fetched whole); post-processing leaves out unphased records before the
  # first block and after the last two of a chromosome, which only changes
  # whatshap stats' variant count, so that is reported
  local lost unphased_lost
  lost="$(comm -23 <(phased_sites "${work}/base.vcf.gz") <(phased_sites "${out}.tmp") | wc -l)"
  [[ "${lost}" -eq 0 ]] || { rm -f "${out}.tmp" "${out}.tmp.tbi"; die "${name}: ${lost} phased variants of ${phaser} missing from MethPhaser's VCF, see ${logf}"; }
  unphased_lost=$(( $("${BCFTOOLS}" index -n "${work}/base.vcf.gz") - $("${BCFTOOLS}" index -n "${out}.tmp") ))
  [[ "${unphased_lost}" -le 0 ]] || warn "${name}: ${unphased_lost} unphased records of ${phaser} not written by MethPhaser"

  # the haplotype MethPhaser gave each read (primary alignments only): HP is
  # relative to the read's block in the run's VCF. A read whatshap had tagged
  # keeps its PS (MethPhaser flips HP with the block), a read MethPhaser
  # assigned has none ('.'), so its span places it; a read on two blocks may
  # be written twice
  local tags="${METRICS_DIR}/${name}.methtags.tsv"
  { printf 'read\thp\tps\tchrom\tstart\tend\n'
    for b in "${work}"/out*.bam; do
      [[ -e "${b}" ]] || continue
      "${SAMTOOLS}" view -F 2304 "${b}" | awk -F'\t' '{
          hp = ""; ps = "."
          for (i = 12; i <= NF; i++) {
            if ($i ~ /^HP:i:/) hp = substr($i, 6)
            else if ($i ~ /^PS:i:/) ps = substr($i, 6)
          }
          if (hp == "") next
          len = 0; cigar = $6
          while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
            op = substr(cigar, RLENGTH, 1)
            if (op ~ /[MDN=X]/) len += substr(cigar, 1, RLENGTH - 1)
            cigar = substr(cigar, RLENGTH + 1)
          }
          print $1 "\t" hp "\t" ps "\t" $3 "\t" $4 "\t" $4 + len - 1 }'
    done
  } >"${tags}.tmp"
  mv -f "${tags}.tmp" "${tags}"

  mv -f "${out}.tmp.tbi" "${out}.tbi"
  mv -f "${out}.tmp" "${out}"
  # the BAMs are large and every number is in the VCF and the tag list now;
  # MethPhaser's per-junction tables in methphaser/ stay
  rm -f "${work}"/*.bam "${work}"/*.bai "${work}/base.vcf.gz" "${work}/base.vcf.gz.tbi"
  log "${name}: done, $(( $(wc -l <"${tags}") - 1 )) tagged reads"
)

# the cell's 5mC reads go once MethPhaser ran on every phaser of MP_ALL; while
# one is missing (or being written somewhere) they are kept
drop_reads() {
  local platform="$1" depth="$2" phaser
  for phaser in "${MP_ALL[@]}"; do
    phaser_applies "${phaser}" "${platform}" || continue
    is_done "$(mp_vcf "${platform}" "${depth}" "mp${phaser}")" || return 0
  done
  rm -f "${METH_DIR}/$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"/reads.5mc.bam{,.bai}
}

log "${SAMPLE} ${REFERENCE_NAME}: ${METH_METHODS[*]} on ${PLATFORMS[*]} x ${METH_DEPTHS[*]}x"
for platform in "${PLATFORMS[@]}"; do
  for depth in "${METH_DEPTHS[@]}"; do
    bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
    [[ -s "${bam}" ]] || { warn "no cell BAM ${bam}"; continue; }
    has_5mc "${bam}" || { log "skip $(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}"): no 5mC tags"; continue; }
    for method in "${METH_METHODS[@]}"; do
      case "${method}" in
        lpmod)
          lpmod_one "${platform}" "${depth}" || warn "lpmod failed on ${platform} ${depth}x"
          ;;
        methphaser)
          mp_todo "${platform}" "${depth}" || { log "skip MethPhaser on ${platform} ${depth}x: nothing to do"; continue; }
          prep_reads "${platform}" "${depth}" || { warn "5mC reads of ${platform} ${depth}x not ready (busy or failed)"; continue; }
          for phaser in "${METH_BASE[@]}"; do
            phaser_applies "${phaser}" "${platform}" || continue
            mp_one "${platform}" "${depth}" "${phaser}" || warn "MethPhaser on ${phaser} failed on ${platform} ${depth}x"
          done
          drop_reads "${platform}" "${depth}"
          ;;
        *)
          die "METH_METHODS: unknown method '${method}' (lpmod, methphaser)"
          ;;
      esac
    done
  done
done
log "meth done"
