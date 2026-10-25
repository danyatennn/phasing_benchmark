#!/usr/bin/env bash

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${PHASE_DIR}" "${LOGS_DIR}/phase" "${TMP_DIR}"

phase_one() {
  set -euo pipefail

  local caller="$1" platform="$2" depth="$3" phaser="$4"
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
  local out="${PHASE_DIR}/${name}.vcf.gz"

  is_done "${out}" && { log "skip ${name}"; return 0; }
  phaser_applies "${phaser}" "${platform}" || return 0

  local lock="${out}.lock"
  mkdir "${lock}" 2>/dev/null || { log "busy elsewhere ${name} ($(cat "${lock}/owner" 2>/dev/null))"; return 0; }
  printf '%s %s %s\n' "$(hostname)" "$$" "$(date '+%F %H:%M')" >"${lock}/owner"
  trap "rm -f '${lock}/owner'; rmdir '${lock}' 2>/dev/null" EXIT INT TERM HUP

  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local logf="${LOGS_DIR}/phase/${name}.log"

  # truthindel / truthall: experiment A (bin/make_truth_variants.sh), the truth
  # SNPs plus other het variants; LongPhase and HapCUT2 are told to use the
  # indels, WhatsHap (2.x) and HiPhase phase every variant they are given
  local vcf_in indels=0
  case "${caller}" in
    truth)      vcf_in="${TRUTH_UNPHASED}" ;;
    truthindel) vcf_in="${TRUTH_INDEL_UNPHASED}"; indels=1 ;;
    truthall)   vcf_in="${TRUTH_ALL_UNPHASED}"; indels=1 ;;
    *)          vcf_in="${VCF_DIR}/$(vcf_name "${platform}" "${depth}" "${caller}").het_snps.vcf.gz" ;;
  esac
  [[ -s "${vcf_in}" ]] || die "no phaser input for ${name}: ${vcf_in}"

  log "${name}"

  case "${phaser}" in
    whatshap)
      "${WHATSHAP}" phase \
          --output "${out}" \
          --reference "${REF}" \
          --sample "${SAMPLE}" \
          --ignore-read-groups \
          --internal-downsampling "${WHATSHAP_DOWNSAMPLING}" \
          "${vcf_in}" "${bam}" \
        >"${logf}" 2>&1 \
        || die "whatshap failed for ${name}, see ${logf}"
      "${TABIX}" -f -p vcf "${out}"
      ;;

    margin)
      local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
      "${BCFTOOLS}" view "${vcf_in}" >"${tmp}/in.vcf"

      # no -p unless MARGIN_DEPTH is set: Margin downsamples to maxDepth
      local depth_cap=()
      [[ -n "${MARGIN_DEPTH}" ]] && depth_cap=(-p "${MARGIN_DEPTH}")

      "${MARGIN}" phase \
          "${bam}" "${REF}" "${tmp}/in.vcf" "$(margin_params "${platform}")" \
          -t "${PHASE_THREADS}" \
          ${depth_cap[@]+"${depth_cap[@]}"} \
          -o "${tmp}/out" \
          -M \
        >"${logf}" 2>&1 \
        || { rm -rf "${tmp}"; die "margin failed for ${name}, see ${logf}"; }

      "${BGZIP}" -f -c "${tmp}/out.phased.vcf" >"${out}"
      "${TABIX}" -f -p vcf "${out}"
      rm -rf "${tmp}"
      ;;

    hapcut2)
      local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
      "${BCFTOOLS}" view "${vcf_in}" >"${tmp}/in.vcf"

      local hairs_indels=()
      [[ "${indels}" == 1 ]] && hairs_indels=(--indels 1)
      "${HAPCUT2_BIN}/extractHAIRS" \
          "$(hapcut2_flag "${platform}")" 1 \
          ${hairs_indels[@]+"${hairs_indels[@]}"} \
          --bam "${bam}" \
          --VCF "${tmp}/in.vcf" \
          --ref "${REF}" \
          --out "${tmp}/frags" \
        >"${logf}" 2>&1 \
        || { rm -rf "${tmp}"; die "extractHAIRS failed for ${name}, see ${logf}"; }

      "${HAPCUT2_BIN}/HAPCUT2" \
          --fragments "${tmp}/frags" \
          --VCF "${tmp}/in.vcf" \
          --output "${tmp}/hap" \
          --outvcf 1 \
        >>"${logf}" 2>&1 \
        || { rm -rf "${tmp}"; die "hapcut2 failed for ${name}, see ${logf}"; }

      local hap_vcf; hap_vcf="$(ls "${tmp}"/hap.phased.* 2>/dev/null | head -1)"
      [[ -n "${hap_vcf}" ]] \
        || { rm -rf "${tmp}"; die "hapcut2 wrote no phased VCF for ${name}, see ${logf}"; }

      "${BGZIP}" -f -c "${hap_vcf}" >"${out}"
      "${TABIX}" -f -p vcf "${out}"
      rm -rf "${tmp}"
      ;;

    longphase)
      local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"
      "${BCFTOOLS}" view "${vcf_in}" >"${tmp}/in.vcf"

      local lp_indels=()
      [[ "${indels}" == 1 ]] && lp_indels=(--indels)
      "${LONGPHASE}" phase \
          -s "${tmp}/in.vcf" \
          -b "${bam}" \
          -r "${REF}" \
          -t "${PHASE_THREADS}" \
          -o "${tmp}/out" \
          "$(longphase_flag "${platform}")" \
          ${lp_indels[@]+"${lp_indels[@]}"} \
        >"${logf}" 2>&1 \
        || { rm -rf "${tmp}"; die "longphase failed for ${name}, see ${logf}"; }

      "${BGZIP}" -f -c "${tmp}/out.vcf" >"${out}"
      "${TABIX}" -f -p vcf "${out}"
      rm -rf "${tmp}"
      ;;

    hiphase)
      # PacBio only (phaser_applies); default settings. It reads the bgzipped,
      # indexed input as is and writes a bgzipped VCF with its own index.
      local tmp; tmp="$(mktemp -d -p "${TMP_DIR}")"

      "${HIPHASE}" \
          --bam "${bam}" \
          --vcf "${vcf_in}" \
          --reference "${REF}" \
          --sample-name "${SAMPLE}" \
          --ignore-read-groups \
          --threads "${PHASE_THREADS}" \
          --output-vcf "${tmp}/out.vcf.gz" \
        >"${logf}" 2>&1 \
        || { rm -rf "${tmp}"; die "hiphase failed for ${name}, see ${logf}"; }

      mv -f "${tmp}/out.vcf.gz" "${out}"
      "${TABIX}" -f -p vcf "${out}"
      rm -rf "${tmp}"
      ;;

    *)
      die "unknown phaser '${phaser}'"
      ;;
  esac
}

export -f phase_one vcf_name sub_bam bam_key longphase_flag margin_params \
          hapcut2_flag phaser_applies is_done log warn die
export PHASE_DIR VCF_DIR LOGS_DIR DATA_DIR REF SAMPLE REFERENCE_NAME REGION_NAME \
       TRUTH_UNPHASED TRUTH_INDEL_UNPHASED TRUTH_ALL_UNPHASED WHATSHAP LONGPHASE MARGIN MARGIN_PARAMS_DIR HAPCUT2_BIN HIPHASE \
       BCFTOOLS BGZIP TABIX PHASE_THREADS TMP_DIR WHATSHAP_DOWNSAMPLING MARGIN_DEPTH

log "phasing: ${#CALLERS[@]} callers x ${#PLATFORMS[@]} platforms x ${#DEPTHS[@]} depths x ${#PHASERS[@]} phasers"

joblog="${LOGS_DIR}/parallel_phase_${RUN_TAG}.log"
keep_joblog "${joblog}"

"${PARALLEL}" \
    --joblog "${joblog}" \
    -j "${PHASE_JOBS}" \
    phase_one {1} {2} {3} {4} \
    ::: "${CALLERS[@]}" \
    ::: "${PLATFORMS[@]}" \
    ::: "${DEPTHS[@]}" \
    ::: "${PHASERS[@]}" \
  || true

# a run whose caller never produced a VCF fails here; let the rest finish
failed="$(awk 'NR > 1 && $7 != 0 { n++ } END { print n + 0 }' "${joblog}")"
[[ "${failed}" -eq 0 ]] || warn "${failed} runs failed, logs in ${LOGS_DIR}/phase"

log "phasing done"
