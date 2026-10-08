#!/usr/bin/env bash

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${PHASE_DIR}" "${LOGS_DIR}/phase"

phase_one() {
  set -euo pipefail

  local caller="$1" platform="$2" depth="$3" phaser="$4"
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}" "${phaser}")"
  local out="${PHASE_DIR}/${name}.vcf.gz"

  is_done "${out}" && { log "skip ${name}"; return 0; }

  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local logf="${LOGS_DIR}/phase/${name}.log"

  local vcf_in
  if [[ "${caller}" == truth ]]; then
    vcf_in="${TRUTH_UNPHASED}"
  else
    vcf_in="${VCF_DIR}/$(vcf_name "${platform}" "${depth}" "${caller}").het_snps.vcf.gz"
  fi

  log "${name}"

  case "${phaser}" in
    whatshap)
      "${WHATSHAP}" phase \
          --output "${out}" \
          --reference "${REF}" \
          --sample "${SAMPLE}" \
          --ignore-read-groups \
          "${vcf_in}" "${bam}" \
        >"${logf}" 2>&1 \
        || die "whatshap failed for ${name}, see ${logf}"
      "${TABIX}" -f -p vcf "${out}"
      ;;

    margin)
      local tmp; tmp="$(mktemp -d)"
      "${BCFTOOLS}" view "${vcf_in}" >"${tmp}/in.vcf"

      "${MARGIN}" phase \
          "${bam}" "${REF}" "${tmp}/in.vcf" "$(margin_params "${platform}")" \
          -t "${PHASE_THREADS}" \
          -o "${tmp}/out" \
          -M \
        >"${logf}" 2>&1 \
        || { rm -rf "${tmp}"; die "margin failed for ${name}, see ${logf}"; }

      "${BGZIP}" -f -c "${tmp}/out.phased.vcf" >"${out}"
      "${TABIX}" -f -p vcf "${out}"
      rm -rf "${tmp}"
      ;;

    hapcut2)
      local tmp; tmp="$(mktemp -d)"
      "${BCFTOOLS}" view "${vcf_in}" >"${tmp}/in.vcf"

      "${HAPCUT2_BIN}/extractHAIRS" \
          "$(hapcut2_flag "${platform}")" 1 \
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

      "${BGZIP}" -f -c "${tmp}/hap.phased.vcf" >"${out}"
      "${TABIX}" -f -p vcf "${out}"
      rm -rf "${tmp}"
      ;;

    longphase)
      local tmp; tmp="$(mktemp -d)"
      "${BCFTOOLS}" view "${vcf_in}" >"${tmp}/in.vcf"

      "${LONGPHASE}" phase \
          -s "${tmp}/in.vcf" \
          -b "${bam}" \
          -r "${REF}" \
          -t "${PHASE_THREADS}" \
          -o "${tmp}/out" \
          "$(longphase_flag "${platform}")" \
        >"${logf}" 2>&1 \
        || { rm -rf "${tmp}"; die "longphase failed for ${name}, see ${logf}"; }

      "${BGZIP}" -f -c "${tmp}/out.vcf" >"${out}"
      "${TABIX}" -f -p vcf "${out}"
      rm -rf "${tmp}"
      ;;
  esac
}

export -f phase_one vcf_name sub_bam bam_key longphase_flag margin_params \
          hapcut2_flag is_done log warn die
export PHASE_DIR VCF_DIR LOGS_DIR DATA_DIR REF SAMPLE REFERENCE_NAME REGION_NAME \
       TRUTH_UNPHASED WHATSHAP LONGPHASE MARGIN MARGIN_PARAMS_DIR HAPCUT2_BIN \
       BCFTOOLS BGZIP TABIX PHASE_THREADS

log "phasing: ${#CALLERS[@]} callers x ${#PLATFORMS[@]} platforms x ${#DEPTHS[@]} depths x ${#PHASERS[@]} phasers"

"${PARALLEL}" \
    --joblog "${LOGS_DIR}/parallel_phase.log" \
    -j "${PHASE_JOBS}" \
    phase_one {1} {2} {3} {4} \
    ::: "${CALLERS[@]}" \
    ::: "${PLATFORMS[@]}" \
    ::: "${DEPTHS[@]}" \
    ::: "${PHASERS[@]}"

log "phasing done"
