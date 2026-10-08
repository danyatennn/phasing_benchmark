#!/usr/bin/env bash
# vcf/<platform>_<ref>_<sample>_<depth>x_<region>_<caller>.vcf.gz - all calls
# vcf/<platform>_<ref>_<sample>_<depth>x_<region>_<caller>.het_snps.vcf.gz - phaser input

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${VCF_DIR}" "${CLAIR3_DIR}" "${LOGS_DIR}/call"

call_one() {
  set -euo pipefail

  local caller="$1" platform="$2" depth="$3"
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}")"
  local out="${VCF_DIR}/${name}.vcf.gz"
  local het="${VCF_DIR}/${name}.het_snps.vcf.gz"

  is_done "${het}" && { log "skip ${name}"; return 0; }

  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local logf="${LOGS_DIR}/call/${name}.log"
  local work="${CLAIR3_DIR}/${name}"

  log "${name}: calling"
  rm -rf "${work}"
  mkdir -p "${work}"

  case "${caller}" in
    clair3)
      "${SING}" exec --nv --cleanenv --env TMPDIR=/tmp \
          -B "${DATA_DIR},$(dirname "${REF}"),$(dirname "${CONFIDENT_BED}"),${CLAIR3_DIR}" \
          "${CLAIR3_SIF}" \
          /opt/bin/run_clair3.sh \
            --bam_fn="${bam}" \
            --ref_fn="${REF}" \
            --threads="${CLAIR3_THREADS}" \
            --platform="$(clair3_platform "${platform}")" \
            --model_path="/opt/models/$(clair3_model "${platform}")" \
            --sample_name="${SAMPLE}" \
            --bed_fn="${CONFIDENT_BED}" \
            --output="${work}" \
            --remove_intermediate_dir \
            --use_gpu \
          >"${logf}" 2>&1 \
        || die "clair3 failed for ${name}, see ${logf}"
      cp "${work}/merge_output.vcf.gz" "${out}"
      ;;

    deepvariant)
      "${SING}" exec --cleanenv --env TMPDIR=/tmp \
          -B "${DATA_DIR},$(dirname "${REF}"),$(dirname "${CONFIDENT_BED}"),${CLAIR3_DIR}" \
          "${DEEPVARIANT_SIF}" \
          /opt/deepvariant/bin/run_deepvariant \
            --model_type="$(deepvariant_model "${platform}")" \
            --ref="${REF}" \
            --reads="${bam}" \
            --regions="${CONFIDENT_BED}" \
            --sample_name="${SAMPLE}" \
            --output_vcf="${work}/output.vcf.gz" \
            --num_shards="${CALL_THREADS}" \
            --intermediate_results_dir="${work}/tmp" \
          >"${logf}" 2>&1 \
        || die "deepvariant failed for ${name}, see ${logf}"
      cp "${work}/output.vcf.gz" "${out}"
      ;;

    longcalld)
      # --region-file keeps longcallD inside the same regions as the others
      "${LONGCALLD}" call \
          -t "${CALL_THREADS}" \
          "$(longcalld_preset "${platform}")" \
          --region-file "${CONFIDENT_BED}" \
          "${REF}" "${bam}" \
        2>"${logf}" \
        | "${BGZIP}" -c >"${out}" \
        || die "longcalld failed for ${name}, see ${logf}"
      ;;
  esac

  "${TABIX}" -f -p vcf "${out}"

  "${BCFTOOLS}" view -f PASS -m2 -M2 -v snps -g het -Oz -o "${het}" "${out}"
  "${TABIX}" -f -p vcf "${het}"

  log "${name}: $("${BCFTOOLS}" index -n "${het}") het SNPs"
}

export -f call_one vcf_name sub_bam bam_key is_done log warn die \
          clair3_platform clair3_model deepvariant_model longcalld_preset
export VCF_DIR CLAIR3_DIR LOGS_DIR DATA_DIR REF CONFIDENT_BED SAMPLE \
       REFERENCE_NAME REGION_NAME SING CLAIR3_SIF DEEPVARIANT_SIF LONGCALLD \
       CLAIR3_THREADS CALL_THREADS BCFTOOLS TABIX BGZIP \
       CLAIR3_MODEL_ont4k CLAIR3_MODEL_ont5k CLAIR3_MODEL_hifi CLAIR3_MODEL_revio

for caller in ${CALL_ONLY:-${VARIANT_CALLERS[@]}}; do
  jobs="$(call_jobs "${caller}")"
  log "${caller}: ${#PLATFORMS[@]} platforms x ${#DEPTHS[@]} depths, ${jobs} at a time"
  "${PARALLEL}" \
      --joblog "${LOGS_DIR}/parallel_call_${caller}.log" \
      -j "${jobs}" \
      call_one "${caller}" {1} {2} \
      ::: "${PLATFORMS[@]}" \
      ::: "${DEPTHS[@]}"
done

log "calling done"
