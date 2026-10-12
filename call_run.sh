#!/usr/bin/env bash
# vcf/<platform>_<ref>_<sample>_<depth>x_<region>_<caller>.vcf.gz - all calls
# vcf/<platform>_<ref>_<sample>_<depth>x_<region>_<caller>.het_snps.vcf.gz - phaser input

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${VCF_DIR}" "${CLAIR3_DIR}" "${LOGS_DIR}/call" "${SINGULARITY_TMPDIR}" "${SCRATCH_DIR}"

call_one() {
  set -euo pipefail

  local caller="$1" platform="$2" depth="$3"
  local name; name="$(vcf_name "${platform}" "${depth}" "${caller}")"
  local out="${VCF_DIR}/${name}.vcf.gz"
  local het="${VCF_DIR}/${name}.het_snps.vcf.gz"

  is_done "${het}" && { log "skip ${name}"; return 0; }

  local lock="${het}.lock"
  mkdir "${lock}" 2>/dev/null || { log "busy elsewhere ${name}"; return 0; }
  trap "rmdir '${lock}' 2>/dev/null" EXIT INT TERM HUP

  local bam; bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
  local logf="${LOGS_DIR}/call/${name}.log"
  local work="${CLAIR3_DIR}/${name}"
  local scratch="${SCRATCH_DIR}/${name}"

  local gpu=()
  [[ -n "${GPU_DEVICE}" ]] && gpu=(--env "CUDA_VISIBLE_DEVICES=${GPU_DEVICE}")

  log "${name}: calling"
  rm -rf "${work}" "${scratch}"
  mkdir -p "${work}" "${scratch}"

  case "${caller}" in
    clair3)
      "${SING}" exec --nv --cleanenv ${gpu[@]+"${gpu[@]}"} --env TMPDIR="${scratch}" \
          -B "${DATA_DIR},$(dirname "${REF}"),$(dirname "${CONFIDENT_BED}"),${CLAIR3_DIR},${SCRATCH_DIR}" \
          "${CLAIR3_SIF}" \
          /opt/bin/run_clair3.sh \
            --bam_fn="${bam}" \
            --ref_fn="${REF}" \
            --threads="${CLAIR3_THREADS}" \
            --platform="$(clair3_platform "${platform}")" \
            --model_path="/opt/models/$(clair3_model "${platform}")" \
            --sample_name="${SAMPLE}" \
            --bed_fn="${CONFIDENT_BED}" \
            --output="${scratch}" \
            --remove_intermediate_dir \
            --use_gpu \
          >"${logf}" 2>&1 \
        || { rm -rf "${scratch}"; die "clair3 failed for ${name}, see ${logf}"; }
      cp "${scratch}/merge_output.vcf.gz" "${out}"
      rm -rf "${scratch}"
      ;;

    deepvariant)
      local nv=()
      [[ "${DEEPVARIANT_USE_GPU}" == 1 ]] && nv=(--nv)
      "${SING}" exec ${nv[@]+"${nv[@]}"} --cleanenv ${gpu[@]+"${gpu[@]}"} --env TMPDIR="${scratch}" \
          -B "${DATA_DIR},$(dirname "${REF}"),$(dirname "${CONFIDENT_BED}"),${CLAIR3_DIR},${SCRATCH_DIR}" \
          "${DEEPVARIANT_SIF}" \
          /opt/deepvariant/bin/run_deepvariant \
            --model_type="$(deepvariant_model "${platform}")" \
            --ref="${REF}" \
            --reads="${bam}" \
            --regions="${CONFIDENT_BED}" \
            --sample_name="${SAMPLE}" \
            --output_vcf="${scratch}/output.vcf.gz" \
            --num_shards="${CALL_THREADS}" \
            --intermediate_results_dir="${scratch}/dv" \
          >"${logf}" 2>&1 \
        || { rm -rf "${scratch}"; die "deepvariant failed for ${name}, see ${logf}"; }
      cp "${scratch}/output.vcf.gz" "${out}"
      rm -rf "${scratch}"
      ;;

    longcalld)
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

  if [[ "$("${BCFTOOLS}" query -l "${out}")" != "${SAMPLE}" ]]; then
    log "${name}: renaming sample $("${BCFTOOLS}" query -l "${out}") -> ${SAMPLE}"
    printf '%s\n' "${SAMPLE}" >"${work}/sample.txt"
    "${BCFTOOLS}" reheader -s "${work}/sample.txt" -o "${out}.tmp" "${out}"
    mv -f "${out}.tmp" "${out}"
    "${TABIX}" -f -p vcf "${out}"
  fi

  "${BCFTOOLS}" view -f PASS -m2 -M2 -v snps -g het -Oz -o "${het}" "${out}"
  "${TABIX}" -f -p vcf "${het}"

  log "${name}: $("${BCFTOOLS}" index -n "${het}") het SNPs"
}

export -f call_one vcf_name sub_bam bam_key is_done log warn die \
          clair3_platform clair3_model deepvariant_model longcalld_preset
export VCF_DIR CLAIR3_DIR LOGS_DIR DATA_DIR REF CONFIDENT_BED SAMPLE \
       REFERENCE_NAME REGION_NAME SING CLAIR3_SIF DEEPVARIANT_SIF LONGCALLD \
       CLAIR3_THREADS CALL_THREADS DEEPVARIANT_USE_GPU GPU_DEVICE SCRATCH_DIR BCFTOOLS TABIX BGZIP \
       SINGULARITY_TMPDIR SINGULARITY_CACHEDIR \
       CLAIR3_MODEL_ont4k CLAIR3_MODEL_ont5k CLAIR3_MODEL_hifi CLAIR3_MODEL_revio

for caller in ${CALL_ONLY:-${VARIANT_CALLERS[@]}}; do
  jobs="$(call_jobs "${caller}")"
  log "${caller}: ${#PLATFORMS[@]} platforms x ${#DEPTHS[@]} depths, ${jobs} at a time"
  joblog="${LOGS_DIR}/parallel_call_${caller}_${RUN_TAG}.log"
  keep_joblog "${joblog}"

  "${PARALLEL}" \
      --joblog "${joblog}" \
      -j "${jobs}" \
      call_one "${caller}" {1} {2} \
      ::: "${PLATFORMS[@]}" \
      ::: "${DEPTHS[@]}" \
    || true

  failed="$(awk 'NR > 1 && $7 != 0 { n++ } END { print n + 0 }' "${joblog}")"
  [[ "${failed}" -eq 0 ]] \
    || warn "${caller}: ${failed} cells failed, logs in ${LOGS_DIR}/call"
done

log "calling done"
