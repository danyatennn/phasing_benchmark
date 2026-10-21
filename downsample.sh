#!/usr/bin/env bash

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

require_vars SAMPLE SAMTOOLS MOSDEPTH PARALLEL
require_cmds "${SAMTOOLS}" "${MOSDEPTH}" "${PARALLEL}"

mkdir -p "${DATA_DIR}"
rm -rf "${TMP_DEPTH_DIR}"
mkdir -p "${TMP_DEPTH_DIR}"

make_one() {
    set -euo pipefail
    local platform="$1"
    local depth="$2"
    local src
    src="$(source_bam "${platform}")"

    require_file "${src}" "source BAM for ${platform}"
    require_file "${src}.bai" "source BAM index for ${platform}"
    
    local src_depth
    src_depth="$(source_depth "${platform}")"

    [[ -n "${src_depth}" ]] || die "no measured source depth for ${platform}"

    local key
    key="${platform}_${REFERENCE_NAME}_${SAMPLE}_${depth}x_wgs"

    local out_dir
    out_dir="${DATA_DIR}/${key}"

    local out
    out="${out_dir}/input.bam"

    mkdir -p "${out_dir}"

    # Finished = the BAM and an index written after it
    if [[ -s "${out}" && -s "${out}.bai" && ! "${out}" -nt "${out}.bai" ]]; then
        local kept
        kept="$(awk -F'\t' -v p="${platform}" -v d="${depth}" \
            'NR > 1 && $1 == p && $2 == d { print; exit }' "${ACHIEVED_TSV}" 2>/dev/null || true)"
        if [[ -z "${kept}" ]]; then
            kept="$(printf '%s\t%s\t%s\t%s' "${platform}" "${depth}" "NA" "$(measure_depth "${out}")")"
        fi
        log "${key}: exists, kept"
        printf '%s\n' "${kept}" > "${TMP_DEPTH_DIR}/${platform}_${depth}.tsv"
        return 0
    fi

    log "${key}: starting"

    rm -f "${out}" "${out}.bai"
    # written under a temporary name and moved into place once indexed
    local tmp_bam="${out_dir}/input.tmp.bam"
    rm -f "${tmp_bam}" "${tmp_bam}.bai"

    local frac
    frac=$(awk -v target="${depth}" -v source="${src_depth}" \
        'BEGIN {
            if (source <= 0) exit 1
            printf "%.6f", target/source
        }' || die "${key}: bad source depth '${src_depth}' for ${platform}")

    if awk -v f="${frac}" 'BEGIN { exit !(f >= 0.999999) }'
    then
        warn "${key}: source depth is ${src_depth}x; using full BAM"
        frac="1.000000"

        "${SAMTOOLS}" view -@ "${THREADS_PER_JOB}" -b -x HP -x PS -x PC -o "${tmp_bam}" "${src}"
        "${SAMTOOLS}" index -@ "${THREADS_PER_JOB}" "${tmp_bam}"
    else
        local sarg
        sarg=$(awk -v seed="${SEED}" -v frac="${frac}" \
            'BEGIN { printf "%d.%06d",
                seed,
                int(frac * 1000000 + 0.5)
            }')
        log "${key}: source=${src_depth}x target=${depth}x fraction=${frac}"
        "${SAMTOOLS}" view -@ "${THREADS_PER_JOB}" -b -x HP -x PS -x PC -s "${sarg}" -o "${tmp_bam}" "${src}"
        "${SAMTOOLS}" index -@ "${THREADS_PER_JOB}" "${tmp_bam}"
    fi
    # the BAM first, then its index, so the index is never older than the BAM
    mv -f "${tmp_bam}" "${out}"
    mv -f "${tmp_bam}.bai" "${out}.bai"
    local got
    got="$(measure_depth "${out}")"
    log "${key}: achieved ${got}x"
    printf '%s\t%s\t%s\t%s\n' "${platform}" "${depth}" "${frac}" "${got}" > "${TMP_DEPTH_DIR}/${platform}_${depth}.tsv"
}


export -f make_one source_bam sample_path source_depth measure_depth require_file log warn die

export MOSDEPTH DATA_DIR REFERENCE_NAME T2T_SRC_DIR SAMPLE SAMPLES_TSV \
       SAMTOOLS THREADS_PER_JOB PARALLEL_JOBS TMP_DEPTH_DIR SEED \
       DEPTH_TSV ACHIEVED_TSV ONT_HG002_V14_5KHZ_NO_ALT_BAM \
       ONT_HG002_dorado_NO_ALT_BAM PB_HG002_20K_NO_ALT_BAM HG002_REVIO_BAM

export MOSDEPTH_THREADS="${DEPTH_THREADS}"

log "whole-genome downsampling against ${REFERENCE_NAME}"
require_file "${DEPTH_TSV}" "source depth table"
"${PARALLEL}" \
    --joblog "${DOWNSAMPLE_JOBLOG}" \
    -j "${PARALLEL_JOBS}" \
    make_one {1} {2} \
    ::: "${PLATFORMS[@]}" \
    ::: "${DEPTHS[@]}"


printf 'platform\tdepth\tfraction\tachieved_depth\n' > "${ACHIEVED_TSV}"

find "${TMP_DEPTH_DIR}" -type f -name '*.tsv' -print0 | sort -z \
| xargs -0 -r cat >> "${ACHIEVED_TSV}"

rm -rf "${TMP_DEPTH_DIR}"

log "downsampling done"
log "depth summary: ${ACHIEVED_TSV}"