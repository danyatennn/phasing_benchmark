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

    [[ -n "${src_depth}" ]] \
        || die "no measured source depth for ${platform}"

    local key
    key="${platform}_${REFERENCE_NAME}_${SAMPLE}_${depth}x_wgs"

    local out_dir
    out_dir="${DATA_DIR}/${key}"

    local out
    out="${out_dir}/input.bam"

    mkdir -p "${out_dir}"

    log "${key}: starting"

    rm -f "${out}" "${out}.bai"

    local frac
    frac=$(awk \
        -v target="${depth}" \
        -v source="${src_depth}" \
        'BEGIN {
            if (source <= 0) exit 1
            printf "%.6f", target/source
        }' \
        || die "${key}: bad source depth '${src_depth}' for ${platform}")

    if awk -v f="${frac}" 'BEGIN { exit !(f >= 0.999999) }'
    then
        warn "${key}: source depth is ${src_depth}x; using full BAM"
        ln -s "${src}" "${out}"
        ln -s "${src}.bai" "${out}.bai"
        frac="1.000000"
    else
        local sarg
        sarg=$(awk \
            -v seed="${SEED}" \
            -v frac="${frac}" \
            'BEGIN { printf "%d.%06d",
                seed,
                int(frac * 1000000 + 0.5)
            }')

        log "${key}: source=${src_depth}x target=${depth}x fraction=${frac}"

        "${SAMTOOLS}" view \
            -@ "${THREADS_PER_JOB}" \
            -b \
            -s "${sarg}" \
            -o "${out}" \
            "${src}"

        "${SAMTOOLS}" index \
            -@ "${THREADS_PER_JOB}" \
            "${out}"
    fi

    local got
    got="$(measure_depth "${out}")"

    log "${key}: achieved ${got}x"

    printf '%s\t%s\t%s\t%s\n' \
        "${platform}" \
        "${depth}" \
        "${frac}" \
        "${got}" \
        > "${TMP_DEPTH_DIR}/${platform}_${depth}.tsv"
}


export -f make_one
export -f source_bam
export -f source_depth
export -f measure_depth
export -f require_file
export -f log
export -f warn
export -f die

export MOSDEPTH
export MOSDEPTH_THREADS="${DEPTH_THREADS}"
export DATA_DIR
export REFERENCE_NAME
export SAMPLE
export SAMTOOLS
export THREADS_PER_JOB
export PARALLEL_JOBS
export TMP_DEPTH_DIR
export SEED
export DEPTH_TSV
export ONT_HG002_V14_5KHZ_NO_ALT_BAM
export ONT_HG002_dorado_NO_ALT_BAM
export PB_HG002_20K_NO_ALT_BAM
export HG002_REVIO_BAM

echo "[INFO] Whole-genome downsampling"
require_file "${DEPTH_TSV}" "source depth table"
"${PARALLEL}" \
    --joblog "${DATA_DIR}/parallel_downsample.log" \
    -j "${PARALLEL_JOBS}" \
    make_one {1} {2} \
    ::: "${PLATFORMS[@]}" \
    ::: "${DEPTHS[@]}"


printf 'platform\tdepth\tfraction\tachieved_depth\n' \
    > "${ACHIEVED_TSV}"

find "${TMP_DEPTH_DIR}" \
    -type f \
    -name '*.tsv' \
    -print0 \
| sort -z \
| xargs -0 -r cat \
>> "${ACHIEVED_TSV}"

rm -rf "${TMP_DEPTH_DIR}"


log "downsampling done"
log "depth summary: ${ACHIEVED_TSV}"