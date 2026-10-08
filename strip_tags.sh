#!/usr/bin/env bash

# Removes the HP/PS/PC haplotype tags from the downsampled BAMs
set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

strip_one() {
    set -euo pipefail

    local platform="$1"
    local depth="$2"

    local key="${platform}_${REFERENCE_NAME}_${SAMPLE}_${depth}x_wgs"
    local bam="${DATA_DIR}/${key}/input.bam"
    local tmp="${bam}.tmp"

    require_file "${bam}" "downsampled BAM"

    log "${key}: stripping HP/PS/PC"

    rm -f "${tmp}" "${tmp}.bai"

    "${SAMTOOLS}" view \
        -@ "${THREADS_PER_JOB}" \
        -b \
        -x HP -x PS -x PC \
        -o "${tmp}" \
        "${bam}"

    "${SAMTOOLS}" index -@ "${THREADS_PER_JOB}" "${tmp}"

    mv -f "${tmp}"     "${bam}"
    mv -f "${tmp}.bai" "${bam}.bai"

    log "${key}: done"
}

export -f strip_one log die require_file
export SAMTOOLS DATA_DIR REFERENCE_NAME SAMPLE THREADS_PER_JOB

log "stripping HP/PS/PC from ${#PLATFORMS[@]} x ${#DEPTHS[@]} BAMs"

"${PARALLEL}" \
    --joblog "${DATA_DIR}/parallel_strip_tags.log" \
    -j "${PARALLEL_JOBS}" \
    strip_one {1} {2} \
    ::: "${PLATFORMS[@]}" \
    ::: "${DEPTHS[@]}"

log "HP/PS/PC removed from all BAMs"
