#!/usr/bin/env bash
# One-off preparation of the GIAB/GA4GH stratification BEDs.
#
#   bash bin/prepare_strata.sh                      GRCh38
#   REFERENCE_NAME=T2T bash bin/prepare_strata.sh   T2T
#
# Extracts the tarballs, writes regions/strata.tsv, and for the current
# reference produces
#   regions/prepared/<ref>_<stratum>.bed          stratum intersected with the
#                                                 confident BED - the region all
#                                                 the metrics are computed over
#   regions/truth/<ref>_<stratum>.phased.vcf.gz   truth restricted to it
#
# The two references do not ship the same strata: GRCh38 has GenomeSpecific and
# Ancestry, CHM13 has rDNA and Telomere. strata.tsv is the union, and whatever
# the current reference is missing is reported and skipped.

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

mkdir -p "${STRATA_DIR}" "${PREPARED_DIR}" "${STRATA_TRUTH_DIR}"

TIER1="notinalldifficultregions alldifficultregions segdups lowmappabilityall
       AllTandemRepeatsandHomopolymers_slop5 satellites_slop5 MHC chrX_nonPAR"

TIER2="SimpleRepeat_homopolymer_4to6_slop5 SimpleRepeat_homopolymer_7to11_slop5
       SimpleRepeat_homopolymer_ge12_slop5 SimpleRepeat_homopolymer_ge21_slop5
       AllTandemRepeats_le50bp_slop5 AllTandemRepeats_51to200bp_slop5
       AllTandemRepeats_201to10000bp_slop5 AllTandemRepeats_ge10001bp_slop5
       gclt25orgt65_slop50"

tier_of() {
  for s in ${TIER1}; do [[ "$1" == "${s}" ]] && { echo 1; return; }; done
  for s in ${TIER2}; do [[ "$1" == "${s}" ]] && { echo 2; return; }; done
  echo 3
}

# the archives are 1.3 and 0.9 GB, so unpack once
for ref in GRCh38 CHM13; do
  [[ -d "${STRATA_DIR}/${ref}@all" ]] && continue
  log "unpacking ${ref} stratifications"
  tar xzf "${STRATA_DIR}/genome-stratifications-${ref}@all.tar.gz" -C "${STRATA_DIR}"
done

if [[ ! -s "${STRATA_LIST}" ]]; then
  log "building $(basename "${STRATA_LIST}")"
  {
    printf 'tier\tcategory\tstratum\n'
    for ref in GRCh38 CHM13; do
      # GenomeSpecific ships strata for HG001 through HG007; only ours matter
      find "${STRATA_DIR}/${ref}@all" -name "${ref}_*.bed.gz" \
        | awk '!/GenomeSpecific/ || /HG002/' \
        | sed "s|${STRATA_DIR}/${ref}@all/||; s|/${ref}_|\t|; s|\.bed\.gz$||"
    done | sort -u -k2,2 | while IFS=$'\t' read -r category stratum; do
      printf '%s\t%s\t%s\n' "$(tier_of "${stratum}")" "${category}" "${stratum}"
    done | sort -k1,1n -k3,3
  } >"${STRATA_LIST}"
fi
log "strata in list: $(( $(wc -l <"${STRATA_LIST}") - 1 ))"

prepare_one() {
  set -euo pipefail

  local stratum="$1"
  local bed; bed="$(stratum_bed "${stratum}")"
  local truth; truth="$(stratum_truth "${stratum}")"

  [[ -s "${truth}.tbi" ]] && { log "skip ${stratum}"; return 0; }

  local src
  src="$(find "${STRATA_DIR}/${STRATA_REF}@all" -name "${STRATA_REF}_${stratum}.bed.gz" | head -1)"
  [[ -n "${src}" ]] || { warn "${STRATA_REF}: нет страты ${stratum}"; return 0; }

  gzip -dc "${src}" >"${bed}.raw"
  bed_intersect "${CONFIDENT_BED}" "${bed}.raw" >"${bed}"
  rm -f "${bed}.raw"

  "${BCFTOOLS}" view -T "${bed}" -Oz -o "${truth}" "${TRUTH_PHASED}"
  "${TABIX}" -f -p vcf "${truth}"

  log "  $(printf '%-52s %12s bp  %9s het SNP' \
      "${stratum}" "$(bed_bp "${bed}")" "$("${BCFTOOLS}" index -n "${truth}")")"
}

export -f prepare_one stratum_bed stratum_truth bed_intersect bed_bp log warn die
export STRATA_DIR STRATA_REF PREPARED_DIR STRATA_TRUTH_DIR REFERENCE_NAME \
       CONFIDENT_BED TRUTH_PHASED BCFTOOLS TABIX

log "preparing ${REFERENCE_NAME} (strata named ${STRATA_REF})"
strata_names \
  | "${PARALLEL}" --joblog "${LOGS_DIR}/parallel_prepare_strata_${REFERENCE_NAME}.log" \
        -j "${REGIONS_JOBS}" prepare_one {} \
  || true

log "prepared: $(ls "${PREPARED_DIR}/${REFERENCE_NAME}"_*.bed 2>/dev/null | wc -l) strata"
