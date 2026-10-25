PROJECT_ROOT="/autofs/bal36/dten/phasing_benchmark"

# One GIAB sample per run. HG002 keeps the layout the grid was built in; every
# other sample gets the same tree under samples/<sample>/, so nothing it writes
# can land among the HG002 results. Its BAMs, truth and BED come from
# samples.tsv (sample_path in lib.sh).
SAMPLE="${SAMPLE:-HG002}"
SAMPLES_TSV="${PROJECT_ROOT}/samples.tsv"
if [[ "${SAMPLE}" == HG002 ]]; then
  SAMPLE_ROOT="${PROJECT_ROOT}"
else
  SAMPLE_ROOT="${PROJECT_ROOT}/samples/${SAMPLE}"
fi

DATA_DIR="${SAMPLE_ROOT}/data"

DERIVED_DIR="${SAMPLE_ROOT}/derived"
if [[ "${SAMPLE}" == HG002 ]]; then
  DEPTHS=(${PHASING_DEPTHS:-10 20 30 40 50})
  PLATFORMS=(${PHASING_PLATFORMS:-ont4k ont5k hifi revio})
else
  # the other samples are benchmarked at 30x, on whatever platforms they have
  DEPTHS=(${PHASING_DEPTHS:-30})
  PLATFORMS=(${PHASING_PLATFORMS:-$(sample_platforms "${SAMPLE}")})
  for _p in ont4k ont5k hifi revio; do
    export "SAMPLE_BAM_${_p}=$(sample_path "${SAMPLE}" "bam_${_p}")"
  done
fi
PHASING_BIN="/autofs/bal36/dten/miniforge3/envs/phasing/bin"
BCFTOOLS="${PHASING_BIN}/bcftools"
TABIX="${PHASING_BIN}/tabix"
BGZIP="${PHASING_BIN}/bgzip"

SEED=1
PARALLEL_JOBS=6
THREADS_PER_JOB=6

DEPTH_THREADS=4

REFERENCE_NAME="${REFERENCE_NAME:-GRCh38}"
if [[ "${REFERENCE_NAME}" == "GRCh38" ]]; then
  _ref_sfx=""
else
  _ref_sfx="_${REFERENCE_NAME}"
fi

DEPTH_TSV="${DATA_DIR}/source_depth${_ref_sfx}.tsv"
ACHIEVED_TSV="${DATA_DIR}/achieved_depth${_ref_sfx}.tsv"
TMP_DEPTH_DIR="${DATA_DIR}/.achieved_depth${_ref_sfx}_tmp"
DOWNSAMPLE_JOBLOG="${DATA_DIR}/parallel_downsample${_ref_sfx}.log"

CLAIR3_DIR="${SAMPLE_ROOT}/clair3"
LOGS_DIR="${SAMPLE_ROOT}/logs"
VCF_DIR="${SAMPLE_ROOT}/vcf"
PHASE_DIR="${SAMPLE_ROOT}/phased"
METRICS_DIR="${SAMPLE_ROOT}/metrics"

# The paper is about the autosomes (chr1-22): the v4.2.1 truth of the other
# samples has nothing else, and the assemblies are scored on them only. With
# STRATA_AUTOSOMES=1 the strata are cut to chr1-22 and everything built from
# them goes to *_auto directories, next to the whole-genome ones, which are
# left as they are.
STRATA_AUTOSOMES="${STRATA_AUTOSOMES:-1}"
if [[ "${STRATA_AUTOSOMES}" == 1 ]]; then
  _strata_sfx="_auto"
else
  _strata_sfx=""
fi
METRICS_STRATA_DIR="${SAMPLE_ROOT}/metrics_strata${_strata_sfx}"

# the stratification archives and the list of strata are the same for every
# sample; what they are cut to (the sample's confident BED) is not
REGIONS_DIR="${PROJECT_ROOT}/regions"
STRATA_DIR="${REGIONS_DIR}/strata"
STRATA_LIST="${REGIONS_DIR}/strata.tsv"
PREPARED_DIR="${SAMPLE_ROOT}/regions/prepared${_strata_sfx}"
STRATA_TRUTH_DIR="${SAMPLE_ROOT}/regions/truth${_strata_sfx}"

if [[ "${REFERENCE_NAME}" == "GRCh38" ]]; then
  STRATA_REF="GRCh38"
else
  STRATA_REF="CHM13"
fi

REGIONS_JOBS="${REGIONS_JOBS:-12}"
REGIONS_SHARD="${REGIONS_SHARD:-}"
STRATA_ONLY="${STRATA_ONLY:-}"
STRATA_TIER="${STRATA_TIER:-}"
TMP_DIR="${SAMPLE_ROOT}/tmp"


SCRATCH_DIR="${SCRATCH_DIR:-${TMP_DIR}}"
export SINGULARITY_TMPDIR="${SINGULARITY_TMPDIR:-${SCRATCH_DIR}/singularity}"
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-${SCRATCH_DIR}/singularity_cache}"

SING="/autofs/bal33/zxzheng/env/miniconda2/envs/singularity-env/bin/singularity"
CLAIR3_SIF="${PROJECT_ROOT}/clair3_gpu.sif"
CLAIR3_THREADS="${CLAIR3_THREADS:-32}"

if [[ "${REFERENCE_NAME}" == "GRCh38" ]]; then
  REF="${REF_GRCH38}"
  CONFIDENT_BED="${CONFIDENT_BED:-${HG002_GRCH38_V5Q_BED}}"
  TRUTH_SRC_VCF="${TRUTH_SRC_VCF:-${HG002_GRCH38_V5Q_TRUTH_VCF}}"
else
  REF="${REF_CHM13}"
  CONFIDENT_BED="${CONFIDENT_BED:-${HG002_CHM13_V5Q_BED}}"
  TRUTH_SRC_VCF="${TRUTH_SRC_VCF:-${HG002_CHM13_V5Q_TRUTH_VCF}}"
fi

# Any other sample takes its truth from samples.tsv, never the HG002 one: an
# empty path stops make_truth.sh and the callers with a clear message.
if [[ "${SAMPLE}" != HG002 ]]; then
  CONFIDENT_BED="$(sample_path "${SAMPLE}" "bed_${REFERENCE_NAME}")"
  TRUTH_SRC_VCF="$(sample_path "${SAMPLE}" "truth_${REFERENCE_NAME}")"
  # the BED is read as plain text everywhere (awk, Clair3, longcallD); a
  # gzipped one from samples.tsv is unpacked once into the sample's derived/
  if [[ "${CONFIDENT_BED}" == *.gz && -s "${CONFIDENT_BED}" ]]; then
    _plain_bed="${DERIVED_DIR}/$(basename "${CONFIDENT_BED%.gz}")"
    if [[ ! -s "${_plain_bed}" ]]; then
      mkdir -p "${DERIVED_DIR}"
      gzip -dc "${CONFIDENT_BED}" >"${_plain_bed}.$$" && mv -f "${_plain_bed}.$$" "${_plain_bed}"
    fi
    CONFIDENT_BED="${_plain_bed}"
  fi
fi

REGION_NAME="${REGION_NAME:-wgs}"

CLAIR3_MODEL_ont4k="${CLAIR3_MODEL_ont4k:-r1041_e82_400bps_sup_v410}"
CLAIR3_MODEL_ont5k="${CLAIR3_MODEL_ont5k:-r1041_e82_400bps_sup_v500}"
CLAIR3_MODEL_hifi="${CLAIR3_MODEL_hifi:-hifi_sequel2}"
CLAIR3_MODEL_revio="${CLAIR3_MODEL_revio:-hifi_revio}"

VARIANT_CALLERS=(clair3 deepvariant longcalld)
CALLERS=(${PHASING_CALLERS:-truth "${VARIANT_CALLERS[@]}"})
# hiphase runs on the PacBio platforms only (phaser_applies in lib.sh)
PHASERS=(${PHASING_PHASERS:-whatshap longphase margin hapcut2 hiphase})

RUN_TAG="${REFERENCE_NAME}"
RUN_TAG="${RUN_TAG}_$(IFS=-; printf '%s' "${PLATFORMS[*]}")"
RUN_TAG="${RUN_TAG}_$(IFS=-; printf '%s' "${DEPTHS[*]}")"
RUN_TAG="${RUN_TAG}_$(IFS=-; printf '%s' "${CALLERS[*]}")"

CALL_THREADS="${CALL_THREADS:-32}"

# a hung caller blocks its host until someone notices, which cost us seven hours
# on bal34. On timeout the cell fails, the stage carries on, the lock is released
# by the trap, and the next run picks the cell up again.
CALL_TIMEOUT="${CALL_TIMEOUT:-12h}"

CLAIR3_JOBS="${CLAIR3_JOBS:-1}"
GPU_DEVICE="${GPU_DEVICE:-}"

DEEPVARIANT_USE_GPU="${DEEPVARIANT_USE_GPU:-1}"
if [[ "${DEEPVARIANT_USE_GPU}" == 1 ]]; then
  DEEPVARIANT_JOBS="${DEEPVARIANT_JOBS:-1}"
else
  DEEPVARIANT_JOBS="${DEEPVARIANT_JOBS:-2}"
fi
LONGCALLD_JOBS="${LONGCALLD_JOBS:-2}"

SOFTWARE_DIR="${PROJECT_ROOT}/software"

WHATSHAP="${PHASING_BIN}/whatshap"
LONGPHASE="${SOFTWARE_DIR}/longphase_linux-x64"
LONGCALLD="${SOFTWARE_DIR}/longcallD-v0.0.11_x64-linux/longcallD"
if [[ "${DEEPVARIANT_USE_GPU}" == 1 ]]; then
  DEEPVARIANT_SIF="${SOFTWARE_DIR}/deepvariant_1.10.0-gpu.sif"
else
  DEEPVARIANT_SIF="${SOFTWARE_DIR}/deepvariant_1.10.0.sif"
fi
MARGIN="${SOFTWARE_DIR}/margin/build/margin"
MARGIN_PARAMS_DIR="${SOFTWARE_DIR}/margin/params/phase"
HIFIASM="${SOFTWARE_DIR}/hifiasm/hifiasm"
HAPCUT2_BIN="/autofs/bal36/dten/miniforge3/envs/hapcut2/bin"
HIPHASE="${HIPHASE:-/autofs/bal36/dten/miniforge3/envs/hiphase/bin/hiphase}"

# meth_run.sh: methylation-aware phasing on the cells whose reads carry 5mC
# (MM/ML tags), truth input. MethPhaser from its GitHub main (0.0.3 on bioconda
# writes a broken VCF); its env also has modkit.
METH_ENV_BIN="${METH_ENV_BIN:-/autofs/bal36/dten/miniforge3/envs/methphaser/bin}"
METHPHASER_DIR="${METHPHASER_DIR:-${SOFTWARE_DIR}/methphaser}"
METH_DIR="${SAMPLE_ROOT}/meth"
METH_THREADS="${METH_THREADS:-16}"

WHATSHAP_DOWNSAMPLING="${WHATSHAP_DOWNSAMPLING:-15}"

# Margin downsamples the reads to maxDepth from its params file (32 in
# allParams.phase_vcf.*.json) unless -p overrides it. Empty: no -p, Margin as
# its authors set it up. A number: passed as -p. Runs up to October 2026 used 30.
MARGIN_DEPTH="${MARGIN_DEPTH:-}"

PHASE_JOBS="${PHASE_JOBS:-6}"
PHASE_THREADS="${PHASE_THREADS:-8}"
EVAL_JOBS="${EVAL_JOBS:-6}"

HAPLOTAG_REGION="${HAPLOTAG_REGION:-}"

# Read-level accuracy (readacc_run.sh, readacc_anchor.sh):
#   READACC_MIN_SNPS   a block of fewer phased het SNPs orients nothing, and
#                      its reads count as untagged. Margin writes a phase set
#                      even for a single SNP; a read tagged by it is scored on
#                      one allele, not on phase. whatshap stats counts blocks
#                      of two or more, so 2 matches it.
#   READACC_AUTOSOMES  1: reads and blocks on chr1-22 only.
# With 1 and 0 readacc_run.sh writes the original <run>.readacc.tsv; any other
# setting goes to its own file, see readacc_suffix in lib.sh.
READACC_MIN_SNPS="${READACC_MIN_SNPS:-2}"
READACC_AUTOSOMES="${READACC_AUTOSOMES:-1}"

# per-read tables (bin/read_tables.sh, bin/allele_matrix.sh) and what is
# derived from them
READS_DIR="${SAMPLE_ROOT}/reads"
METRICS_READS_DIR="${SAMPLE_ROOT}/metrics_reads"
SWITCHBED_DIR="${SAMPLE_ROOT}/metrics_switchbed"

# python with pysam (the phasing env, which whatshap needs anyway) and the one
# the figures are drawn with
PHASING_PYTHON="${PHASING_BIN}/python"
PLOTS_PYTHON="/autofs/bal36/dten/miniforge3/envs/plots/bin/python"

CHROM_LENGTHS="${DATA_DIR}/chrom_lengths_${REFERENCE_NAME}.tsv"

TRUTH_PHASED="${TRUTH_PHASED:-${DERIVED_DIR}/truth${_ref_sfx}.het_snps.phased.vcf.gz}"
TRUTH_UNPHASED="${TRUTH_UNPHASED:-${DERIVED_DIR}/truth${_ref_sfx}.het_snps.unphased.vcf.gz}"
# experiment A (bin/make_truth_variants.sh): the same truth het SNPs plus
#   truthindel  biallelic het indels
#   truthall    every other het variant too: 1/2 sites split into two records,
#               MNPs and complex variants atomized
# phased (to score the added variants) and unphased (the phasers' input)
TRUTH_INDEL_PHASED="${DERIVED_DIR}/truth${_ref_sfx}.indel.phased.vcf.gz"
TRUTH_INDEL_UNPHASED="${DERIVED_DIR}/truth${_ref_sfx}.indel.unphased.vcf.gz"
TRUTH_ALL_PHASED="${DERIVED_DIR}/truth${_ref_sfx}.all.phased.vcf.gz"
TRUTH_ALL_UNPHASED="${DERIVED_DIR}/truth${_ref_sfx}.all.unphased.vcf.gz"

MINIMAP2="/autofs/bal36/dten/bin/minimap2-2.31_x64-linux/minimap2"

T2T_NAME="T2T"
T2T_SRC_DIR="${DATA_DIR}/t2t_source"
# the minimap2 index of CHM13 is the same for every sample
MMI_DIR="${PROJECT_ROOT}/data/t2t_index"

ALIGN_THREADS="${ALIGN_THREADS:-32}"
SORT_THREADS="${SORT_THREADS:-8}"
SORT_MEM="${SORT_MEM:-4G}"

ASM_DIR="${SAMPLE_ROOT}/asm"
HIFIASM="${SOFTWARE_DIR}/hifiasm/hifiasm"
ASM_THREADS="${ASM_THREADS:-48}"

if [[ "${SAMPLE}" == HG002 ]]; then
  ASM_PLATFORMS=(${ASM_PLATFORMS:-hifi revio ont4k ont5k})
  ASM_DEPTHS=(${ASM_DEPTHS:-20 30 50})
else
  ASM_PLATFORMS=(${ASM_PLATFORMS:-${PLATFORMS[*]}})
  ASM_DEPTHS=(${ASM_DEPTHS:-30})
fi

ASM_REGION="${ASM_REGION:-}"

# An assembly does not depend on the reference, so each one is built once, from
# the reads of the ASM_READS_REF cell, and dipcall then places it on either
# reference.
ASM_READS_REF="${ASM_READS_REF:-GRCh38}"

# hifiasm      reads only
# hifiasmtrio  hifiasm phased with the parental k-mers
# verkkotrio   verkko phased with the parental k-mers, PacBio cells only
ASM_CALLERS=(${ASM_CALLERS:-hifiasm hifiasmtrio verkkotrio})

# Parents' reads for the trio assemblies: BAM, CRAM or FASTQ files separated by
# spaces - Illumina R1 and R2, or long accurate reads (HiFi). Only reads are
# needed, no VCF and no parental assembly. Default: pat_reads and mat_reads of
# the sample in samples.tsv.
PAT_READS="${PAT_READS:-$(sample_path "${SAMPLE}" pat_reads)}"
MAT_READS="${MAT_READS:-$(sample_path "${SAMPLE}" mat_reads)}"
TRIO_DIR="${DATA_DIR}/trio"
TRIO_BIN="/autofs/bal36/dten/miniforge3/envs/asmtrio/bin"   # yak, meryl, merqury
TRIO_THREADS="${TRIO_THREADS:-32}"
TRIO_MEMORY="${TRIO_MEMORY:-200}"                             # GB, for meryl

VERKKO_BIN="/autofs/bal36/dten/miniforge3/envs/verkko/bin"
VERKKO_THREADS="${VERKKO_THREADS:-64}"
VERKKO_MEMORY="${VERKKO_MEMORY:-400}"                         # GB

# dipcall runs without a PAR file: assemblies are compared on autosomes only
DIPCALL_BIN="/autofs/bal36/dten/miniforge3/envs/asmeval/bin"
DIPCALL_THREADS="${DIPCALL_THREADS:-16}"

# The diploid HG002 Q100 assembly, maternal and paternal chromosomes in one
# FASTA (chrN_MATERNAL, chrN_PATERNAL), for the parental-haplotype truth of
# bin/matpat_truth.sh
Q100_FASTA="${Q100_FASTA:-${PROJECT_ROOT}/data/q100/hg002v1.1.fasta.gz}"
