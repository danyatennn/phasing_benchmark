PROJECT_ROOT="/autofs/bal36/dten/phasing_benchmark"
DATA_DIR="${PROJECT_ROOT}/data"

DERIVED_DIR="${PROJECT_ROOT}/derived"
DEPTHS=(${PHASING_DEPTHS:-10 20 30 40 50})
PLATFORMS=(${PHASING_PLATFORMS:-ont4k ont5k hifi revio})
PHASING_BIN="/autofs/bal36/dten/miniforge3/envs/phasing/bin"
BCFTOOLS="${PHASING_BIN}/bcftools"
TABIX="${PHASING_BIN}/tabix"
BGZIP="${PHASING_BIN}/bgzip"
SAMPLE="HG002"

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

CLAIR3_DIR="${PROJECT_ROOT}/clair3"
LOGS_DIR="${PROJECT_ROOT}/logs"
VCF_DIR="${PROJECT_ROOT}/vcf"
PHASE_DIR="${PROJECT_ROOT}/phased"
METRICS_DIR="${PROJECT_ROOT}/metrics"
TMP_DIR="${PROJECT_ROOT}/tmp"

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

REGION_NAME="${REGION_NAME:-wgs}"

CLAIR3_MODEL_ont4k="${CLAIR3_MODEL_ont4k:-r1041_e82_400bps_sup_v410}"
CLAIR3_MODEL_ont5k="${CLAIR3_MODEL_ont5k:-r1041_e82_400bps_sup_v500}"
CLAIR3_MODEL_hifi="${CLAIR3_MODEL_hifi:-hifi_sequel2}"
CLAIR3_MODEL_revio="${CLAIR3_MODEL_revio:-hifi_revio}"

VARIANT_CALLERS=(clair3 deepvariant longcalld)
CALLERS=(${PHASING_CALLERS:-truth "${VARIANT_CALLERS[@]}"})
PHASERS=(${PHASING_PHASERS:-whatshap longphase margin hapcut2})

CALL_THREADS="${CALL_THREADS:-32}"

CLAIR3_JOBS="${CLAIR3_JOBS:-1}"
DEEPVARIANT_JOBS="${DEEPVARIANT_JOBS:-2}"
LONGCALLD_JOBS="${LONGCALLD_JOBS:-2}"

SOFTWARE_DIR="${PROJECT_ROOT}/software"

WHATSHAP="${PHASING_BIN}/whatshap"
LONGPHASE="${SOFTWARE_DIR}/longphase_linux-x64"
LONGCALLD="${SOFTWARE_DIR}/longcallD-v0.0.11_x64-linux/longcallD"
DEEPVARIANT_SIF="${SOFTWARE_DIR}/deepvariant_1.10.0.sif"
MARGIN="${SOFTWARE_DIR}/margin/build/margin"
MARGIN_PARAMS_DIR="${SOFTWARE_DIR}/margin/params/phase"
HIFIASM="${SOFTWARE_DIR}/hifiasm/hifiasm"
HAPCUT2_BIN="/autofs/bal36/dten/miniforge3/envs/hapcut2/bin"

WHATSHAP_DOWNSAMPLING="${WHATSHAP_DOWNSAMPLING:-15}"
MARGIN_DEPTH="${MARGIN_DEPTH:-30}"

PHASE_JOBS="${PHASE_JOBS:-6}"
PHASE_THREADS="${PHASE_THREADS:-8}"
EVAL_JOBS="${EVAL_JOBS:-6}"

HAPLOTAG_REGION="${HAPLOTAG_REGION:-}"

CHROM_LENGTHS="${DATA_DIR}/chrom_lengths_${REFERENCE_NAME}.tsv"

TRUTH_PHASED="${TRUTH_PHASED:-${DERIVED_DIR}/truth${_ref_sfx}.het_snps.phased.vcf.gz}"
TRUTH_UNPHASED="${TRUTH_UNPHASED:-${DERIVED_DIR}/truth${_ref_sfx}.het_snps.unphased.vcf.gz}"

MINIMAP2="/autofs/bal36/dten/bin/minimap2-2.31_x64-linux/minimap2"

T2T_NAME="T2T"
T2T_SRC_DIR="${DATA_DIR}/t2t_source"
MMI_DIR="${DATA_DIR}/t2t_index"

ALIGN_THREADS="${ALIGN_THREADS:-32}"
SORT_THREADS="${SORT_THREADS:-8}"
SORT_MEM="${SORT_MEM:-4G}"

ASM_DIR="${PROJECT_ROOT}/asm"
HIFIASM="${SOFTWARE_DIR}/hifiasm/hifiasm"
ASM_THREADS="${ASM_THREADS:-48}"

ASM_PLATFORMS=(${ASM_PLATFORMS:-hifi revio ont4k ont5k})
ASM_DEPTHS=(${ASM_DEPTHS:-20 30 50})

ASM_REGION="${ASM_REGION:-}"

DIPCALL_BIN="/autofs/bal36/dten/miniforge3/envs/asmeval/bin"
DIPCALL_PAR="${DIPCALL_PAR:-}"
DIPCALL_THREADS="${DIPCALL_THREADS:-16}"
