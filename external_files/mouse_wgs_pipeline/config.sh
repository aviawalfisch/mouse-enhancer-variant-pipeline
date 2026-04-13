#!/bin/bash
# ==============================================================================
# Mouse WGS Germline Variant Calling Pipeline - Configuration
# ==============================================================================
# This file contains all configurable parameters for the pipeline.
# Source this file from all step scripts.
# ==============================================================================

# ==============================================================================
# USER MUST EDIT - REQUIRED PARAMETERS
# ==============================================================================
# TODO: Set these paths before running the pipeline

# FASTQ input directory (contains sample subdirectories like MSDS01/, MSDS02/)
FASTQ_DIR=""  # TODO: e.g., "/data/raw/mouse_wgs"

# Project output directory (all results will be written here)
PROJECT_DIR=""  # TODO: e.g., "/scratch/wgs_project"

# Reference FASTA (GRCm38) - REQUIRED, must exist
REF_FA=""  # TODO: e.g., "/refs/GRCm38/GRCm38.fa"

# ==============================================================================
# USER MUST EDIT - RESOURCE PARAMETERS
# ==============================================================================
# TODO: Adjust based on your system/cluster resources

THREADS=4                    # CPU cores for alignment/calling
RAM="16G"                    # Max RAM per job (format: 16G or 16000M)
JAVA_HEAP="12G"              # Java heap for GATK (should be < RAM - 4G)
TEMP_DIR="/tmp"              # Fast scratch for sorting/GenomicsDB

# ==============================================================================
# USER MUST EDIT - ENVIRONMENT STRATEGY
# ==============================================================================
# TODO: Uncomment and modify ONE of the following blocks

# Option 1: Conda environment
# CONDA_ENV="bio"
# eval "$(conda shell.bash hook)"
# conda activate "$CONDA_ENV"

# Option 2: Module system
# module load bwa/0.7.17
# module load samtools/1.17
# module load gatk/4.4.0.0
# module load bcftools/1.17

# Option 3: Tools already in PATH (no action needed)

# ==============================================================================
# OPTIONAL PARAMETERS - DEFAULTS SHOULD WORK
# ==============================================================================

# Interval mode: all | autosomes | custom
INTERVAL_MODE="all"
# Custom intervals file (only used if INTERVAL_MODE=custom)
INTERVALS=""  # e.g., "/refs/targets.bed"

# BQSR settings (disabled by default - no mouse known-sites)
ENABLE_BQSR=false
KNOWN_SITES=""  # e.g., "/refs/mgp_known_snps.vcf.gz"

# Contig aliasing (fail-fast by default)
ALLOW_CONTIG_ALIAS=false

# Reference download (disabled by default)
DOWNLOAD_REF=false
REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"

# Joint calling method: genomicsdb | combine
JOINT_CALLER="genomicsdb"

# Annotation tool: snpeff | skip
ANNOTATOR="snpeff"

# ==============================================================================
# HARD FILTER THRESHOLDS (CONFIGURABLE)
# ==============================================================================
# SNP filters
SNP_QD_THRESHOLD="2.0"
SNP_FS_THRESHOLD="60.0"
SNP_MQ_THRESHOLD="40.0"
SNP_SOR_THRESHOLD="3.0"
SNP_MQRANKSUM_THRESHOLD="-12.5"
SNP_READPOSRANKSUM_THRESHOLD="-8.0"

# INDEL filters
INDEL_QD_THRESHOLD="2.0"
INDEL_FS_THRESHOLD="200.0"
INDEL_SOR_THRESHOLD="10.0"
INDEL_READPOSRANKSUM_THRESHOLD="-20.0"

# Genotype filters
GT_DP_THRESHOLD="10"
GT_GQ_THRESHOLD="10"

# ==============================================================================
# DERIVED PATHS (DO NOT EDIT)
# ==============================================================================
LOGS_DIR="${PROJECT_DIR}/logs"
ALIGNED_LANES_DIR="${PROJECT_DIR}/aligned_lanes"
PROCESSED_BAM_DIR="${PROJECT_DIR}/processed_bam"
GVCFS_DIR="${PROJECT_DIR}/gvcfs"
GENOMICSDB_DIR="${PROJECT_DIR}/genomicsdb"
JOINT_DIR="${PROJECT_DIR}/joint"
FILTERED_DIR="${PROJECT_DIR}/filtered"
ANNOTATED_DIR="${PROJECT_DIR}/annotated"
QC_DIR="${PROJECT_DIR}/qc"
MANIFEST="${PROJECT_DIR}/samples.tsv"
MASTER_LOG="${LOGS_DIR}/master.log"

# Reference derived paths
REF_FAI="${REF_FA}.fai"
REF_DICT="${REF_FA%.fa}.dict"
REF_BWA="${REF_FA}.bwt"

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

log_info() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1"
    echo "$msg"
    [[ -d "$LOGS_DIR" ]] && echo "$msg" >> "$MASTER_LOG"
}

log_warn() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $1"
    echo "$msg" >&2
    [[ -d "$LOGS_DIR" ]] && echo "$msg" >> "$MASTER_LOG"
}

log_error() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1"
    echo "$msg" >&2
    [[ -d "$LOGS_DIR" ]] && echo "$msg" >> "$MASTER_LOG"
    exit 1
}

check_required_var() {
    local var_name="$1"
    local var_value="${!var_name}"
    if [[ -z "$var_value" ]]; then
        log_error "Required variable $var_name is not set. Edit config.sh and set it."
    fi
}

check_file_exists() {
    local filepath="$1"
    local description="$2"
    if [[ ! -f "$filepath" ]]; then
        log_error "$description not found: $filepath"
    fi
}

check_dir_exists() {
    local dirpath="$1"
    local description="$2"
    if [[ ! -d "$dirpath" ]]; then
        log_error "$description not found: $dirpath"
    fi
}

check_command() {
    local cmd="$1"
    if ! command -v "$cmd" &> /dev/null; then
        log_error "Required command '$cmd' not found in PATH. Check your environment."
    fi
}

ensure_dir() {
    local dirpath="$1"
    mkdir -p "$dirpath" || log_error "Failed to create directory: $dirpath"
}

# Atomic write: write to temp file then move
atomic_mv() {
    local src="$1"
    local dst="$2"
    mv "$src" "$dst" || log_error "Failed to move $src to $dst"
}

# Check if output exists and is non-empty (for skip-if-exists)
output_exists() {
    local filepath="$1"
    [[ -f "$filepath" && -s "$filepath" ]]
}

# Get contigs from reference .fai
get_ref_contigs() {
    cut -f1 "$REF_FAI"
}

# Get autosomes (works with chr or no-chr naming)
get_autosomes() {
    get_ref_contigs | grep -E '^(chr)?[0-9]+$' | sort -V
}

# Validate contig compatibility between two files
validate_contigs() {
    local file1="$1"
    local file1_contigs="$2"
    local file2="$3"
    local file2_contigs="$4"
    
    local diff_result
    diff_result=$(diff <(echo "$file1_contigs" | sort) <(echo "$file2_contigs" | sort) 2>&1)
    
    if [[ -n "$diff_result" ]]; then
        if [[ "$ALLOW_CONTIG_ALIAS" == "true" ]]; then
            log_warn "Contig mismatch between $file1 and $file2 (aliasing enabled)"
            log_warn "Differences: $diff_result"
        else
            log_error "Contig mismatch between $file1 and $file2. Set ALLOW_CONTIG_ALIAS=true to override (not recommended)."
        fi
    fi
}

# Validate all required inputs before heavy compute
validate_config() {
    log_info "Validating configuration..."
    
    # Required variables
    check_required_var "FASTQ_DIR"
    check_required_var "PROJECT_DIR"
    check_required_var "REF_FA"
    
    # Required files/directories
    check_dir_exists "$FASTQ_DIR" "FASTQ directory"
    check_file_exists "$REF_FA" "Reference FASTA"
    
    # Check reference indices
    if [[ ! -f "$REF_FAI" ]]; then
        log_warn "Reference .fai index not found. Will create in setup step."
    fi
    if [[ ! -f "$REF_DICT" ]]; then
        log_warn "Reference .dict not found. Will create in setup step."
    fi
    if [[ ! -f "$REF_BWA" ]]; then
        log_warn "BWA index not found. Will create in setup step."
    fi
    
    # Check required tools
    check_command "bwa"
    check_command "samtools"
    check_command "gatk"
    check_command "bcftools"
    
    # Optional tools
    if [[ "$ANNOTATOR" == "snpeff" ]] && ! command -v snpEff &> /dev/null; then
        log_warn "snpEff not found. Annotation will be skipped."
        ANNOTATOR="skip"
    fi
    
    log_info "Configuration validated successfully."
}

# ==============================================================================
# END OF CONFIG
# ==============================================================================
