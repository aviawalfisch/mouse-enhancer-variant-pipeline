#!/bin/bash
# ==============================================================================
# 01_setup.sh - Validate reference, create indices, generate sample manifest
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ==============================================================================
# MAIN
# ==============================================================================
log_info "Step 01: Setup starting..."

# Validate required config
validate_config

# Create output directories
ensure_dir "$PROJECT_DIR"
ensure_dir "$LOGS_DIR"
ensure_dir "${LOGS_DIR}/setup"
ensure_dir "$ALIGNED_LANES_DIR"
ensure_dir "$PROCESSED_BAM_DIR"
ensure_dir "$GVCFS_DIR"
ensure_dir "$JOINT_DIR"
ensure_dir "$FILTERED_DIR"
ensure_dir "$ANNOTATED_DIR"
ensure_dir "$QC_DIR"

# ==============================================================================
# REFERENCE INDICES
# ==============================================================================
log_info "Checking reference indices..."

# Download reference if enabled and missing
if [[ "$DOWNLOAD_REF" == "true" ]] && [[ ! -f "$REF_FA" ]]; then
    log_info "Downloading reference from $REF_URL..."
    wget -O "${REF_FA}.gz" "$REF_URL"
    gunzip "${REF_FA}.gz"
fi

check_file_exists "$REF_FA" "Reference FASTA"

# Create .fai index
if [[ ! -f "$REF_FAI" ]]; then
    log_info "Creating reference .fai index..."
    samtools faidx "$REF_FA"
fi

# Create .dict
if [[ ! -f "$REF_DICT" ]]; then
    log_info "Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R "$REF_FA" -O "$REF_DICT"
fi

# Create BWA index
if [[ ! -f "$REF_BWA" ]]; then
    log_info "Creating BWA index (this may take a while)..."
    bwa index "$REF_FA"
fi

log_info "Reference indices verified."

# ==============================================================================
# SAMPLE DISCOVERY
# ==============================================================================
log_info "Discovering samples from $FASTQ_DIR..."

# Manifest header
echo -e "sample_id\tlane_id\tread1_path\tread2_path" > "$MANIFEST"

sample_count=0
lane_count=0

# Find all R1 files and derive R2
while IFS= read -r -d '' r1_file; do
    # Derive R2 path
    r2_file="${r1_file/_1.fq.gz/_2.fq.gz}"
    
    if [[ ! -f "$r2_file" ]]; then
        log_error "Missing R2 for R1: $r1_file (expected: $r2_file)"
    fi
    
    # Check files are non-empty
    if [[ ! -s "$r1_file" ]]; then
        log_error "Empty R1 file: $r1_file"
    fi
    if [[ ! -s "$r2_file" ]]; then
        log_error "Empty R2 file: $r2_file"
    fi
    
    # Extract sample ID from directory name
    sample_dir=$(dirname "$r1_file")
    sample_id=$(basename "$sample_dir")
    
    # Extract lane ID from filename
    filename=$(basename "$r1_file")
    
    # Try to extract lane from L## pattern
    if [[ "$filename" =~ _L([0-9]+)_ ]]; then
        lane_id="L${BASH_REMATCH[1]}"
    else
        # Fallback: use chunk ID from filename prefix
        prefix="${filename%_1.fq.gz}"
        # Create deterministic chunk ID
        lane_id="chunk$(printf '%03d' $((lane_count + 1)))"
    fi
    
    # Write to manifest
    echo -e "${sample_id}\t${lane_id}\t${r1_file}\t${r2_file}" >> "$MANIFEST"
    
    ((lane_count++)) || true
    
done < <(find "$FASTQ_DIR" -name "*_1.fq.gz" -print0 | sort -z)

# Count unique samples
sample_count=$(tail -n +2 "$MANIFEST" | cut -f1 | sort -u | wc -l)

if [[ $lane_count -eq 0 ]]; then
    log_error "No FASTQ pairs found in $FASTQ_DIR. Expected *_1.fq.gz and *_2.fq.gz files."
fi

log_info "Found $sample_count samples with $lane_count lane/chunk pairs."
log_info "Manifest written to: $MANIFEST"

# ==============================================================================
# VALIDATE READ GROUP UNIQUENESS
# ==============================================================================
log_info "Validating read group uniqueness..."

rg_ids=$(tail -n +2 "$MANIFEST" | awk -F'\t' '{print $1"."$2}')
rg_unique=$(echo "$rg_ids" | sort -u | wc -l)
rg_total=$(echo "$rg_ids" | wc -l)

if [[ "$rg_unique" -ne "$rg_total" ]]; then
    log_error "Duplicate read group IDs detected! Check sample/lane combinations."
fi

log_info "Read group uniqueness validated: $rg_total unique RGs."

# ==============================================================================
# GENERATE INTERVAL LIST
# ==============================================================================
log_info "Generating interval list (mode: $INTERVAL_MODE)..."

INTERVALS_FILE="${PROJECT_DIR}/intervals.list"

case "$INTERVAL_MODE" in
    all)
        get_ref_contigs > "$INTERVALS_FILE"
        ;;
    autosomes)
        get_autosomes > "$INTERVALS_FILE"
        ;;
    custom)
        if [[ -z "$INTERVALS" ]] || [[ ! -f "$INTERVALS" ]]; then
            log_error "INTERVAL_MODE=custom but INTERVALS file not found: $INTERVALS"
        fi
        cp "$INTERVALS" "$INTERVALS_FILE"
        ;;
    *)
        log_error "Unknown INTERVAL_MODE: $INTERVAL_MODE (expected: all, autosomes, custom)"
        ;;
esac

interval_count=$(wc -l < "$INTERVALS_FILE")
log_info "Interval list created with $interval_count intervals: $INTERVALS_FILE"

# ==============================================================================
# SUMMARY
# ==============================================================================
log_info "Setup complete!"
log_info "  Samples: $sample_count"
log_info "  Lanes/chunks: $lane_count"
log_info "  Intervals: $interval_count"
log_info "  Manifest: $MANIFEST"
log_info "  Reference: $REF_FA"
