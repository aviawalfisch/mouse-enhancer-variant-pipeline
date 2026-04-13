#!/bin/bash
# ==============================================================================
# Mouse WGS Germline Variant Calling Pipeline - Main Dispatcher
# ==============================================================================
# Usage: ./mouse_wgs_pipeline.sh <command> [options]
# Commands: setup, align, dedup, call, joint, filter, annotate, evaluate, all
# Options: --force (re-run even if outputs exist)
#          --sample <id> (run only for specific sample)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ==============================================================================
# USAGE
# ==============================================================================
usage() {
    cat << EOF
Mouse WGS Germline Variant Calling Pipeline

Usage: $0 <command> [options]

Commands:
  setup      Validate reference, create indices, generate sample manifest
  align      Align FASTQs per lane/chunk (BWA-MEM)
  dedup      Merge lanes, mark duplicates, index, run QC
  call       Run HaplotypeCaller in GVCF mode per sample
  joint      Consolidate GVCFs and joint-genotype cohort
  filter     Apply hard filters to SNPs and INDELs
  annotate   Annotate variants with snpEff (if available)
  evaluate   Generate QC reports and summary statistics
  all        Run all steps in order

Options:
  --force         Re-run step even if outputs exist
  --sample <id>   Run only for specific sample (align, dedup, call)
  --help          Show this help message

Examples:
  $0 setup                    # First: validate and generate manifest
  $0 align                    # Align all samples
  $0 align --sample MSDS01    # Align only MSDS01
  $0 all                      # Run entire pipeline
  $0 call --force             # Re-run calling even if GVCFs exist

EOF
    exit 0
}

# ==============================================================================
# PARSE ARGUMENTS
# ==============================================================================
COMMAND=""
FORCE=false
SAMPLE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        setup|align|dedup|call|joint|filter|annotate|evaluate|all)
            COMMAND="$1"
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --help|-h)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

if [[ -z "$COMMAND" ]]; then
    echo "Error: No command specified."
    usage
fi

export FORCE
export SAMPLE

# ==============================================================================
# RUN STEP
# ==============================================================================
run_step() {
    local step_name="$1"
    local step_script="${SCRIPT_DIR}/${step_name}.sh"
    
    if [[ ! -f "$step_script" ]]; then
        log_error "Step script not found: $step_script"
    fi
    
    log_info "=========================================="
    log_info "Starting step: $step_name"
    log_info "=========================================="
    
    bash "$step_script"
    local exit_code=$?
    
    if [[ $exit_code -eq 0 ]]; then
        log_info "Step $step_name completed successfully."
    else
        log_error "Step $step_name failed with exit code $exit_code"
    fi
    
    return $exit_code
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
log_info "Pipeline started: command=$COMMAND force=$FORCE sample=$SAMPLE"

case "$COMMAND" in
    setup)
        run_step "01_setup"
        ;;
    align)
        run_step "02_align"
        ;;
    dedup)
        run_step "03_dedup"
        ;;
    call)
        run_step "04_call"
        ;;
    joint)
        run_step "05_joint"
        ;;
    filter)
        run_step "06_filter"
        ;;
    annotate)
        run_step "07_annotate"
        ;;
    evaluate)
        run_step "08_evaluate"
        ;;
    all)
        run_step "01_setup"
        run_step "02_align"
        run_step "03_dedup"
        run_step "04_call"
        run_step "05_joint"
        run_step "06_filter"
        run_step "07_annotate"
        run_step "08_evaluate"
        log_info "=========================================="
        log_info "All steps completed successfully!"
        log_info "=========================================="
        ;;
esac

log_info "Pipeline finished: command=$COMMAND"
