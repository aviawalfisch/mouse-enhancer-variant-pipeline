#!/usr/bin/env python3
"""
scripts/run_variant_pipeline.py

Main entry point for the Mouse mm39 Noncoding Variant Analysis Pipeline.
Orchestrates the conversion of variant-enhancer hits to ranked candidate gene lists.
"""

import argparse
import os
import subprocess
import sys
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Ensure project root is in path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

def run_cmd(cmd: str):
    """Executes a shell command and checks for errors."""
    logger.info(f"Executing: {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}")
        sys.exit(e.returncode)

def main():
    parser = argparse.ArgumentParser(
        description="Run the Mouse mm39 Noncoding Variant Analysis Pipeline (Variant -> Enhancer -> Gene Mapping)."
    )
    # Required inputs
    parser.add_argument("--input", required=True, help="Path to raw variant-enhancer hits file (e.g., DNV-06-CL3-35.enhancer_hits.txt).")
    parser.add_argument("--gtf", required=True, help="Path to GENCODE GTF file (e.g., gencode.vM32.annotation.gtf.gz).")
    
    # Optional inputs/outputs
    parser.add_argument("--bigwig", help="Optional path to BigWig file (e.g., atac_mm39.bw) to calculate enhancer activity.")
    parser.add_argument("--outdir", default="results", help="Base directory for final results (default: results).")
    parser.add_argument("--interim", default="data/interim", help="Directory for intermediate files (default: data/interim).")
    parser.add_argument("--refdir", default="data/reference", help="Directory for reference processing (default: data/reference).")

    args = parser.parse_args()

    # Create directories
    os.makedirs(args.interim, exist_ok=True)
    os.makedirs(args.refdir, exist_ok=True)
    os.makedirs(args.outdir, exist_ok=True)

    # Derived paths
    variants_tsv = os.path.join(args.interim, "parsed_variants_with_enhancers.tsv")
    enhancers_bed = os.path.join(args.interim, "unique_enhancers.bed")
    tss_bed = os.path.join(args.refdir, "mm39_tss.bed")
    nearest_gene_tsv = os.path.join(args.interim, "enhancer_nearest_gene.tsv")
    activity_tsv = os.path.join(args.interim, "enhancer_activity.tsv")
    annotated_variants = os.path.join(args.interim, "annotated_variants.tsv")
    enhancer_summary = os.path.join(args.interim, "enhancer_summary.tsv")
    
    # Final Output Paths
    final_ranked = os.path.join(args.outdir, "unique_candidate_genes_ranked.tsv")
    final_curated = os.path.join(args.outdir, "enhancer_candidate_genes_curated.tsv")
    final_gene_summary = os.path.join(args.outdir, "unique_candidate_genes.tsv")

    # Pipeline Steps
    logger.info("Starting Mouse Variant Analysis Pipeline...")

    # Step 1: Parsing
    logger.info("Step 1: Parsing raw enhancer hits...")
    run_cmd(f"python src/parsing/variant_parser.py --input {args.input} --out-variants {variants_tsv} --out-enhancers {enhancers_bed}")

    # Step 2: TSS Generation
    logger.info("Step 2: Extracting TSS from GTF...")
    run_cmd(f"python src/annotation/tss_generator.py --gtf {args.gtf} --out {tss_bed}")

    # Step 3: Mapping (Requires bedtools)
    logger.info("Step 3: Mapping enhancers to nearest genes (requires bedtools)...")
    run_cmd(f"bedtools closest -a {enhancers_bed} -b {tss_bed} -d > {nearest_gene_tsv}")

    # Step 4: Activity Calculation (Conditional)
    activity_arg = ""
    if args.bigwig and os.path.exists(args.bigwig) and os.path.getsize(args.bigwig) > 0:
        logger.info("Step 4: Calculating enhancer activity from BigWig...")
        run_cmd(f"python src/utils/bigwig_utils.py --bw {args.bigwig} --bed {enhancers_bed} --out {activity_tsv}")
        activity_arg = f"--activity {activity_tsv}"
    else:
        logger.warning("Step 4: Skipping activity calculation (BigWig missing or empty).")

    # Step 5: Annotation
    logger.info("Step 5: Attaching genes and scoring enhancers...")
    run_cmd(f"python src/annotation/gene_mapper.py --nearest {nearest_gene_tsv} {activity_arg} --variants {variants_tsv} --out {annotated_variants}")

    # Step 6: Summarization and Ranking
    logger.info("Step 6: Ranking candidate genes and generating final tables...")
    run_cmd(f"python src/summarization/gene_summarizer.py --annotated {annotated_variants} --out-enhancer-summary {enhancer_summary} --out-gene-summary {final_gene_summary} --out-ranked {final_ranked} --out-curated {final_curated}")

    logger.info(f"Pipeline complete! Final tables available in: {args.outdir}/tables/")

if __name__ == "__main__":
    main()
