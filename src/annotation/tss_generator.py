import pandas as pd
import gzip
from typing import List
import os

def generate_tss_bed(gtf_file: str) -> pd.DataFrame:
    """
    Extracts TSS (Transcription Start Sites) from a GTF file.

    Args:
        gtf_file: Path to the GENCODE GTF file (gzipped or plain).

    Returns:
        A DataFrame containing TSS bed information: [chrom, start, end, gene_name].
    """
    rows = []
    
    # Check if gzipped
    if gtf_file.endswith(".gz"):
        f = gzip.open(gtf_file, "rt")
    else:
        f = open(gtf_file, "r")
        
    for line in f:
        if line.startswith("#"): 
            continue
        parts = line.strip().split("\t")
        if parts[2] != "gene": 
            continue

        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        info = parts[8]

        gene_name = None
        for field in info.split(";"):
            field = field.strip()
            if field.startswith("gene_name"):
                gene_name = field.split('"')[1]

        if gene_name is None: 
            continue

        # For TSS, we use start-1 for + strand and end-1 for - strand
        tss = start-1 if strand == "+" else end-1
        rows.append([chrom, tss, tss+1, gene_name])
        
    f.close()
    
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_name"])
    # Sort for bedtools compatibility
    df = df.sort_values(["chrom", "start", "end"])
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate TSS BED from GENCODE GTF.")
    parser.add_argument("--gtf", required=True, help="Input GTF file (can be .gz).")
    parser.add_argument("--out", required=True, help="Output TSS BED file.")
    args = parser.parse_args()

    tss_df = generate_tss_bed(args.gtf)
    tss_df.to_csv(args.out, sep="\t", header=False, index=False)
    
    print(f"Extracted {len(tss_df)} TSS entries.")
