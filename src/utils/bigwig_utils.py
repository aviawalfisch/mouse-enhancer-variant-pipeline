import pandas as pd
import pyBigWig
from typing import Optional
import os

def calculate_bigwig_mean(
    bw_path: str, 
    bed_path: str, 
    id_col_idx: int = 3
) -> pd.DataFrame:
    """
    Calculates the mean BigWig signal over BED regions.

    Args:
        bw_path: Path to the BigWig file.
        bed_path: Path to the BED file.
        id_col_idx: 0-based index of the column to use as the ID in the BED file.

    Returns:
        A DataFrame with [id, chr, start, end, mean_signal].
    """
    bw = pyBigWig.open(bw_path)
    rows = []

    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
            # Use specified ID or default to range
            if len(parts) > id_col_idx:
                name = parts[id_col_idx]
            else:
                name = f"{chrom}:{start}-{end}"

            try:
                mean_val = bw.stats(chrom, start, end, type="mean")[0]
                if mean_val is None:
                    mean_val = 0.0
            except RuntimeError:
                mean_val = 0.0

            rows.append([name, chrom, start, end, mean_val])

    bw.close()
    
    df = pd.DataFrame(rows, columns=["enhancer_id", "chr", "start", "end", "mean_signal"])
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calculate mean BigWig signal over BED.")
    parser.add_argument("--bw", required=True, help="Input BigWig file.")
    parser.add_argument("--bed", required=True, help="Input BED file.")
    parser.add_argument("--out", required=True, help="Output TSV file.")
    args = parser.parse_args()

    result_df = calculate_bigwig_mean(args.bw, args.bed)
    result_df.to_csv(args.out, sep="\t", index=False)
    
    print(f"Calculated mean signal for {len(result_df)} regions.")
