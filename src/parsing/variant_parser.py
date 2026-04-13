import pandas as pd
from typing import List, Tuple
import os

def parse_enhancer_hits(
    input_file: str, 
    n_samples: int = 6
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parses the enhancer hits file and extracts variant and enhancer information.

    Args:
        input_file: Path to the raw enhancer hits file.
        n_samples: Number of sample columns in the input file.

    Returns:
        A tuple of (variants_df, unique_enhancers_df).
    """
    # Fixed VCF-like columns before sample columns
    vcf_base = [
        "chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"
    ]

    # Sample columns
    sample_cols = [f"sample_{i+1}" for i in range(n_samples)]

    # Annotation columns at the end
    annot_cols = ["enh_chr", "enh_start", "enh_end", "annotation"]

    all_cols = vcf_base + sample_cols + annot_cols

    df = pd.read_csv(
        input_file,
        sep="\t",
        header=None,
        names=all_cols,
        dtype=str
    )

    # Basic cleanup
    df["chr"] = df["chr"].apply(lambda x: x if str(x).startswith("chr") else f"chr{x}")
    df["pos"] = df["pos"].astype(int)
    df["enh_chr"] = df["enh_chr"].apply(lambda x: x if str(x).startswith("chr") else f"chr{x}")
    df["enh_start"] = df["enh_start"].astype(int)
    df["enh_end"] = df["enh_end"].astype(int)

    # Generate Variant ID
    df["variant_id"] = (
        df["chr"].str.replace("^chr", "", regex=True) + ":" +
        df["pos"].astype(str) + ":" +
        df["ref"] + ">" + df["alt"]
    )

    # Keep a clean output table for variants
    variants_df = df[[
        "variant_id", "chr", "pos", "ref", "alt",
        "qual", "filter", "info",
        "enh_chr", "enh_start", "enh_end", "annotation"
    ]].copy()

    # Create unique enhancer BED
    enhancers_df = (
        df[["enh_chr", "enh_start", "enh_end", "annotation"]]
        .drop_duplicates()
        .sort_values(["enh_chr", "enh_start", "enh_end"])
        .copy()
    )

    enhancers_df["enhancer_id"] = (
        enhancers_df["enh_chr"] + ":" +
        enhancers_df["enh_start"].astype(str) + "-" +
        enhancers_df["enh_end"].astype(str)
    )

    return variants_df, enhancers_df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Parse enhancer hits into variants and enhancers.")
    parser.add_argument("--input", required=True, help="Input enhancer hits file.")
    parser.add_argument("--out-variants", required=True, help="Output variants TSV.")
    parser.add_argument("--out-enhancers", required=True, help="Output enhancers BED.")
    args = parser.parse_args()

    v_df, e_df = parse_enhancer_hits(args.input)
    v_df.to_csv(args.out_variants, sep="\t", index=False)
    
    # BED format: chr, start, end, name
    e_df[["enh_chr", "enh_start", "enh_end", "enhancer_id"]].to_csv(
        args.out_enhancers, sep="\t", index=False, header=False
    )
    
    print(f"Processed {len(v_df)} variants and {len(e_df)} unique enhancers.")
