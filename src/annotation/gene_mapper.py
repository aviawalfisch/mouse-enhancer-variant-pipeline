import pandas as pd
from typing import Optional
import os

def map_and_score(
    nearest_bed_output: str, 
    activity_file: Optional[str] = None
) -> pd.DataFrame:
    """
    Scores enhancer-gene pairs based on distance and optionally activity.

    Args:
        nearest_bed_output: Output from 'bedtools closest'.
        activity_file: Optional path to activity TSV with mean_signal.

    Returns:
        A DataFrame with scored enhancer-gene pairs.
    """
    nearest = pd.read_csv(nearest_bed_output, sep="\t", header=None)
    nearest.columns = [
        "enh_chr", "enh_start", "enh_end", "enhancer_id",
        "gene_chr", "gene_start", "gene_end", "gene_name", "distance"
    ]

    if activity_file:
        activity = pd.read_csv(activity_file, sep="\t")
        df = nearest.merge(
            activity,
            on=["enhancer_id"],
            suffixes=("", "_act"),
            how="left"
        )
        # Use mean_signal for scoring
        df["distance_adj"] = df["distance"].replace(0, 1)
        df["score"] = df["mean_signal"] / df["distance_adj"]
    else:
        # Distance-only scoring
        df = nearest.copy()
        df["distance_adj"] = df["distance"].replace(0, 1)
        df["score"] = 1 / df["distance_adj"]

    return df

def attach_to_variants(
    variants_df_path: str, 
    scored_enhancers_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Attaches gene annotations back to the variant list.

    Args:
        variants_df_path: Path to the parsed variants TSV.
        scored_enhancers_df: Scored enhancer-gene DataFrame.

    Returns:
        A merged DataFrame.
    """
    variants = pd.read_csv(variants_df_path, sep="\t")
    
    # Merge on enhancer coordinates to ensure correct mapping
    final = variants.merge(
        scored_enhancers_df,
        on=["enh_chr", "enh_start", "enh_end"],
        how="left"
    )
    
    return final

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Map and score enhancer-gene pairs.")
    parser.add_argument("--nearest", required=True, help="Bedtools closest output.")
    parser.add_argument("--activity", help="Optional activity TSV.")
    parser.add_argument("--variants", required=True, help="Parsed variants TSV.")
    parser.add_argument("--out", required=True, help="Output final annotation TSV.")
    args = parser.parse_args()

    scored_df = map_and_score(args.nearest, args.activity)
    final_df = attach_to_variants(args.variants, scored_df)
    final_df.to_csv(args.out, sep="\t", index=False)
    
    print(f"Annotated {len(final_df)} variants with genes.")
