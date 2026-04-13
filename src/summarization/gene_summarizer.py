import pandas as pd
from typing import List
import os

def summarize_enhancers(annotation_df_path: str) -> pd.DataFrame:
    """
    Summarizes enhancers by grouping variants and mapped genes.
    """
    df = pd.read_csv(annotation_df_path, sep="\t")
    
    # Required columns: enhancer_id, gene_name, distance, variant_id, pos
    summary = (
        df.groupby(["enhancer_id", "gene_name", "distance"])
          .agg(
              n_variants=("variant_id", "count"),
              variants=("variant_id", lambda x: ",".join(x)),
              min_pos=("pos", "min"),
              max_pos=("pos", "max")
          )
          .reset_index()
    )
    return summary

def summarize_unique_genes(enhancer_summary_df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarizes unique genes with multiple enhancers and variants.
    """
    gene_summary = (
        enhancer_summary_df.groupby("gene_name")
          .agg(
              n_enhancers=("enhancer_id", "nunique"),
              n_variants=("n_variants", "sum"),
              min_distance=("distance", "min")
          )
          .reset_index()
    )
    return gene_summary

def rank_genes(gene_summary_df: pd.DataFrame) -> pd.DataFrame:
    """
    Ranks genes based on the number of enhancers, variants, and distance.
    """
    df = gene_summary_df.copy()
    
    # Priority score: Higher enhancers, higher variants, lower distance
    df["priority_score"] = (
        df["n_enhancers"] * 100 +
        df["n_variants"] * 10 +
        (1 / df["min_distance"].replace(0, 1)) * 1000
    )

    df = df.sort_values(
        ["n_enhancers", "n_variants", "min_distance"], 
        ascending=[False, False, True]
    )
    return df

def filter_curated(enhancer_summary_df: pd.DataFrame, max_dist: int = 50000) -> pd.DataFrame:
    """
    Filters for curated genes and specific distance threshold.
    """
    df = enhancer_summary_df.copy()
    
    # Identify predicted/unnamed genes
    df["is_predicted"] = df["gene_name"].str.match(r"^(Gm\d+|.*Rik$|LOC)", na=False)
    
    filtered = df[(df["distance"] <= max_dist) & (~df["is_predicted"])]
    return filtered

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Summarize and rank candidate genes.")
    parser.add_argument("--annotated", required=True, help="Annotated variants TSV.")
    parser.add_argument("--out-enhancer-summary", required=True, help="Output enhancer summary TSV.")
    parser.add_argument("--out-gene-summary", required=True, help="Output gene summary TSV.")
    parser.add_argument("--out-ranked", required=True, help="Output ranked genes TSV.")
    parser.add_argument("--out-curated", required=True, help="Output curated genes TSV.")
    args = parser.parse_args()

    enh_summary = summarize_enhancers(args.annotated)
    enh_summary.to_csv(args.out_enhancer_summary, sep="\t", index=False)
    
    curated = filter_curated(enh_summary)
    curated.to_csv(args.out_curated, sep="\t", index=False)
    
    gene_summary = summarize_unique_genes(curated)
    gene_summary.to_csv(args.out_gene_summary, sep="\t", index=False)
    
    ranked = rank_genes(gene_summary)
    ranked.to_csv(args.out_ranked, sep="\t", index=False)
    
    print(f"Ranked {len(ranked)} candidate genes.")
