import pandas as pd

# Read nearest-gene output from bedtools closest
nearest = pd.read_csv("enhancer_nearest_gene.tsv", sep="\t", header=None)

# Adjust these column names if your file has a slightly different structure
nearest.columns = [
    "enh_chr", "enh_start", "enh_end", "enhancer_id",
    "gene_chr", "gene_start", "gene_end", "gene_name", "distance"
]

activity = pd.read_csv("enhancer_activity.tsv", sep="\t")

df = nearest.merge(
    activity,
    left_on=["enhancer_id", "enh_chr", "enh_start", "enh_end"],
    right_on=["enhancer_id", "chr", "start", "end"],
    how="left"
)

df["distance"] = df["distance"].replace(0, 1)
df["score"] = df["mean_signal"] / df["distance"]

out = df[[
    "enhancer_id", "enh_chr", "enh_start", "enh_end",
    "gene_name", "distance", "mean_signal", "score"
]].copy()

out.to_csv("enhancer_gene_scored.tsv", sep="\t", index=False)
print(f"Saved {len(out)} rows to enhancer_gene_scored.tsv")