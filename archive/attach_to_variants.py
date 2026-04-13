import pandas as pd

variants = pd.read_csv("parsed_variants_with_enhancers.tsv", sep="\t")
enh = pd.read_csv("enhancer_gene_scored.tsv", sep="\t")

final = variants.merge(
    enh,
    left_on=["enh_chr", "enh_start", "enh_end"],
    right_on=["enh_chr", "enh_start", "enh_end"],
    how="left"
)

final.to_csv("final_variant_annotation.tsv", sep="\t", index=False)
print(f"Saved {len(final)} rows to final_variant_annotation.tsv")