import pandas as pd
import os

INPUT_FILE = "DNV-06-CL3-35.enhancer_hits.txt"
OUTPUT_VARIANTS = "parsed_variants_with_enhancers.tsv"
OUTPUT_ENHANCERS = "unique_enhancers.bed"
#raise Exception("TEST")

print("Running script...")
print("cwd:", os.getcwd())
print("script:", __file__)
print("input exists:", os.path.exists(INPUT_FILE))

INPUT_FILE = "DNV-06-CL3-35.enhancer_hits.txt"
OUTPUT_VARIANTS = "parsed_variants_with_enhancers.tsv"
OUTPUT_ENHANCERS = "unique_enhancers.bed"

# Number of samples in the file
N_SAMPLES = 6

# Fixed VCF-like columns before sample columns
VCF_BASE = [
    "chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"
]

# Sample columns
sample_cols = [f"sample_{i+1}" for i in range(N_SAMPLES)]

# Annotation columns at the end
ANNOT_COLS = ["enh_chr", "enh_start", "enh_end", "annotation"]

all_cols = VCF_BASE + sample_cols + ANNOT_COLS

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    header=None,
    names=all_cols,
    dtype=str
)

# Basic cleanup
df["chr"] = df["chr"].astype(str)
df["pos"] = df["pos"].astype(int)
df["enh_chr"] = df["enh_chr"].astype(str)
df["enh_start"] = df["enh_start"].astype(int)
df["enh_end"] = df["enh_end"].astype(int)

# Variant ID
df["variant_id"] = (
    "chr" + df["chr"].str.replace("^chr", "", regex=True) + ":" +
    df["pos"].astype(str) + ":" +
    df["ref"] + ">" + df["alt"]
)

# Normalize chromosome names
df["chr"] = df["chr"].apply(lambda x: x if str(x).startswith("chr") else f"chr{x}")
df["enh_chr"] = df["enh_chr"].apply(lambda x: x if str(x).startswith("chr") else f"chr{x}")

# Keep a clean output table
out = df[[
    "variant_id", "chr", "pos", "ref", "alt",
    "qual", "filter", "info",
    "enh_chr", "enh_start", "enh_end", "annotation"
]].copy()

out.to_csv(OUTPUT_VARIANTS, sep="\t", index=False)

# Create unique enhancer BED
enh = (
    df[["enh_chr", "enh_start", "enh_end", "annotation"]]
    .drop_duplicates()
    .sort_values(["enh_chr", "enh_start", "enh_end"])
    .copy()
)

enh["enhancer_id"] = (
    enh["enh_chr"] + ":" +
    enh["enh_start"].astype(str) + "-" +
    enh["enh_end"].astype(str)
)

bed = enh[["enh_chr", "enh_start", "enh_end", "enhancer_id"]]
bed.to_csv(OUTPUT_ENHANCERS, sep="\t", index=False, header=False)

print(f"Saved parsed variants to: {OUTPUT_VARIANTS}")
print(f"Saved unique enhancers BED to: {OUTPUT_ENHANCERS}")
print(f"Total variants: {len(out)}")
print(f"Unique enhancers: {len(bed)}")