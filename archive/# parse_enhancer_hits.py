# parse_enhancer_hits.py
# Comments in English

import pandas as pd
from pathlib import Path
import os

# --------------------------------
# Paths (robust to working dir)
# --------------------------------
BASE_DIR = Path(__file__).resolve().parent

INPUT_FILE = BASE_DIR / "DNV-06-CL3-35.enhancer_hits.txt"
OUTPUT_VARIANTS = BASE_DIR / "parsed_variants_with_enhancers.tsv"
OUTPUT_ENHANCERS = BASE_DIR / "unique_enhancers.bed"

print("Running script...")
print("cwd:", os.getcwd())
print("script:", __file__)
print("input path:", INPUT_FILE)
print("input exists:", INPUT_FILE.exists())

if not INPUT_FILE.exists():
    raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

# --------------------------------
# Config
# --------------------------------
N_SAMPLES = 6

VCF_BASE = [
    "chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"
]

sample_cols = [f"sample_{i+1}" for i in range(N_SAMPLES)]

ANNOT_COLS = ["enh_chr", "enh_start", "enh_end", "annotation"]

all_cols = VCF_BASE + sample_cols + ANNOT_COLS

# --------------------------------
# Load data
# --------------------------------
df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    header=None,
    names=all_cols,
    dtype=str
)

# --------------------------------
# Cleanup types
# --------------------------------
df["chr"] = df["chr"].astype(str)
df["pos"] = df["pos"].astype(int)
df["enh_chr"] = df["enh_chr"].astype(str)
df["enh_start"] = df["enh_start"].astype(int)
df["enh_end"] = df["enh_end"].astype(int)

# --------------------------------
# Variant ID
# --------------------------------
df["variant_id"] = (
    "chr" + df["chr"].str.replace("^chr", "", regex=True) + ":" +
    df["pos"].astype(str) + ":" +
    df["ref"] + ">" + df["alt"]
)

# --------------------------------
# Normalize chromosome names
# --------------------------------
df["chr"] = df["chr"].apply(lambda x: x if str(x).startswith("chr") else f"chr{x}")
df["enh_chr"] = df["enh_chr"].apply(lambda x: x if str(x).startswith("chr") else f"chr{x}")

# --------------------------------
# Output: variants
# --------------------------------
out = df[[
    "variant_id", "chr", "pos", "ref", "alt",
    "qual", "filter", "info",
    "enh_chr", "enh_start", "enh_end", "annotation"
]].copy()

out.to_csv(OUTPUT_VARIANTS, sep="\t", index=False)

# --------------------------------
# Output: unique enhancers BED
# --------------------------------
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

# --------------------------------
# Summary
# --------------------------------
print(f"Saved parsed variants to: {OUTPUT_VARIANTS}")
print(f"Saved unique enhancers BED to: {OUTPUT_ENHANCERS}")
print(f"Total variants: {len(out)}")
print(f"Unique enhancers: {len(bed)}")