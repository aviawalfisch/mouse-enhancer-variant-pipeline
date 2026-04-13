# make_tss_bed.py
# Comments in English

import pandas as pd
from pathlib import Path
import os
import re
import gzip

# --------------------------------
# Paths (robust to working dir)
# --------------------------------
BASE_DIR = Path(__file__).resolve().parent

GTF_FILE = BASE_DIR / "gencode.vM32.annotation.gtf.gz"
OUTPUT_BED = BASE_DIR / "mm39_tss.bed"

print("Running script...")
print("cwd:", os.getcwd())
print("script:", __file__)
print("gtf path:", GTF_FILE)
print("gtf exists:", GTF_FILE.exists())
print("gtf is file:", GTF_FILE.is_file())

if not GTF_FILE.exists():
    raise FileNotFoundError(f"GTF file not found: {GTF_FILE}")

if not GTF_FILE.is_file():
    raise IsADirectoryError(f"Expected a file but found a directory: {GTF_FILE}")

# --------------------------------
# Helper to extract GTF attributes
# --------------------------------
def extract_gtf_attribute(attr_string: str, key: str):
    """
    Extract a value from a GTF attributes field.
    Example: key='gene_name' from 'gene_name "Sox9"; gene_id "..."'
    """
    pattern = rf'{re.escape(key)}\s+"([^"]+)"'
    match = re.search(pattern, attr_string)
    return match.group(1) if match else None

# --------------------------------
# Open plain or gzipped GTF
# --------------------------------
def open_gtf(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")

# --------------------------------
# Parse GTF and collect TSS entries
# --------------------------------
rows = []

with open_gtf(GTF_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")

        if len(parts) < 9:
            continue

        feature_type = parts[2]
        if feature_type != "gene":
            continue

        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        info = parts[8]

        gene_name = extract_gtf_attribute(info, "gene_name")
        if gene_name is None:
            continue

        # BED is 0-based, half-open
        if strand == "+":
            tss_start = start - 1
        elif strand == "-":
            tss_start = end - 1
        else:
            continue

        tss_end = tss_start + 1

        # Normalize chromosome names
        chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"

        rows.append([chrom, tss_start, tss_end, gene_name, strand])

# --------------------------------
# Build output table
# --------------------------------
df = pd.DataFrame(rows, columns=["chr", "start", "end", "gene", "strand"])

df = (
    df.drop_duplicates()
      .sort_values(["chr", "start", "end", "gene"])
      .reset_index(drop=True)
)

# Standard 4-column BED
bed = df[["chr", "start", "end", "gene"]].copy()

# --------------------------------
# Save BED
# --------------------------------
bed.to_csv(OUTPUT_BED, sep="\t", header=False, index=False)

# --------------------------------
# Summary
# --------------------------------
print(f"Saved TSS BED to: {OUTPUT_BED}")
print(f"Total TSS entries: {len(bed)}")
print(f"Unique genes: {bed['gene'].nunique()}")