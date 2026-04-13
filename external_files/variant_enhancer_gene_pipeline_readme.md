# Variant → Enhancer → Gene Mapping Pipeline (Mouse, mm39)

## Overview

This repository contains a full, reproducible pipeline for mapping noncoding variants to candidate genes via enhancer regions in mouse (mm39).

The core logic:

variant → enhancer → nearest gene

This provides an initial regulatory interpretation of noncoding variants.

---

## Biological Background

Enhancers are noncoding DNA regions that regulate gene expression by binding transcription factors and interacting with promoters via chromatin looping. They can act over long genomic distances (kb–Mb scale).

Therefore, variants located in enhancers may influence gene expression rather than protein sequence.

---

## Input Data

### Primary file

DNV-06-CL3-35.enhancer_hits.txt

This file already contains:
- variant positions
- enhancer overlaps

Meaning the pipeline starts from:

variant → enhancer

---

## Reference Data

Genome build: mm39

Gene annotation: GENCODE M32

Download:

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz

---

## Dependencies

### System

- Python ≥ 3.10
- bedtools

### Python

pip install pandas

---

## Pipeline Steps

### 1. Parse enhancer hits

Script: parse_enhancer_hits.py

Output:
parsed_variants_with_enhancers.tsv

---

### 2. Extract unique enhancers

Output:
unique_enhancers.bed

---

### 3. Generate TSS BED file

Script: make_tss_bed.py

```python
import pandas as pd

gtf_file = "gencode.vM32.annotation.gtf"
rows = []

with open(gtf_file) as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split("\t")
        if parts[2] != "gene": continue

        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        info = parts[8]

        gene_name = None
        for field in info.split(";"):
            if "gene_name" in field:
                gene_name = field.split('"')[1]

        if gene_name is None: continue

        tss = start-1 if strand == "+" else end-1
        rows.append([chrom, tss, tss+1, gene_name])

pd.DataFrame(rows).to_csv("mm39_tss.bed", sep="\t", header=False, index=False)
```

---

### 4. Map enhancers to nearest genes

```bash
bedtools closest -a unique_enhancers.bed -b mm39_tss.bed -d > enhancer_nearest_gene.tsv
```

---

### 5. Attach genes to variants

Script: attach_to_variants.py

```python
import pandas as pd

variants = pd.read_csv("parsed_variants_with_enhancers.tsv", sep="\t")
enh = pd.read_csv("enhancer_nearest_gene.tsv", sep="\t", header=None)

enh.columns = ["enh_chr","enh_start","enh_end","enhancer_id","gene_chr","gene_start","gene_end","gene_name","distance"]

final = variants.merge(
    enh,
    left_on=["enh_chr","enh_start","enh_end"],
    right_on=["enh_chr","enh_start","enh_end"],
    how="left"
)

final.to_csv("final_variant_annotation.tsv", sep="\t", index=False)
```

---

### 6. Extract simplified mapping

```python
import pandas as pd

df = pd.read_csv("final_variant_annotation.tsv", sep="\t")
clean = df[["variant_id","chr","pos","ref","alt","enhancer_id","gene_name","distance"]]
clean.to_csv("variant_to_nearest_gene.tsv", sep="\t", index=False)
```

---

### 7. Summarize enhancers

```python
import pandas as pd

df = pd.read_csv("variant_to_nearest_gene.tsv", sep="\t")

summary = (
    df.groupby(["enhancer_id","gene_name","distance"])
      .agg(
          n_variants=("variant_id","count"),
          variants=("variant_id", lambda x: ",".join(x)),
          min_pos=("pos","min"),
          max_pos=("pos","max")
      )
      .reset_index()
)

summary.to_csv("enhancer_candidate_genes.tsv", sep="\t", index=False)
```

---

### 8. Filter genes

```python
import pandas as pd

df = pd.read_csv("enhancer_candidate_genes.tsv", sep="\t")

df["is_predicted"] = df["gene_name"].str.match(r"^(Gm\\d+|.*Rik$|LOC)", na=False)

filtered = df[df["distance"] <= 50000]
filtered.to_csv("enhancer_candidate_genes_filtered.tsv", sep="\t", index=False)

curated = filtered[~filtered["is_predicted"]]
curated.to_csv("enhancer_candidate_genes_curated.tsv", sep="\t", index=False)
```

---

### 9. Unique gene summary

```python
import pandas as pd

df = pd.read_csv("enhancer_candidate_genes_curated.tsv", sep="\t")

gene_summary = (
    df.groupby("gene_name")
      .agg(
          n_enhancers=("enhancer_id","nunique"),
          n_variants=("n_variants","sum"),
          min_distance=("distance","min")
      )
      .reset_index()
)

gene_summary.to_csv("unique_candidate_genes.tsv", sep="\t", index=False)
```

---

### 10. Rank genes

```python
import pandas as pd

df = pd.read_csv("unique_candidate_genes.tsv", sep="\t")

df["priority_score"] = (
    df["n_enhancers"]*100 +
    df["n_variants"]*10 +
    (1/df["min_distance"].replace(0,1))*1000
)

df = df.sort_values(["n_enhancers","n_variants","min_distance"], ascending=[False,False,True])

df.to_csv("unique_candidate_genes_ranked.tsv", sep="\t", index=False)
```

---

## Key Outputs

- variant_to_nearest_gene.tsv
- enhancer_candidate_genes.tsv
- unique_candidate_genes.tsv
- unique_candidate_genes_ranked.tsv

---

## Limitations

- nearest-gene assumption only
- no chromatin activity
- no 3D genome contacts

---

## Comparison to ABC Model

ABC model defines:

ABC score = enhancer activity × 3D contact

Paper:
https://www.nature.com/articles/s41588-019-0538-0

Differences:

| Feature | This pipeline | ABC |
|--------|--------------|-----|
| Enhancer activity | ❌ | ✔ |
| 3D contact | ❌ | ✔ |
| Cell-type specificity | ❌ | ✔ |
| Nearest gene | ✔ | ❌ |

---

## Next Steps

- integrate ATAC-seq
- run ABC pipeline
- perform GO enrichment

