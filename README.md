# Mouse Noncoding Variant Analysis Pipeline (mm39)

## 🧬 Biological Context

This pipeline analyzes noncoding variants identified in **mouse tail tissue** that overlap with known or predicted enhancer regions.

### What are Enhancers?
Enhancers are noncoding regulatory DNA elements that control gene expression. Unlike promoters, which are at the transcription start site, enhancers can be located tens or hundreds of kilobases away. They interact with their target gene's promoter via **chromatin looping**.

### Why Enhancer Variants Matter?
Variants in enhancers (e.g., SNPs, Indels) do not change the protein sequence directly. Instead, they may:
1. Alter the binding affinity of transcription factors (TFs).
2. Disrupt the regulatory activity of the enhancer.
3. Lead to the dysregulation (overexpression or silencing) of target genes, contributing to phenotypic changes or disease.

### Pipeline Logic
This project implements a baseline mapping strategy:
**Variant → Enhancer → Nearest Gene (TSS)**

---

## 🔬 Comparison to the ABC Model

The **Activity-by-Contact (ABC) model** is the current state-of-the-art for predicting enhancer–gene interactions.

**ABC Score = Enhancer Activity × 3D Contact**

| Feature | This Pipeline | ABC Model |
| :--- | :--- | :--- |
| **Enhancer Activity** | (Optional) ATAC signal | H3K27ac ChIP-seq + ATAC-seq |
| **3D Proximity** | Nearest TSS (Linear distance) | Hi-C data or contact models |
| **Specificity** | General proximity | Cell-type specific |

**Reference:** Fulco et al. "Activity-by-contact model of enhancer–promoter interactions from thousands of CRISPR perturbations." *Nature Genetics* (2019). [DOI: 10.1038/s41588-019-0538-0](https://doi.org/10.1038/s41588-019-0538-0)

This project represents a **baseline approximation**. Moving to a full ABC implementation is a natural next step for higher-confidence mapping.

---

## 🧱 Repository Structure

```text
.
├── README.md               # Project documentation
├── data/
│   ├── interim/            # Intermediate pipeline files
│   ├── processed/          # Placeholder for final processed datasets
│   ├── raw/                # Original input files (DNV hits, BigWigs)
│   └── reference/          # Reference annotations (GTF, TSS BED)
├── environment.yml         # Conda environment definition
├── requirements.txt        # Python dependencies
├── results/
│   ├── reports/            # Placeholder for plots and summaries
│   └── tables/             # Final TSV outputs (Ranked lists)
├── scripts/
│   └── run_variant_pipeline.py  # Entry-point script
├── src/                    # Modular Python source code
│   ├── __init__.py
│   ├── annotation/         # TSS and Gene mapping logic
│   ├── parsing/            # Variant/Enhancer parsing logic
│   ├── summarization/      # Ranking and filtering logic
│   └── utils/              # BigWig and general utilities
└── tests/                  # Unit tests
    └── test_tss_generator.py
```

---

## 🧪 Provenance Tracking: Data Summary

The following metrics were computed from the primary analysis of `DNV-06-CL3-35.enhancer_hits.txt`:

- **Number of Variants:** 42
- **Number of Unique Enhancers:** 24
- **Number of Candidate Genes (Total):** 20
- **Number of Curated Candidate Genes:** 6

---

## 📊 Output Classification

### Final Outputs (`results/tables/`)
- `unique_candidate_genes_ranked.tsv`: **Primary result.** Genes ranked by enhancer density, variant count, and TSS proximity.
- `enhancer_candidate_genes_curated.tsv`: List of high-confidence enhancers and their nearest curated genes (dist < 50kb).

### Supporting Outputs (`data/interim/`)
- `variant_to_nearest_gene.tsv`: Direct mapping of each variant to its closest gene.
- `enhancer_summary.tsv`: Statistics grouped at the enhancer level.

### Intermediate Files (`data/interim/`)
- `enhancer_nearest_gene.tsv`: Raw output from `bedtools closest`.
- `parsed_variants_with_enhancers.tsv`: Cleaned version of the raw input hits.

---

## 📌 Interpretation: Key Findings

- **Enhancer Hotspots:** The gene **Mir101c** was identified as a potential hotspot, associated with **5 unique enhancers** and **7 variants**.
- **Proximity High-Confidence:** The gene **Atp5g2** was found to be the closest to an identified enhancer, with a distance of only **266 bp**.
- **Candidate Genes:** Other curated genes of interest include **Chd7**, **Myo16**, and **Tmem120b**, all associated with variants in potential regulatory elements.

---

## ⚠️ Important Limitations

1. **Nearest Gene Approximation:** Enhancers often skip the nearest gene to regulate a more distant one. This pipeline only captures the closest TSS.
2. **Activity Modeling:** While distance and (optional) signal are used, this is not a full regulatory activity model.
3. **Chromatin State:** This analysis does not currently incorporate cell-type-specific H3K27ac or full ATAC-seq peaks unless specifically provided.
4. **3D Genome:** Linear distance is used as a proxy for physical contact, which does not account for chromatin loops or TAD boundaries.

---

## 📂 External Data & Prerequisites

To keep the repository lightweight, large genomic data files are not included in the version control. Before running the pipeline, ensure the following files are available:

1.  **Reference Annotation:**
    -   Download the GENCODE M32 (mm39) GTF: [gencode.vM32.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz)
    -   Place it in: `data/reference/`

2.  **Chromatin Accessibility (Optional):**
    -   Provide a BigWig (`.bw`) file for enhancer activity scoring (e.g., ATAC-seq or DNase-seq).
    -   Place it in: `data/raw/`

3.  **Variant Hits:**
    -   The input file `DNV-06-CL3-35.enhancer_hits.txt` should be in `data/raw/`.

---

## 🚀 Reproducibility: How to Run

### Installation
```bash
# Using Conda
conda env create -f environment.yml
conda activate mouse-variant-pipeline

# Or using Pip
pip install -r requirements.txt
```

### Usage
The entire analysis can be reproduced using the single entry-point script:

```bash
python scripts/run_variant_pipeline.py \
  --input data/raw/DNV-06-CL3-35.enhancer_hits.txt \
  --gtf data/reference/gencode.vM32.annotation.gtf.gz \
  --bigwig data/raw/atac_mm39.bw \
  --outdir results/
```

### Testing
To run the automated test suite, use:
```bash
PYTHONPATH=. pytest tests/
```

---

## 🧭 Next Steps

- **Integrate ATAC-seq:** Use higher-resolution chromatin accessibility peaks to refine enhancer boundaries.
- **Run ABC Model:** Apply the full Activity-by-Contact logic (DOI: 10.1038/s41588-019-0538-0) for higher-confidence mapping.
- **Gene Ontology (GO) Enrichment:** Run enrichment analysis on the curated gene list to identify shared biological pathways.
- **Cell-type Specificity:** Incorporate enhancer databases (e.g., ENCODE, SCREEN) for specific mouse tissues.
