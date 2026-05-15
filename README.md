# TCGA-BRCA Differential Expression & Pathway Analysis

Bulk RNA-seq analysis comparing breast cancer tumor vs. normal tissue.  
**80 samples** (40 tumor, 40 normal) from TCGA-BRCA · DESeq2 · fgsea Hallmark pathways

---

## Pipeline

```
TCGA download (barcodes filtered — ~300MB, not 5GB)
        ↓
QC (library size, PCA, sample distances)
        ↓
DESeq2 (raw counts, Wald test, tumor vs normal)
        ↓
fgsea (Hallmark pathways, ranked gene list)
        ↓
Visualization (volcano, heatmap, enrichment plots)
```

---

## Key Results

### Upregulated in Tumor

| Pathway | NES | padj |
|---|---|---|
| HALLMARK_E2F_TARGETS | +2.22 | 2.1e-19 |
| HALLMARK_G2M_CHECKPOINT | +2.17 | 1.5e-15 |
| HALLMARK_MITOTIC_SPINDLE | +1.94 | 2.5e-07 |
| HALLMARK_MYC_TARGETS_V1 | +1.90 | 9.6e-07 |
| HALLMARK_MTORC1_SIGNALING | +1.77 | 7.3e-05 |

### Downregulated in Tumor

| Pathway | NES | padj |
|---|---|---|
| HALLMARK_ADIPOGENESIS | −1.99 | 7.8e-12 |
| HALLMARK_MYOGENESIS | −1.67 | 1.5e-04 |
| HALLMARK_FATTY_ACID_METABOLISM | −1.54 | 1.2e-02 |

**Story:** Tumor cells are locked in proliferation mode (E2F, MYC, mTORC1 active) while losing normal breast tissue identity (adipogenesis, myogenesis suppressed) — consistent with known BRCA biology.

---

## Requirements

```r
BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "fgsea"
))

install.packages(c("msigdbr", "data.table", "ggplot2", "pheatmap"))
```

---

## Usage

```r
source("tcga_brca.R")
```

- `set.seed(42)` used throughout for reproducibility
- Raw counts only (STAR - Counts, unstranded) — not TPM/FPKM
- Reference level: Solid Tissue Normal
- Filter: genes with ≥ 10 counts in at least 1 sample
- fgsea ranking: `log2FC × −log10(pvalue + 1e−300)`

---

## Project Structure

```
tcga-brca-deseq2/
├── tcga_brca.R          # full pipeline
├── README.md
└── .gitignore           # excludes GDCdata/ and large files
```

---

## What I Learned

- Bulk RNA-seq pipeline end-to-end on real cancer data
- Why raw counts matter for DESeq2 (not normalized data)
- Confounding, overlap, rank deficiency in experimental design
- Wald test, dispersion estimation, multiple testing correction
- fgsea ranked gene list logic and NES interpretation
- Git + SSH setup for version control

---

*Simran · TCGA RNA-seq learning project · 2026*
