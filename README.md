# TCGA-BRCA Differential Expression & Pathway Analysis

RNA-seq analysis of breast cancer tumor vs. normal tissue using TCGA data, DESeq2, and fgsea.

## Overview

This pipeline downloads 40 primary tumor and 40 solid tissue normal samples from TCGA-BRCA, runs differential expression analysis, and performs gene set enrichment analysis using Hallmark pathways.

## Pipeline Steps

1. **Data download** — TCGAbiolinks, filtered to 80 samples via barcodes (~300MB vs full 5GB)
2. **QC** — Library size check, PCA, sample distance heatmap
3. **DESeq2** — Raw counts, Wald test, tumor vs normal
4. **fgsea** — Hallmark pathways, ranked by log2FC × −log10(p-value)

## Key Results (TCGA-BRCA, n=80)

### Upregulated in Tumor
| Pathway | NES | padj |
|---|---|---|
| HALLMARK_E2F_TARGETS | 2.22 | 2.1e-19 |
| HALLMARK_G2M_CHECKPOINT | 2.17 | 1.5e-15 |
| HALLMARK_MITOTIC_SPINDLE | 1.94 | 2.5e-07 |
| HALLMARK_MYC_TARGETS_V1 | 1.90 | 9.6e-07 |
| HALLMARK_MTORC1_SIGNALING | 1.77 | 7.3e-05 |

### Downregulated in Tumor
| Pathway | NES | padj |
|---|---|---|
| HALLMARK_ADIPOGENESIS | −1.99 | 7.8e-12 |
| HALLMARK_MYOGENESIS | −1.67 | 1.5e-04 |
| HALLMARK_FATTY_ACID_METABOLISM | −1.54 | 1.2e-02 |

**Biological interpretation:** Tumor samples show constitutive activation of proliferation programs (E2F, MYC, mTORC1) alongside loss of normal breast tissue identity (adipogenesis, myogenesis) — consistent with known BRCA hallmarks.

## Requirements

```r
BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "fgsea"
))
install.packages(c("msigdbr", "ggplot2", "pheatmap"))
```

## Usage

```r
source("tcga_brca.R")
# Downloads data, runs DESeq2 + fgsea, outputs results
# set.seed(42) used throughout for reproducibility
```

## Project Structure

```
tcga-brca-deseq2/
├── tcga_brca.R       # full pipeline
├── README.md
└── .gitignore        # excludes GDCdata/ and large files
```

## Notes

- Raw counts used (STAR - Counts, unstranded) — not TPM/FPKM
- Reference level: Solid Tissue Normal (log2FC = Tumor relative to Normal)
- Filtering: genes with ≥10 counts in at least 1 sample retained
- fgsea ranking metric: log2FC × −log10(pvalue + 1e−300)

## Author

Simran | TCGA RNA-seq learning project
