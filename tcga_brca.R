BiocManager::install("SummarizedExperiment", force=TRUE)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

BiocManager::install("GenomicRanges", force=TRUE)
BiocManager::install("S4Vectors", force=TRUE)
library(GenomicRanges)
library(S4Vectors)

R.version$version.string

BiocManager::install(c(
  "S4Vectors",
  "GenomicRanges", 
  "SummarizedExperiment",
  "DESeq2",
  "TCGAbiolinks"
), force= TRUE, update = TRUE)

R.version$version.string

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
). #Tcga data query banai



metadata <- getResults(query)

tumor_meta  <- metadata[metadata$sample_type == "Primary Tumor", ]
normal_meta <- metadata[metadata$sample_type == "Solid Tissue Normal", ]

nrow(tumor_meta)   # verify
nrow(normal_meta)  # verify

set.seed(42)

tumor_40  <- tumor_meta[sample(nrow(tumor_meta),  40), ]
normal_40 <- normal_meta[sample(nrow(normal_meta), 40), ]

# Combine karo
selected_80 <- rbind(tumor_40, normal_40)

# Verify
nrow(selected_80)
table(selected_80$sample_type)

selected_barcodes <- selected_80$cases
# ya agar cases nahi hai toh:
selected_barcodes <- selected_80$barcode

head(selected_barcodes)  # dekhlo kaisa dikhta hai
length(selected_barcodes) # 80 hone chahiye

query_80 <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  barcode = selected_barcodes    # ← yahan filter ho gaya!
)

# Verify — 80 hone chahiye 
nrow(getResults(query_80))

GDCdownload(
  query_80,
  method = "api",
  files.per.chunk = 10
)

se <- GDCprepare(query_80)

# Verify
dim(se)
table(se$sample_type)   # 40 tumor, 40 normal

# Counts nikalo
counts_mat <- assay(se, "unstranded")

# colData banao
col_info <- data.frame(
  sample_type = se$sample_type,
  row.names = colnames(se)
)

col_info$sample_type <- factor(
  col_info$sample_type,
  levels = c("Solid Tissue Normal", "Primary Tumor")
)

# DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = col_info,
  design = ~ sample_type
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 1
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds) 

# Results
res <- results(dds,
               contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"),
               alpha = 0.05
)

summary(res)
res_clean <- res[!is.na(res$padj), ]
nrow(res_clean)  # kitne genes bache

# Significant = padj < 0.05 AND |log2FC| > 1
res_sig <- res_clean[
  res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1,
]

nrow(res_sig)

# Tumor mein ZYADA expressed
upregulated <- res_sig[res_sig$log2FoldChange > 1, ]
upregulated <- upregulated[order(upregulated$log2FoldChange, decreasing = TRUE), ]

# Tumor mein KAM expressed  
downregulated <- res_sig[res_sig$log2FoldChange < -1, ]
downregulated <- downregulated[order(downregulated$log2FoldChange), ]

# Count
nrow(upregulated)
nrow(downregulated)

# Top genes dekho
head(upregulated, 10)
head(downregulated, 10)

# rowData mein gene info hota hai
gene_info <- rowData(se)
colnames(gene_info)   # dekho kya available hai
head(gene_info)

# Mapping banao: ENSEMBL ID → Gene Symbol
id_to_symbol <- setNames(
  gene_info$gene_name,
  rownames(gene_info)
)

# Res_clean mein symbol add karo
res_clean$symbol <- id_to_symbol[rownames(res_clean)]

# Verify
head(res_clean[, c("log2FoldChange", "padj", "symbol")])

BiocManager::install("fgsea")
install.packages("msigdbr")   # pathway gene sets ke liye

library(fgsea)
library(msigdbr)


# Ranking metric = sign(log2FC) * -log10(pvalue)
# Yeh direction aur significance dono capture karta hai

ranked_genes <- res_clean$log2FoldChange * -log10(res_clean$pvalue + 1e-300)
# 1e-300 add kiya taaki log10(0) se bachein

names(ranked_genes) <- res_clean$symbol

# NA symbols hatao
ranked_genes <- ranked_genes[!is.na(names(ranked_genes))]

# Sort karo — ZAROORI hai fgsea ke liye
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Verify
head(ranked_genes, 5)
tail(ranked_genes, 5)
length(ranked_genes)


# Hallmark pathways — beginner ke liye best
# 50 well-defined biological pathways hain
pathways_h <- msigdbr(
  species = "Homo sapiens",
  collection = "H"          # H = Hallmark
)

# fgsea format mein convert karo — named list chahiye
pathway_list <- split(
  pathways_h$gene_symbol,
  pathways_h$gs_name
)

# Dekho kitne pathways hain
length(pathway_list)        # 50 hone chahiye
pathway_list[[1]][1:10]     # ek pathway ke genes dekho


set.seed(42)   # results reproducible honge

fgsea_res <- fgsea(
  pathways = pathway_list,
  stats    = ranked_genes,
  minSize  = 15,     # pathway mein minimum genes
  maxSize  = 500     # pathway mein maximum genes
)

# Results dekho
head(fgsea_res[order(fgsea_res$padj), ], 20)

# Duplicates check karo
sum(duplicated(names(ranked_genes)))

# Duplicates hatao — sirf pehla wala rakho
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

# Verify
sum(duplicated(names(ranked_genes)))  # 0 hona chahiye


# Step 1: Duplicates hatao
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

# Step 2: fgsea dobara
set.seed(42)
fgsea_res <- fgsea(
  pathways = pathway_list,
  stats    = ranked_genes,
  minSize  = 15,
  maxSize  = 500
)

# Step 3: Results dekho
head(fgsea_res[order(fgsea_res$padj), ], 20)

# Significant pathways
sig_pathways <- fgsea_res[fgsea_res$padj < 0.05, ]
sig_pathways <- sig_pathways[order(sig_pathways$NES, decreasing = TRUE), ]

# Upregulated pathways (NES > 0 = enriched in Tumor)
up_pathways   <- sig_pathways[sig_pathways$NES > 0, ]

# Downregulated pathways (NES < 0 = depleted in Tumor)
down_pathways <- sig_pathways[sig_pathways$NES < 0, ]

nrow(up_pathways)
nrow(down_pathways)

