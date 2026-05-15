library(ggplot2)
library(pheatmap)
library(DESeq2)
library(fgsea)

install.packages("pheatmap")

res_df <- as.data.frame(res_clean)
colnames(res_df)
res_df$symbol<- id_to_symbol[rownames(res_df)]

res_df$color <- "Not Significant"
res_df$color[res_df$log2FoldChange >1 & res_df$padj<0.05] <- "Upregulated"
res_df$color[res_df$log2FoldChange < -1 & res_df$padj <0.05] <- "Downregulated"
table(res_df$color)

# To label the Top genes 
top_genes <- res_df[res_df$color != "Not significant", ]
top_genes <- top_genes[order(top_genes$padj), ][1:15, ]

# Plot
volcano <- ggplot(res_df,
                  aes(x = log2FoldChange,
                      y = -log10(padj),
                      color = color)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c(
    "Upregulated"     = "firebrick",
    "Downregulated"   = "steelblue",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             alpha = 0.4) +
  geom_text(data = top_genes,
            aes(label = symbol),
            size = 2.5,
            hjust = -0.1,
            check_overlap = TRUE) +
  theme_minimal() +
  labs(
    title    = "Tumor vs Normal — TCGA BRCA",
    x        = "log2 Fold Change",
    y        = "-log10(adjusted p-value)",
    color    = ""
  ) +
  theme(legend.position = "bottom")

volcano

# Save it
ggsave("volcano_brca.png",
       plot   = volcano,
       width  = 8,
       height = 6,
       dpi    = 300)

#Heatmap
# Top 50 significant genes lo
top50 <- res_sig[order(res_sig$padj), ][1:50, ]
top50_genes <- rownames(top50)

# Normalized counts (vst already in QC)
# If vst wasnot done:
vsd <- vst(dds, blind = FALSE)

# Counts matrix subsetting
heatmap_mat <- assay(vsd)[top50_genes, ]

# Gene symbols as row names
rownames(heatmap_mat) <- id_to_symbol[rownames(heatmap_mat)]

# Making the Sample annotation 
annotation_col <- data.frame(
  Type = col_info$sample_type,
  row.names = rownames(col_info)
)

# Colors
ann_colors <- list(
  Type = c(
    "Primary Tumor"       = "firebrick",
    "Solid Tissue Normal" = "steelblue"
  )
)

# Plot
pheatmap(heatmap_mat,
         annotation_col  = annotation_col,
         annotation_colors = ann_colors,
         scale           = "row",      # gene-wise z-score
         show_colnames   = FALSE,      # 80 sample names — too many
         fontsize_row    = 7,
         cluster_cols    = TRUE,       # samples cluster
         cluster_rows    = TRUE,       # genes cluster 
         main            = "Top 50 DE Genes — TCGA BRCA",
         filename        = "heatmap_brca.png",
         width           = 10,
         height          = 8
)

pheatmap(heatmap_mat,
         annotation_col = annotation_col,
         scale          = "row",
         show_colnames  = FALSE,
         fontsize_row   = 7,
         main           = "Top 50 DE Genes"
)

# fgsea — minimum version
plotEnrichment(
  pathway_list[["HALLMARK_E2F_TARGETS"]],
  ranked_genes)

# Top pathway enrichment plot
plotEnrichment(
  pathway_list[["HALLMARK_E2F_TARGETS"]],
  ranked_genes
) +
  labs(
    title = "HALLMARK_E2F_TARGETS",
    x     = "Rank in gene list",
    y     = "Enrichment score"
  )

ggsave("fgsea_e2f.png", width = 7, height = 4, dpi = 300)

# Table plot for all significant pathways
top20_paths <- sig_pathways[
  order(abs(sig_pathways$NES), decreasing = TRUE),
]$pathway[1:20]

png("fgsea_table.png", width = 1200, height = 800, res = 120)
plotGseaTable(
  pathway_list[top20_paths],
  ranked_genes,
  fgsea_res,
  gseaParam = 0.5
)
dev.off()

pheatmap(heatmap_mat,
         annotation_col = annotation_col,
         scale          = "row",
         show_colnames  = FALSE)
