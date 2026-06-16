#!/usr/bin/env Rscript
# PCA plot for sample clustering quality check
# Input:  04_deseq2/results/vsd.rds (from 04b_deseq2.r)
# Output: 04_deseq2/results/PCA_plot.{png,pdf}
# Run from: toxa13-wheat-rna/

library(DESeq2)
library(ggplot2)

cat("Loading VST data...\n")
vsd <- readRDS("04_deseq2/results/vsd.rds")

cat("=== PCA Plot ===\n\n")

pca_data    <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -0.8, hjust = 0.5, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 10),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(
    name   = "Treatment",
    values = c("3KO" = "#E41A1C", "SN15" = "#377EB8", "Tween" = "#4DAF4A")
  )

ggsave("04_deseq2/results/PCA_plot.png", p, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("04_deseq2/results/PCA_plot.pdf", p, width = 8, height = 6, bg = "white")

cat("Saved: 04_deseq2/results/PCA_plot.png, 04_deseq2/results/PCA_plot.pdf\n\n")
cat("PCA coordinates:\n")
print(pca_data)
