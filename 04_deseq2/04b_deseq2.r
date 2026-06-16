#!/usr/bin/env Rscript
# DESeq2 differential expression analysis
# Input:  config/metadata.csv + 03_kallisto/kallisto_quant/*/abundance.tsv
# Output: 04_deseq2/results/txi.rds, dds.rds, vsd.rds
#         04_deseq2/results/deseq_results_{SN15_vs_Tween,3KO_vs_Tween,3KO_vs_SN15}.csv
# Run from: toxa13-wheat-rna/
# Requires: 04a_metadata.r to have been run first

library(DESeq2)
library(tximport)
library(readr)
library(dplyr)

dir.create("04_deseq2/results", showWarnings = FALSE)

cat("=== Loading Kallisto Data ===\n\n")

metadata <- read.csv("config/metadata.csv")
metadata$condition <- factor(
  metadata$condition, levels = c("Tween", "SN15", "3KO")
)

cat("Samples:\n")
print(table(metadata$condition))
cat("\n")

files <- file.path(metadata$path, "abundance.tsv")
names(files) <- metadata$sample

missing <- !file.exists(files)
if (any(missing)) {
  stop("Missing abundance.tsv:\n", paste(files[missing], collapse = "\n"))
}

cat("Importing kallisto quantifications...\n")
txi <- tximport(files, type = "kallisto", txOut = TRUE)
cat("Transcripts imported:", nrow(txi$counts), "\n\n")

saveRDS(txi, "04_deseq2/results/txi.rds")

cat("=== DESeq2 Analysis ===\n\n")

dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)

cat("Transcripts before filtering:", nrow(dds), "\n")
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
cat("Transcripts after filtering (>=10 reads in >=2 samples):",
    nrow(dds), "\n\n")

cat("Running DESeq2...\n")
dds <- DESeq(dds)
cat("DESeq2 complete\n\n")

cat("Running variance stabilising transformation...\n")
vsd <- vst(dds, blind = FALSE)

saveRDS(dds, "04_deseq2/results/dds.rds")
saveRDS(vsd, "04_deseq2/results/vsd.rds")
cat("Saved: 04_deseq2/results/txi.rds, dds.rds, vsd.rds\n\n")

cat("=== Extracting Results ===\n\n")

comparisons <- list(
  SN15_vs_Tween = c("SN15", "Tween"),
  "3KO_vs_Tween" = c("3KO",  "Tween"),
  "3KO_vs_SN15" = c("3KO",  "SN15")
)

for (comp_name in names(comparisons)) {
  num <- comparisons[[comp_name]][1]
  den <- comparisons[[comp_name]][2]
  cat("Comparison:", comp_name, "(", num, "vs", den, ")\n")

  res    <- results(dds, contrast = c("condition", num, den))
  res_df <- as.data.frame(res)
  res_df$transcript_id <- rownames(res_df)

  outfile <- paste0("04_deseq2/results/deseq_results_", comp_name, ".csv")
  write_csv(res_df, outfile)

  up   <- sum(res_df$log2FoldChange >  2 & res_df$padj < 0.05, na.rm = TRUE)
  down <- sum(res_df$log2FoldChange < -2 & res_df$padj < 0.05, na.rm = TRUE)
  cat("  Up (LFC>2, FDR<0.05):", up, "  Down:", down, "\n")
  cat("  Saved:", outfile, "\n\n")
}

cat("=== DESeq2 Complete ===\n")
