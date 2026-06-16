#!/usr/bin/env Rscript
# Expressed vs differential gene analysis
# Input:  04_deseq2/results/txi.rds, 04_deseq2/results/deseq_results_*.csv
# Output: 04_deseq2/results/summary_*.csv, 04_deseq2/results/gene_list_*.csv
# Run from: toxa13-wheat-rna/

library(readr)
library(dplyr)

cat("=== Expressed vs Differential Genes Analysis ===\n\n")

txi      <- readRDS("04_deseq2/results/txi.rds")
metadata <- read.csv("config/metadata.csv")
counts   <- txi$counts

sn15_samples  <- metadata$sample[metadata$condition == "SN15"]
ko3_samples   <- metadata$sample[metadata$condition == "3KO"]
tween_samples <- metadata$sample[metadata$condition == "Tween"]

expressed <- function(mat, cols) {
  rownames(mat)[rowSums(mat[, cols, drop = FALSE] >= 10) >= 1]
}

genes_in_sn15  <- expressed(counts, sn15_samples)
genes_in_3ko   <- expressed(counts, ko3_samples)
genes_in_tween <- expressed(counts, tween_samples)
genes_passing  <- rownames(counts)[rowSums(counts >= 10) >= 2]
genes_any      <- unique(c(genes_in_sn15, genes_in_3ko, genes_in_tween))

cat("=== Table 1: Gene Expression by Condition ===\n\n")

n <- nrow(counts)
expression_by_condition <- data.frame(
  Condition = c(
    "SN15", "3KO", "Tween",
    "Any condition (union)",
    "Passed DESeq2 filter (>=10 in >=2 samples)"
  ),
  Transcripts_Expressed = c(
    length(genes_in_sn15), length(genes_in_3ko), length(genes_in_tween),
    length(genes_any), length(genes_passing)
  ),
  Percent_of_Reference = round(100 * c(
    length(genes_in_sn15), length(genes_in_3ko), length(genes_in_tween),
    length(genes_any), length(genes_passing)
  ) / n, 1)
)

write_csv(expression_by_condition, "04_deseq2/results/summary_expression_by_condition.csv")
print(expression_by_condition)
cat("\n")

cat("Loading DESeq2 results...\n")

deg_ids <- function(df) {
  df |>
    filter(abs(log2FoldChange) > 2 & padj < 0.05) |>
    pull(transcript_id)
}

deseq_sn15     <- read_csv("04_deseq2/results/deseq_results_SN15_vs_Tween.csv", show_col_types = FALSE)
deseq_3ko      <- read_csv("04_deseq2/results/deseq_results_3KO_vs_Tween.csv",  show_col_types = FALSE)
deseq_sn15_3ko <- read_csv("04_deseq2/results/deseq_results_3KO_vs_SN15.csv",   show_col_types = FALSE)

sn15_deg     <- deg_ids(deseq_sn15)
ko3_deg      <- deg_ids(deseq_3ko)
sn15_3ko_deg <- deg_ids(deseq_sn15_3ko)
all_tested   <- deseq_sn15$transcript_id

cat("=== Table 2: Expressed vs Differential by Comparison ===\n\n")

expr_vs_de <- data.frame(
  Comparison         = c("SN15 vs Tween", "3KO vs Tween", "3KO vs SN15"),
  Genes_Tested       = length(all_tested),
  Genes_Differential = c(length(sn15_deg), length(ko3_deg), length(sn15_3ko_deg)),
  Genes_Non_DE = c(
    length(all_tested) - length(sn15_deg),
    length(all_tested) - length(ko3_deg),
    length(all_tested) - length(sn15_3ko_deg)
  ),
  Percent_DE = round(100 * c(
    length(sn15_deg), length(ko3_deg), length(sn15_3ko_deg)
  ) / length(all_tested), 1)
)

write_csv(expr_vs_de, "04_deseq2/results/summary_expressed_vs_differential.csv")
print(expr_vs_de)
cat("\n")

all_degs           <- unique(c(sn15_deg, ko3_deg, sn15_3ko_deg))
never_differential <- setdiff(all_tested, all_degs)

cat("=== Table 3: Never Differential Across All Comparisons ===\n\n")

never_de_summary <- data.frame(
  Category = c(
    "Total transcripts in reference",
    "Transcripts tested in DESeq2 (passed filter)",
    "Transcripts DE in at least one comparison",
    "Transcripts NEVER differential",
    "Percent never differential (of tested)"
  ),
  Count = c(
    nrow(counts), length(all_tested), length(all_degs),
    length(never_differential),
    round(100 * length(never_differential) / length(all_tested), 1)
  )
)

write_csv(never_de_summary, "04_deseq2/results/summary_never_differential.csv")
print(never_de_summary)
cat("\n")

write_csv(data.frame(transcript_id = genes_in_sn15),      "04_deseq2/results/gene_list_expressed_in_SN15.csv")
write_csv(data.frame(transcript_id = genes_in_3ko),       "04_deseq2/results/gene_list_expressed_in_3KO.csv")
write_csv(data.frame(transcript_id = genes_in_tween),     "04_deseq2/results/gene_list_expressed_in_Tween.csv")
write_csv(data.frame(transcript_id = never_differential), "04_deseq2/results/gene_list_never_differential.csv")
write_csv(data.frame(transcript_id = all_degs),           "04_deseq2/results/gene_list_differential_any_comparison.csv")

cat("Saved all tables and gene lists to 04_deseq2/results/\n")
