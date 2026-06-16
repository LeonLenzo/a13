#!/usr/bin/env Rscript
# GO enrichment analysis using clusterProfiler
# Input:  03_kallisto/reference/ncbi_to_ensembl.tsv
#         03_kallisto/reference/mart_export.txt
#         04_deseq2/results/deseq_results_*.csv
# Output: 05_clusterprofiler/results/GO_enrichment_*.csv
# Run from: toxa13-wheat-rna/
# Run 06_plots/06_go_plot.r after this to generate the figure.
#
# ID mapping chain:
#   NCBI target_id -> LOC gene ID  (ncbi_to_ensembl.tsv, strip "LOC" prefix)
#   Numeric entrez ID -> Ensembl transcript ID  (mart_export.txt col 9 -> col 2)

library(clusterProfiler)
library(readr)
library(dplyr)

dir.create("05_clusterprofiler/results", showWarnings = FALSE)

cat("=== Loading Reference Data ===\n\n")

id_map <- read_tsv("reference/ncbi_to_ensembl.tsv",
                   col_names = c("transcript_id", "loc_gene_id"),
                   show_col_types = FALSE) %>%
  filter(grepl("^LOC", loc_gene_id)) %>%
  mutate(entrez_id = sub("LOC", "", loc_gene_id))

cat("NCBI ID map loaded:", nrow(id_map), "entries\n")

go_raw <- read_tsv("reference/mart_export.txt", show_col_types = FALSE)

entrez_to_ensembl <- go_raw %>%
  select(
    ensembl_transcript = `Transcript stable ID`,
    entrez_id          = `NCBI gene (formerly Entrezgene) ID`
  ) %>%
  filter(!is.na(entrez_id) & entrez_id != "") %>%
  mutate(entrez_id = as.character(entrez_id)) %>%
  distinct()

cat("Entrez -> Ensembl transcript mapping:", nrow(entrez_to_ensembl), "entries\n")

go_clean <- go_raw %>%
  filter(!is.na(`GO term accession`) & `GO term accession` != "") %>%
  select(
    Gene_ID        = `Transcript stable ID`,
    GO_Terms       = `GO term accession`,
    GO_Description = `GO term name`,
    GO_Category    = `GO domain`
  ) %>%
  mutate(
    Ontology = case_when(
      GO_Category == "biological_process" ~ "BP",
      GO_Category == "molecular_function" ~ "MF",
      GO_Category == "cellular_component" ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Ontology)) %>%
  distinct()

cat("GO annotations:", nrow(go_clean), "rows,",
    length(unique(go_clean$Gene_ID)), "unique transcripts\n\n")

go2gene       <- go_clean %>% select(GO_Terms, Gene_ID) %>% distinct()
go2desc       <- go_clean %>% select(GO_Terms, GO_Description, GO_Category, Ontology) %>% distinct()
gene_universe <- unique(go_clean$Gene_ID)

to_ensembl <- function(ncbi_ids) {
  id_map %>%
    filter(transcript_id %in% ncbi_ids) %>%
    inner_join(entrez_to_ensembl, by = "entrez_id", relationship = "many-to-many") %>%
    pull(ensembl_transcript) %>%
    unique() %>%
    intersect(gene_universe)
}

cat("=== GO Enrichment ===\n\n")

log2fc_threshold <- 2
padj_threshold   <- 0.05
comparisons      <- c("SN15_vs_Tween", "3KO_vs_Tween", "3KO_vs_SN15")

run_enricher <- function(genes, direction, universe) {
  if (length(genes) == 0) return(NULL)
  ego <- enricher(gene = genes, TERM2GENE = go2gene, universe = universe,
                  pvalueCutoff = 0.05, qvalueCutoff = 0.1, minGSSize = 3)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  as.data.frame(ego) %>%
    mutate(ID = as.character(ID), Direction = direction) %>%
    left_join(go2desc, by = c("ID" = "GO_Terms")) %>%
    rename(GO_ID = ID)
}

for (comp_name in comparisons) {
  cat("Processing:", comp_name, "\n")

  degs <- read_csv(paste0("04_deseq2/results/deseq_results_", comp_name, ".csv"),
                   show_col_types = FALSE)

  up_ids   <- degs %>% filter(log2FoldChange >  log2fc_threshold & padj < padj_threshold) %>% pull(transcript_id)
  down_ids <- degs %>% filter(log2FoldChange < -log2fc_threshold & padj < padj_threshold) %>% pull(transcript_id)

  up_genes   <- to_ensembl(up_ids)
  down_genes <- to_ensembl(down_ids)

  cat("  Up in GO universe:", length(up_genes), " Down:", length(down_genes), "\n")

  combined <- bind_rows(
    run_enricher(up_genes,   "Up",   gene_universe),
    run_enricher(down_genes, "Down", gene_universe)
  )

  if (nrow(combined) > 0) {
    combined$Dataset <- comp_name
    outfile <- paste0("05_clusterprofiler/results/GO_enrichment_", comp_name, ".csv")
    write_csv(combined, outfile)
    cat("  Enriched GO terms:", nrow(combined), "->", outfile, "\n")
  } else {
    cat("  No significant enrichment found\n")
  }
  cat("\n")
}

cat("=== GO Enrichment Complete ===\n")
cat("Run 06_plots/06_go_plot.r to generate the figure.\n")
