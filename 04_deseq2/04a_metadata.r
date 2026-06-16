#!/usr/bin/env Rscript
# Generate config/metadata.csv for DESeq2 tximport
# Excludes 3KO-022 (structural FASTQ corruption, PCA outlier)
#          Tween-014 (clusters with SN15 in PCA — likely contamination)
# Run from: toxa13-wheat-rna/

library(readr)

samples <- c(
  "inplanta-3KO-020",   "inplanta-3KO-021",   "inplanta-3KO-023",
  "inplanta-SN15-016",  "inplanta-SN15-017",  "inplanta-SN15-018",  "inplanta-SN15-019",
  "inplanta-Tween-012", "inplanta-Tween-013", "inplanta-Tween-015"
)

conditions <- c(
  "3KO",  "3KO",  "3KO",
  "SN15", "SN15", "SN15", "SN15",
  "Tween","Tween","Tween"
)

metadata <- data.frame(
  sample    = samples,
  condition = conditions,
  path      = file.path("03_kallisto/kallisto_quant", samples),
  stringsAsFactors = FALSE
)

missing <- !dir.exists(metadata$path)
if (any(missing)) {
  cat("ERROR: Missing kallisto output directories:\n")
  print(metadata$path[missing])
  stop("Run 03b_kallisto_quant.slurm.sh on Setonix and download to 03_kallisto/kallisto_quant/.")
}

write_csv(metadata, "config/metadata.csv")

cat("Metadata written to: config/metadata.csv\n")
cat("Samples included:\n")
print(table(metadata$condition))
