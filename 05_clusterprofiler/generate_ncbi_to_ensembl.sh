#!/usr/bin/env bash
# Generate ncbi_to_ensembl.tsv from the NCBI CDS FASTA on Setonix
# Output: reference/ncbi_to_ensembl.tsv
# Run from: toxa13-wheat-rna/ after downloading the FASTA
#
# Each FASTA header looks like:
#   >lcl|NC_057794.1_cds_XP_044452222.1_1 [gene=LOC123183831] [protein=...] ...
# We extract: kallisto target_id (after >) and LOC gene ID (from [gene=] tag).
#
# Run on Setonix (file is ~2GB), then rsync the output TSV locally:
#   rsync llenzo@setonix.pawsey.org.au:<project>/reference/ncbi_to_ensembl.tsv \
#         reference/ncbi_to_ensembl.tsv

source config/project.env   # sets CDS_FASTA

grep "^>" "${CDS_FASTA}" | \
  grep "\[gene=LOC" | \
  awk '{
    target = substr($1, 2)
    for (i = 2; i <= NF; i++) {
      if ($i ~ /^\[gene=LOC/) {
        gene = $i
        gsub(/^\[gene=/, "", gene)
        gsub(/\]$/,      "", gene)
        print target "\t" gene
      }
    }
  }' > reference/ncbi_to_ensembl.tsv

echo "Done. $(wc -l < reference/ncbi_to_ensembl.tsv) entries written."
