# a13

RNA-seq analysis of *Triticum aestivum* (cv. Halberd) response to *Parastagonospora nodorum* isolates: SN15 (wild-type), ToxA13/3KO (ToxA, Tox1, Tox3 triple-knockout mutant), and Tween-20 mock inoculation.

## Experimental design

| Condition | Samples |
|-----------|---------|
| Tween (mock control) | 012, 013, 015 |
| SN15 (wild-type) | 016, 017, 018, 019 |
| ToxA13/3KO (ToxA, Tox1, Tox3 triple-knockout) | 020, 021, 023 |
| Excluded | 022 (structural FASTQ corruption), 014 (PCA outlier — clusters with SN15) |

Paired-end 150 bp Illumina. Reference: *T. aestivum* Chinese Spring IWGSC RefSeq v2.1 (GCF_018294505.1).

---

## Repository structure

Steps 01–03 run on the Pawsey Setonix HPC. Steps 04–06 run locally in R.

```
01_bbrepair/          BBRepair read repair + SLURM logs + per-sample stats
02_fastp/             fastp trimming + per-sample QC reports + MultiQC output
03_kallisto/          Kallisto index + quantification + SLURM logs
  ├── kallisto_quant/ pseudoalignment output (one directory per sample)
reference/            BioMart GO annotations, NCBI→Ensembl ID map, GO curation
04_deseq2/            tximport + DESeq2 differential expression
  └── results/        DESeq2 CSVs, RDS objects, PCA plot
05_clusterprofiler/   GO enrichment analysis
  └── results/        GO_enrichment_*.csv (one per comparison)
06_plots/             GO enrichment figure + supporting category scripts
config/               sample list, Setonix paths (project.env), metadata.csv
```

---

## Running the pipeline

### HPC (Setonix) — steps 01–03

Edit `config/project.env` with your paths, then submit in order:

```bash
sbatch 01_bbrepair/01_bbrepair.slurm.sh
sbatch 02_fastp/02_fastp.slurm.sh
sbatch 03_kallisto/03a_kallisto_index.slurm.sh
sbatch 03_kallisto/03b_kallisto_quant.slurm.sh
```

Rsync `03_kallisto/kallisto_quant/` locally before proceeding.

### Local (R) — steps 04–06

Run from the repo root (`toxa13-wheat-rna/`):

```r
Rscript 04_deseq2/04a_metadata.r        # generates config/metadata.csv
Rscript 04_deseq2/04b_deseq2.r         # DESeq2 — slow (~10 min)
Rscript 04_deseq2/04c_pca.r            # PCA quality check
Rscript 04_deseq2/04d_expressed_genes.r
Rscript 05_clusterprofiler/05_go_enrichment.r  # GO enrichment — slow (~5 min)
Rscript 06_plots/06_go_plot.r          # produces final figure
```

---

## Key parameters

| Step | Tool | Key settings |
|------|------|--------------|
| Read repair | BBRepair | `tossbrokenreads qin=33 -Xmx60g` |
| Trimming | fastp | `-5 -3 -M 30` (5′/3′ quality window, mean Q30) |
| Quantification | Kallisto | `-b 100` bootstraps, 18 threads |
| DE testing | DESeq2 | `~ condition`, filter ≥10 reads in ≥2 samples |
| GO enrichment | clusterProfiler | LFC > 2, FDR < 0.05; `pvalueCutoff = 0.05` |
