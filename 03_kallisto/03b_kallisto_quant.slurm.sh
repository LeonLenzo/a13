#!/bin/bash
#SBATCH --job-name=kallisto_quant
#SBATCH --account=fl3
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --array=0-11
#SBATCH --output=logs/kallisto_quant_%A_%a.out
#SBATCH --error=logs/kallisto_quant_%A_%a.err

# Kallisto pseudoalignment: quantify transcript abundances
# Runs as SLURM array job — one task per sample
# Depends on: step 02 (BBRepair) and step 04 (Kallisto index)
# Submit: sbatch --dependency=afterok:<bbrepair_jobid>:<index_jobid> align/05_kallisto_quant.slurm.sh

source config/project.env

SAMPLES_ARRAY=($(cat "${SAMPLES}"))
SAMPLE="${SAMPLES_ARRAY[$SLURM_ARRAY_TASK_ID]}"

echo "Running Kallisto quant for sample: ${SAMPLE}"

module load spack/0.23.1
spack load kallisto
# Check available versions with: spack find kallisto

SAMPLE_OUT="${KALLISTO_OUT}/${SAMPLE}"
mkdir -p "${SAMPLE_OUT}"

kallisto quant \
    -i "${KALLISTO_INDEX}" \
    -o "${SAMPLE_OUT}" \
    -b 100 \
    -t 18 \
    "${FASTP_OUT}/${SAMPLE}-R1.trimmed.fastq.gz" \
    "${FASTP_OUT}/${SAMPLE}-R2.trimmed.fastq.gz"

echo "Kallisto quant complete for ${SAMPLE}"
echo "Outputs: ${SAMPLE_OUT}/abundance.tsv, ${SAMPLE_OUT}/abundance.h5, ${SAMPLE_OUT}/run_info.json"
