#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --account=fl3
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-11
#SBATCH --output=logs/fastp_%A_%a.out
#SBATCH --error=logs/fastp_%A_%a.err

# fastp: adapter trimming and quality filtering on BBRepair output
# Depends on step 01 (BBRepair) — reads repaired fastqs
# Submit: sbatch --dependency=afterok:<bbrepair_jobid> qc/02_fastp.slurm.sh

source config/project.env

SAMPLES_ARRAY=($(cat "${SAMPLES}"))
SAMPLE="${SAMPLES_ARRAY[$SLURM_ARRAY_TASK_ID]}"

echo "Processing sample: ${SAMPLE}"

module load fastp/0.23.4-5dugkew

mkdir -p "${FASTP_OUT}"
mkdir -p logs/fastp_reports

fastp \
    --in1  "${REPAIRED_OUT}/${SAMPLE}-R1.repaired.fastq.gz" \
    --in2  "${REPAIRED_OUT}/${SAMPLE}-R2.repaired.fastq.gz" \
    --out1 "${FASTP_OUT}/${SAMPLE}-R1.trimmed.fastq.gz" \
    --out2 "${FASTP_OUT}/${SAMPLE}-R2.trimmed.fastq.gz" \
    --thread 8 \
    -5 -3 -M 30 \
    --html "logs/fastp_reports/${SAMPLE}.fastp.html" \
    --json "logs/fastp_reports/${SAMPLE}.fastp.json" \
    --verbose

echo "fastp complete for ${SAMPLE}"
