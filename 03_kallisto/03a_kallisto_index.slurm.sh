#!/bin/bash
#SBATCH --job-name=kallisto_index
#SBATCH --account=fl3
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/kallisto_index_%j.out
#SBATCH --error=logs/kallisto_index_%j.err

# Build Kallisto index from IWGSC RefSeq v2.1 CDS FASTA
# Reference: GCF_018294505.1 (Triticum aestivum Chinese Spring)
# Download CDS FASTA from NCBI before running this script (see config/project.env for filename)
# Submit: sbatch reference/04_kallisto_index.slurm.sh

source config/project.env

module load kallisto/0.50.1--h6de1650_2
# Check available versions with: module avail kallisto

mkdir -p "${PROJECT_DIR}/reference"

if [[ ! -f "${CDS_FASTA}" ]]; then
    echo "ERROR: CDS FASTA not found at ${CDS_FASTA}"
    echo "Download from NCBI: GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_cds_from_genomic.fna.gz"
    exit 1
fi

kallisto index \
    -i "${KALLISTO_INDEX}" \
    "${CDS_FASTA}"

echo "Kallisto index written to: ${KALLISTO_INDEX}"
