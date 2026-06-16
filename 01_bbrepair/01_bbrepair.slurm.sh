#!/bin/bash
#SBATCH --job-name=bbrepair
#SBATCH --account=fl3
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --array=0-11
#SBATCH --output=logs/bbrepair_%A_%a.out
#SBATCH --error=logs/bbrepair_%A_%a.err

# BBRepair: fix corrupted reads and remove unpaired mates.
# 128G / -Xmx120g required for 3KO samples (~67G peak usage at 64G).
# 3KO-022 will crash (structural FASTQ corruption beyond BBRepair's parser) —
# that sample is excluded from analysis.
# Submit: sbatch qc/01_bbrepair.slurm.sh

source config/project.env

SAMPLES_ARRAY=($(cat "${SAMPLES}"))
SAMPLE="${SAMPLES_ARRAY[$SLURM_ARRAY_TASK_ID]}"
N_SAMPLES=${#SAMPLES_ARRAY[@]}

echo "Processing sample: ${SAMPLE}"

module load bbmap/38.96--h5c4e2a8_0
# Check available versions with: module avail bbmap

STATS_DIR="${PROJECT_DIR}/qc/bbrepair_stats"
mkdir -p "${REPAIRED_OUT}" "${STATS_DIR}"

XMXG="${BBMAP_XMX:-120g}"

repair.sh \
    in1="${RAW_READS_DIR}/${SAMPLE}-R1.fastq.gz" \
    in2="${RAW_READS_DIR}/${SAMPLE}-R2.fastq.gz" \
    out1="${REPAIRED_OUT}/${SAMPLE}-R1.repaired.fastq.gz" \
    out2="${REPAIRED_OUT}/${SAMPLE}-R2.repaired.fastq.gz" \
    outs="${REPAIRED_OUT}/${SAMPLE}-singletons.fastq.gz" \
    tossbrokenreads tossjunk qin=33 -Xmx${XMXG} \
    2>&1 | tee "${STATS_DIR}/${SAMPLE}.repair.log"

echo "bbrepair complete for ${SAMPLE}"

# ── Per-sample stats (parse repair.sh log) ────────────────────────────────────
INPUT_READS=$(grep "^Input:" "${STATS_DIR}/${SAMPLE}.repair.log" | awk '{print $2}')
PAIRED_READS=$(grep "^Pairs:" "${STATS_DIR}/${SAMPLE}.repair.log" | awk '{print $2}')
DISCARDED=$(( INPUT_READS - PAIRED_READS ))
PCT=$(awk "BEGIN {printf \"%.2f\", (${PAIRED_READS:-0} / ${INPUT_READS:-1}) * 100}")

echo "${SAMPLE},${INPUT_READS},${PAIRED_READS},${DISCARDED},${PCT}" \
    > "${STATS_DIR}/${SAMPLE}.stats.csv"

echo "  Input reads:    ${INPUT_READS}"
echo "  Retained reads: ${PAIRED_READS}  (${PCT}%)"
echo "  Discarded:      ${DISCARDED}"

# ── Aggregate summary (whichever task finishes last) ──────────────────────────
N_DONE=$(ls "${STATS_DIR}/"*.stats.csv 2>/dev/null | wc -l)
if [ "${N_DONE}" -eq "${N_SAMPLES}" ]; then
    SUMMARY="${STATS_DIR}/repair_summary.csv"
    echo "sample,input_reads,paired_reads,discarded_reads,pct_paired_retained" > "${SUMMARY}"
    cat "${STATS_DIR}/"*.stats.csv >> "${SUMMARY}"
    echo ""
    echo "=== Repair summary ==="
    column -t -s ',' "${SUMMARY}"
fi
