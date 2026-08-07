#!/bin/bash
#SBATCH -J pairing
#SBATCH -p short
#SBATCH -t 00:15:00
#SBATCH -c 1
#SBATCH --mem=256M
#SBATCH -o /n/scratch/users/h/hua575/atb_filtered/gene_disruption_paired/logs/pairing_%A_%a.out
#SBATCH -e /n/scratch/users/h/hua575/atb_filtered/gene_disruption_paired/logs/pairing_%A_%a.err
#SBATCH --array=1-149

# Before submitting: run make_chunks.sh and set --array=1-<N> above to match
# the chunk count it reports. Optionally throttle concurrency, e.g. --array=1-149%50

set -euo pipefail

WORKDIR="/n/scratch/users/h/hua575/atb_filtered/gene_disruption_paired"
SCRIPT_DIR="/home/hua575/baymlab/mapping/pairing"

source /home/hua575/miniconda3/etc/profile.d/conda.sh
conda activate ismapper

TASK_ID=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")
CHUNK_FILE="$WORKDIR/chunks/chunk_${TASK_ID}.txt"

python3 "$SCRIPT_DIR/add_pairing_column.py" \
    --zip-list "$CHUNK_FILE" \
    --out-dir "/n/scratch/users/h/hua575/atb_filtered/gene_disruption"
