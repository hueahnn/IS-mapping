#!/bin/bash
# BWA-align an IS element's cd-hit cluster-representative fasta (produced by
# scripts/cluster_and_align_by_is.py, {is_element}.clusters.query.fasta) against
# each of the 6 masked (no_IS_v2, ISfinder-BLAST-masked) E. coli reference
# genomes, one BAM per reference. Mirrors
# coverage_plots/align_to_masked_refs_v2.sbatch's bwa mem + samtools recipe,
# just pointed at an IS element's pooled cluster reps instead of
# pool_reps/merged_reps_min10.fasta.
#
# usage: IS_ELEMENT=IS1203 sbatch align_is_clusters_to_masked_refs_v2.sbatch
# (IS_ELEMENT defaults to IS1203 if not set)
#SBATCH -J align_is_clusters_v2
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -o /n/scratch/users/h/hua575/atb_filtered/is_level/logs/align_%A_%a.out
#SBATCH -e /n/scratch/users/h/hua575/atb_filtered/is_level/logs/align_%A_%a.err
#SBATCH --array=1-6

set -euo pipefail

IS_ELEMENT="${IS_ELEMENT:-IS1203}"
REF_DIR="/home/hua575/baymlab/mapping/ref_genomes/ecoli"
IS_DIR="/n/scratch/users/h/hua575/atb_filtered/is_level/$IS_ELEMENT"
QUERY_FASTA="$IS_DIR/$IS_ELEMENT.clusters.query.fasta"
OUT_DIR="$IS_DIR/bam"
THREADS=8

# the 6 masked reference genomes (GCF_000005845.2 is a duplicate RefSeq
# accession of GCA_000005845.2 and is intentionally excluded)
ACCESSIONS=(GCA_000005845.2 GCA_000163235.1 GCA_000692435.1 GCA_002473875.1 GCA_002966755.1 GCA_900096825.1)
ACC="${ACCESSIONS[$((SLURM_ARRAY_TASK_ID - 1))]}"

mkdir -p "$OUT_DIR"

source /home/hua575/miniconda3/etc/profile.d/conda.sh
conda activate ismapper

REF_FASTA="$REF_DIR/$ACC/${ACC}_no_IS_v2.fasta"
if [ ! -f "$REF_FASTA" ]; then
    echo "ERROR: missing reference fasta for $ACC: $REF_FASTA" >&2
    exit 1
fi
if [ ! -f "$QUERY_FASTA" ]; then
    echo "ERROR: missing cluster-representative query fasta for $IS_ELEMENT: $QUERY_FASTA" >&2
    exit 1
fi

echo "Aligning $QUERY_FASTA to $ACC..."
bwa mem -t "$THREADS" "$REF_FASTA" "$QUERY_FASTA" \
    | samtools sort -@ "$THREADS" -m 512M -o "$OUT_DIR/${ACC}.bam" -
samtools index "$OUT_DIR/${ACC}.bam"

echo "Done: $OUT_DIR/${ACC}.bam"
