# purpose: pool every accession's non-singleton cluster representative
# sequences (each accession's {accession}.reads.tsv, produced by rule
# clip_and_cluster's per-sample clustering step) across ALL IS elements and
# ALL genomes, split only by side (left/right) -- not by is_element -- then
# re-cluster each pooled side with cd-hit-est. This global re-clustering step
# is independent of whichever method rule clip_and_cluster used per-sample
# (cluster_overhangs_edlib.py or the older CD-HIT-based path both produce the
# same reads.tsv contract this pools from) -- it always uses cd-hit-est itself
# at the global level, so the CDHIT_ARGS below are no longer guaranteed to
# match the per-sample method's thresholds the way they did when both stages
# used CD-HIT. This collapses the same flanking sequence recurring across many
# genomes/IS elements into one global cluster, so downstream analysis can
# count how many distinct genes are actually being hit and how many times each
# one is disrupted.
#
# reads.tsv already drops singleton clusters (a read that didn't match
# anything else within its own genome) -- see rule clip_and_cluster -- so
# every representative sequence pooled here already represents a cluster of
# >=2 reads within its own genome. Singletons are not reconsidered at the
# global level.
#
# usage:
#   /home/hua575/miniconda3/envs/cd-hit/bin/python cluster_global_reps.py \
#       --output_dir /home/hua575/baymlab/mapping/atb/global_clusters \
#       --threads 16
#
# scale note: ~148,743 accession dirs, an estimated ~6.5M left-side and ~6.7M
# right-side cluster representatives to pool. Run this via sbatch, not on a
# login node -- pooling alone is I/O-bound over ~150k small files and can take
# well over an hour; cd-hit-est over millions of sequences will take longer
# still. --skip_pooling lets you rerun just the cd-hit-est step (e.g. to try
# different parameters) without repeating the expensive pooling pass.

import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

CD_HIT_EST_BIN = "/home/hua575/miniconda3/envs/cd-hit/bin/cd-hit-est"

# the parameters rule clip_and_cluster's OLD CD-HIT-based path used to use in
# /home/hua575/baymlab/mapping/Snakefile -- kept as this global step's own
# default since it still clusters with cd-hit-est regardless of which method
# produced the per-sample reads.tsv it's pooling from
CDHIT_ARGS = ["-c", "0.9", "-n", "8", "-G", "0", "-aS", "0.8", "-d", "0", "-M", "0"]

SIDES = ["left", "right"]


def pool_representative_seqs(cluster_dir: Path, output_dir: Path) -> dict:
    """
    Stream every accession's {accession}.reads.tsv into two pooled FASTAs, one
    per side, pooling clusters across every is_element and every accession.
    Written incrementally (not accumulated in memory) since the pooled total
    runs into the millions of sequences.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    pooled_paths = {side: output_dir / f"pooled_{side}.fasta" for side in SIDES}
    handles = {side: open(path, "w") for side, path in pooled_paths.items()}

    accession_dirs = sorted(d for d in cluster_dir.iterdir() if d.is_dir())
    n_clusters = {side: 0 for side in SIDES}
    try:
        for i, acc_dir in enumerate(accession_dirs, 1):
            reads_tsv = acc_dir / f"{acc_dir.name}.reads.tsv"
            if not reads_tsv.exists():
                continue
            df = pd.read_csv(reads_tsv, sep="\t", dtype={"pos": str})
            if df.empty:
                continue

            # one row per (is_element, side, cluster_id) -- reads.tsv repeats
            # the representative once per member read
            grouped = df.groupby(
                ["is_element", "side", "cluster_id"], as_index=False
            ).agg(rep_seq=("representative_sequence", "first"))

            for row in grouped.itertuples(index=False):
                header = f"{acc_dir.name}__{row.is_element}__{row.side}__cluster{row.cluster_id}"
                handles[row.side].write(f">{header}\n{row.rep_seq}\n")
                n_clusters[row.side] += 1

            if i % 5000 == 0:
                print(f"[{i}/{len(accession_dirs)}] accessions processed "
                      f"({n_clusters['left']} left, {n_clusters['right']} right "
                      f"clusters pooled so far)", file=sys.stderr, flush=True)
    finally:
        for fh in handles.values():
            fh.close()

    print(f"done pooling: {n_clusters['left']} left-side and "
          f"{n_clusters['right']} right-side cluster representatives from "
          f"{len(accession_dirs)} accession dirs", file=sys.stderr)
    return pooled_paths


def run_cdhit_est(fasta_path: Path, output_dir: Path, side: str, threads: int) -> None:
    out_path = output_dir / f"global_{side}.cdhit"
    cmd = [CD_HIT_EST_BIN, "-i", str(fasta_path), "-o", str(out_path),
           *CDHIT_ARGS, "-T", str(threads)]
    print(f"running: {' '.join(cmd)}", file=sys.stderr, flush=True)
    subprocess.run(cmd, check=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cluster_dir", type=Path,
                    default=Path("/n/scratch/users/h/hua575/atb_filtered/clusters"),
                    help="per-accession overhang-clustering output dir "
                         "(contains {accession}/{accession}.reads.tsv)")
    p.add_argument("--output_dir", type=Path,
                    default=Path("/n/scratch/users/h/hua575/atb_filtered/global_cdhit"),
                    help="where to write pooled_{side}.fasta and "
                         "global_{side}.cdhit(.clstr) "
                         "(default: atb_filtered/global_cdhit, alongside clusters/)")
    p.add_argument("--threads", type=int, default=16,
                    help="passed to cd-hit-est -T")
    p.add_argument("--skip_pooling", action="store_true",
                    help="reuse existing pooled_{side}.fasta in --output_dir "
                         "and only (re-)run cd-hit-est")
    args = p.parse_args()

    if args.skip_pooling:
        pooled_paths = {side: args.output_dir / f"pooled_{side}.fasta" for side in SIDES}
        missing = [str(path) for path in pooled_paths.values() if not path.exists()]
        if missing:
            sys.exit(f"--skip_pooling but missing pooled fasta(s): {missing}")
    else:
        t0 = time.time()
        pooled_paths = pool_representative_seqs(args.cluster_dir, args.output_dir)
        print(f"pooling took {time.time() - t0:.0f}s", file=sys.stderr)

    for side in SIDES:
        run_cdhit_est(pooled_paths[side], args.output_dir, side, args.threads)


if __name__ == "__main__":
    main()
