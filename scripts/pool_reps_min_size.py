# purpose: pool representative sequences for every (is_element, side,
# cluster_id) cluster across ALL accessions in the per-accession cd-hit
# directory (/n/scratch/users/h/hua575/atb_filtered/cd-hit/{accession}/
# {accession}.reads.tsv, produced by combine_cdhit_clusters.py) into a single
# merged fasta, keeping only clusters whose size (number of member reads) is
# >= --min_size. Unlike cluster_global_reps.py this does NOT re-run cd-hit or
# split by side -- it just filters and merges the existing per-accession
# cluster representatives directly.
#
# usage (single chunk):
#   python pool_reps_min_size.py --accession-list chunk.txt --min_size 10 \
#       --out chunk_out.fasta
#
# usage (all accessions, no chunking -- only for small test runs):
#   python pool_reps_min_size.py --cdhit_dir /n/scratch/users/h/hua575/atb_filtered/cd-hit \
#       --min_size 10 --out merged.fasta

import argparse
import sys
from pathlib import Path

import pandas as pd


def process_accession(acc_dir: Path, min_size: int, out_fh) -> int:
    reads_tsv = acc_dir / f"{acc_dir.name}.reads.tsv"
    if not reads_tsv.exists():
        return 0

    df = pd.read_csv(reads_tsv, sep="\t", dtype={"pos": str})
    if df.empty:
        return 0

    grouped = df.groupby(
        ["is_element", "side", "cluster_id"], as_index=False
    ).agg(rep_seq=("representative_sequence", "first"), n_seqs=("read_id", "size"))

    qualifying = grouped[grouped["n_seqs"] >= min_size]
    written = 0
    for row in qualifying.itertuples(index=False):
        header = f"{acc_dir.name}__{row.is_element}__{row.side}__cluster{row.cluster_id}"
        out_fh.write(f">{header}\n{row.rep_seq}\n")
        written += 1
    return written


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--accession-list", help="file with one accession dir path per line")
    grp.add_argument("--cdhit_dir", type=Path, help="process every accession dir found directly under this dir")
    ap.add_argument("--min_size", type=int, default=10, help="minimum cluster size, inclusive (default 10)")
    ap.add_argument("--out", required=True, help="output fasta path")
    args = ap.parse_args()

    if args.accession_list:
        with open(args.accession_list) as fh:
            acc_dirs = [Path(line.strip()) for line in fh if line.strip()]
    else:
        acc_dirs = sorted(d for d in args.cdhit_dir.iterdir() if d.is_dir())

    total_written = 0
    with open(args.out, "w") as out_fh:
        for i, acc_dir in enumerate(acc_dirs, 1):
            try:
                total_written += process_accession(acc_dir, args.min_size, out_fh)
            except Exception as e:
                print(f"ERROR processing {acc_dir}: {e}", file=sys.stderr)
            if i % 2000 == 0:
                print(f"[{i}/{len(acc_dirs)}] accessions processed, {total_written} reps written so far",
                      file=sys.stderr, flush=True)

    print(f"done: {total_written} representative sequences (cluster size >= {args.min_size}) "
          f"written from {len(acc_dirs)} accession dirs", file=sys.stderr)


if __name__ == "__main__":
    main()
