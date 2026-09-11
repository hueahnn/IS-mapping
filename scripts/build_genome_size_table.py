#!/usr/bin/env python3
# purpose: flatten atb/ecoli_atb_filtered.csv (~50 MB, 23 columns) down to the
# two columns blast_clusters_to_ref.py's coverage-fraction join actually needs,
# so that join reads a ~3 MB table instead of re-parsing the full ATB metadata
# once per sample. Same rationale as pairing/build_is_length_table.sh.
#
# Genome_Size is the ATB assembly's total length for that run, which is a far
# better denominator for estimated coverage than a fixed E. coli constant:
# across the 148,756 accessions in atb/ecoli_atb_sra_accessions.txt these span
# 2.97-6.36 Mb (median 5.23 Mb), so assuming MG1655's 4,641,652 bp would be off
# by a median of 11% and by more than 10% for well over half the cohort.
#
# usage:
#   python build_genome_size_table.py \
#       --atb_csv atb/ecoli_atb_filtered.csv \
#       --out atb/genome_size_table.tsv

import argparse
import csv
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--atb_csv", default="atb/ecoli_atb_filtered.csv")
    p.add_argument("--out", default="atb/genome_size_table.tsv")
    args = p.parse_args()

    n_in = n_out = 0
    seen = set()
    with open(args.atb_csv, newline="") as fh, open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["run_accession", "genome_size"])
        for row in csv.DictReader(fh):
            n_in += 1
            run = (row.get("run_accession") or "").strip()
            raw = (row.get("Genome_Size") or "").strip()
            if not run or not raw:
                continue
            try:
                size = int(float(raw))
            except ValueError:
                continue
            if size <= 0 or run in seen:
                continue
            seen.add(run)
            writer.writerow([run, size])
            n_out += 1

    print(f"read {n_in} ATB rows -> wrote {n_out} run_accession/genome_size pairs to {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
