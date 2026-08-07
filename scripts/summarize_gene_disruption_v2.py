# purpose: stream every {accession}.zip in gene_disruption_v2/ (one zip per
# accession, one {is_element}.tsv per IS element that accession had clusters
# for -- see blast_clusters_to_ref.py's write_split_tables), and tally, across
# EVERY accession and EVERY IS element, how many clusters with >= min_size
# member reads (n_seqs) fall into each hit_type category, plus how many
# protein_coding clusters have a paired left/right confirmation (the
# "pairing" column, already computed -- see pairing/add_pairing_column.py).
#
# Memory-safe by construction: only one (accession, is_element) TSV's rows
# are ever held in memory at a time (grouped by (side, cluster_id) into a
# tiny per-cluster dict, then discarded), and only a handful of global
# integer counters persist across the whole scan -- same lesson as
# scripts/cluster_and_align_by_is.py's OOM fix, applied from the start here
# rather than rediscovered the hard way.
#
# usage:
#   python summarize_gene_disruption_v2.py --min_size 10

import argparse
import csv
import io
import sys
import zipfile
from pathlib import Path

DEFAULT_DIR = Path("/n/scratch/users/h/hua575/atb_filtered/gene_disruption_v2")

# a cluster can have multiple hit rows (one per ref_genome it hit, or >1
# overlapping gene at the same locus via gene_rank) -- collapse to one
# category per cluster so counts partition cleanly. Priority: any
# protein_coding row -> protein_coding; else any intergenic -> intergenic;
# else no_hit (already guaranteed unique per cluster -- assemble_table only
# emits a no_hit row when there were no qualifying hits at all).
PRIORITY = ["protein_coding", "intergenic", "no_hit"]


def process_zip(zip_path: Path, min_size: int, counts: dict) -> None:
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if not name.endswith(".tsv"):
                continue
            with zf.open(name) as fh:
                reader = csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t")
                groups = {}  # (side, cluster_id) -> {"n_seqs", "hit_types", "pc_paired"}
                for row in reader:
                    key = (row["side"], row["cluster_id"])
                    g = groups.get(key)
                    if g is None:
                        g = {"n_seqs": int(row["n_seqs"]), "hit_types": set(), "pc_paired": False}
                        groups[key] = g
                    ht = row["hit_type"]
                    g["hit_types"].add(ht)
                    if ht == "protein_coding" and row["pairing"] != "none":
                        g["pc_paired"] = True

                for g in groups.values():
                    if g["n_seqs"] < min_size:
                        continue
                    for cat in PRIORITY:
                        if cat in g["hit_types"]:
                            counts[cat] += 1
                            if cat == "protein_coding" and g["pc_paired"]:
                                counts["protein_coding_paired"] += 1
                            break


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene_disruption_dir", type=Path, default=DEFAULT_DIR)
    p.add_argument("--min_size", type=int, default=10,
                    help="drop clusters with fewer than this many member reads (n_seqs); "
                         "default 10 means 'do not count if < 10, count if >= 10'")
    args = p.parse_args()

    counts = {"protein_coding": 0, "intergenic": 0, "no_hit": 0, "protein_coding_paired": 0}
    n_ok, n_err = 0, 0

    zip_paths = sorted(args.gene_disruption_dir.glob("*.zip"))
    total = len(zip_paths)
    print(f"{total} accession zips found in {args.gene_disruption_dir}", file=sys.stderr)

    for i, zp in enumerate(zip_paths, 1):
        try:
            process_zip(zp, args.min_size, counts)
            n_ok += 1
        except Exception as e:
            n_err += 1
            print(f"ERROR {zp}: {e}", file=sys.stderr)

        if i % 5000 == 0:
            print(f"[{i}/{total}] accessions processed -- running counts: {counts}",
                  file=sys.stderr, flush=True)

    total_clusters = counts["protein_coding"] + counts["intergenic"] + counts["no_hit"]
    print(f"\naccessions processed: {n_ok} ok, {n_err} errors")
    print(f"min_size (n_seqs) filter: >= {args.min_size}")
    print(f"total clusters counted: {total_clusters:,}")
    print(f"  protein_coding: {counts['protein_coding']:,}")
    print(f"  intergenic:     {counts['intergenic']:,}")
    print(f"  no_hit:         {counts['no_hit']:,}")
    print(f"protein_coding clusters that are paired: {counts['protein_coding_paired']:,}")


if __name__ == "__main__":
    main()
