#!/usr/bin/env python3
"""
Count insertion-sequence (IS) insertion events per IS element across all
sequenced E. coli genomes (SRA accessions), using the gene_disruption/*.zip
output of /home/hua575/baymlab/mapping's blast_clusters_to_ref.py +
pairing/add_pairing_column.py pipeline.

Definition of "one insertion event": within a given (sample, is_element),
each (side, cluster_id) is a node. A row's `pairing` value (when not "none")
links that row's (side, cluster_id) to the opposite side's cluster with that
id -- both sides are two pieces of evidence for the SAME physical
transposition junction. Connected components over that graph = distinct
insertion events for that (sample, is_element). This avoids double-counting
a single event once per side, and avoids double-counting across the 6
reference genomes a cluster's flanking sequence happens to BLAST-hit.

Usage:
    python count_is_events.py --zip-dir DIR --out-csv OUT.csv [--workers N] [--limit N]
"""
import argparse
import csv
import os
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
from multiprocessing import Pool


def union_find_find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union_find_union(parent, a, b):
    ra, rb = union_find_find(parent, a), union_find_find(parent, b)
    if ra != rb:
        parent[ra] = rb


def count_events_in_tsv_text(text):
    """
    Given the raw text of one {sample}/{is_element}.tsv, return the number of
    distinct insertion events (connected components over (side, cluster_id)
    nodes linked by the pairing column).
    """
    lines = text.splitlines()
    if not lines:
        return 0
    header = lines[0].split("\t")
    try:
        side_i = header.index("side")
        cid_i = header.index("cluster_id")
        pairing_i = header.index("pairing")
    except ValueError:
        return 0

    nodes = set()
    parent = {}
    edges = []

    for line in lines[1:]:
        if not line:
            continue
        fields = line.split("\t")
        side = fields[side_i]
        cluster_id = fields[cid_i]
        pairing = fields[pairing_i]

        node = (side, cluster_id)
        if node not in nodes:
            nodes.add(node)
            parent[node] = node

        if pairing and pairing != "none":
            other_side = "right" if side == "left" else "left"
            other_node = (other_side, pairing)
            if other_node not in nodes:
                nodes.add(other_node)
                parent[other_node] = other_node
            edges.append((node, other_node))

    for a, b in edges:
        union_find_union(parent, a, b)

    roots = {union_find_find(parent, n) for n in nodes}
    return len(roots)


def process_zip(zip_path):
    """
    Returns dict: is_element -> (n_events_in_this_sample) for one accession zip.
    """
    counts = {}
    try:
        with zipfile.ZipFile(zip_path) as zf:
            for name in zf.namelist():
                if not name.endswith(".tsv"):
                    continue
                is_element = Path(name).stem
                text = zf.read(name).decode("utf-8", errors="replace")
                n_events = count_events_in_tsv_text(text)
                if n_events:
                    counts[is_element] = counts.get(is_element, 0) + n_events
    except (zipfile.BadZipFile, OSError) as e:
        print(f"WARN: failed to read {zip_path}: {e}", file=sys.stderr)
    return Path(zip_path).stem, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip-dir", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--per-sample-out", default=None,
                     help="optional: write a long-format (sample, is_element, n_events) CSV too")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None, help="only process first N zips (debug)")
    args = ap.parse_args()

    zip_dir = Path(args.zip_dir)
    zip_paths = []
    with os.scandir(zip_dir) as it:
        for entry in it:
            if entry.name.endswith(".zip"):
                zip_paths.append(entry.path)
    zip_paths.sort()
    if args.limit:
        zip_paths = zip_paths[: args.limit]

    print(f"Processing {len(zip_paths)} zip files with {args.workers} workers...", file=sys.stderr)

    total_events = defaultdict(int)
    samples_with_element = defaultdict(int)

    per_sample_fh = None
    per_sample_writer = None
    if args.per_sample_out:
        per_sample_fh = open(args.per_sample_out, "w", newline="")
        per_sample_writer = csv.writer(per_sample_fh)
        per_sample_writer.writerow(["sample", "is_element", "n_events"])

    n_done = 0
    with Pool(args.workers) as pool:
        for sample, counts in pool.imap_unordered(process_zip, zip_paths, chunksize=50):
            for is_element, n_events in counts.items():
                total_events[is_element] += n_events
                samples_with_element[is_element] += 1
                if per_sample_writer:
                    per_sample_writer.writerow([sample, is_element, n_events])
            n_done += 1
            if n_done % 5000 == 0:
                print(f"  {n_done}/{len(zip_paths)} zips processed...", file=sys.stderr)

    if per_sample_fh:
        per_sample_fh.close()

    with open(args.out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["is_element", "n_insertion_events", "n_samples_with_element"])
        for is_element in sorted(total_events, key=lambda k: -total_events[k]):
            w.writerow([is_element, total_events[is_element], samples_with_element[is_element]])

    print(f"Done. {len(total_events)} distinct IS elements had >=1 insertion event.", file=sys.stderr)
    print(f"Total insertion events across all elements: {sum(total_events.values())}", file=sys.stderr)
    print(f"Wrote {args.out_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
