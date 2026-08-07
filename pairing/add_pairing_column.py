#!/usr/bin/env python3
"""Add a 'pairing' column to gene_disruption TSVs.

Reads are mapped against a MASKED reference genome (the IS sequence itself is
not present at the insertion site), so a true left/right overhang pair should
land back-to-back at the junction, with only a small gap (or overlap, e.g.
from a target-site duplication) between them.

For every IS-element TSV inside a genome's zip, rows are grouped by
(ref_genome, contig). Within each group, a left-side overhang and a
right-side overhang on the SAME hit_strand are a valid candidate pair if
abs(right.junction_pos - left.junction_pos) < max_gap (overlap allowed,
i.e. the difference may be negative).

Matching is one-to-one: among all candidate (left, right) pairs in a group,
the globally smallest gap is assigned first, then the next smallest among
remaining unassigned rows, and so on. Rows with no assigned partner (or no
position data, e.g. no_hit rows) get pairing = "none". A matched row's
pairing value is the cluster_id of the row it was paired with.
"""
import argparse
import csv
import io
import os
import sys
import zipfile

DEFAULT_MAX_GAP = 50


def parse_float(s):
    return float(s) if s else None


def compute_pairing(rows, max_gap):
    """Mutate rows in place, setting a 'pairing' key on every row dict."""
    for r in rows:
        r["pairing"] = "none"

    groups = {}
    for idx, r in enumerate(rows):
        ref_genome = r.get("ref_genome")
        contig = r.get("contig")
        junction = parse_float(r.get("junction_pos"))
        strand = r.get("hit_strand")
        side = r.get("side")
        if not ref_genome or not contig or junction is None or not strand:
            continue
        if side not in ("left", "right"):
            continue
        key = (ref_genome, contig)
        groups.setdefault(key, {"left": [], "right": []})
        groups[key][side].append((idx, junction, strand))

    candidates = []
    for sides in groups.values():
        for li, ljunction, lstrand in sides["left"]:
            for ri, rjunction, rstrand in sides["right"]:
                if rstrand != lstrand:
                    continue
                gap = abs(rjunction - ljunction)
                if gap < max_gap:
                    candidates.append((gap, li, ri))

    candidates.sort(key=lambda x: x[0])
    assigned_left = set()
    assigned_right = set()
    for gap, li, ri in candidates:
        if li in assigned_left or ri in assigned_right:
            continue
        assigned_left.add(li)
        assigned_right.add(ri)
        rows[li]["pairing"] = rows[ri]["cluster_id"]
        rows[ri]["pairing"] = rows[li]["cluster_id"]


def process_zip(zip_path, max_gap, out_dir):
    accession = os.path.basename(zip_path)
    out_path = os.path.join(out_dir, accession)

    with zipfile.ZipFile(zip_path, "r") as zin:
        names = [n for n in zin.namelist() if n.endswith(".tsv")]
        out_buf = io.BytesIO()
        with zipfile.ZipFile(out_buf, "w", zipfile.ZIP_DEFLATED) as zout:
            for name in names:
                raw = zin.read(name).decode()
                reader = csv.DictReader(io.StringIO(raw), delimiter="\t")
                fieldnames = list(reader.fieldnames)
                if "pairing" in fieldnames:
                    fieldnames.remove("pairing")
                fieldnames.insert(fieldnames.index("cluster_id") + 1, "pairing")
                rows = list(reader)

                compute_pairing(rows, max_gap)

                out_io = io.StringIO()
                writer = csv.DictWriter(out_io, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
                zout.writestr(name, out_io.getvalue())

    with open(out_path, "wb") as f:
        f.write(out_buf.getvalue())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--zip", help="single genome zip to process")
    grp.add_argument("--zip-list", help="file with one zip path per line")
    ap.add_argument("--max-gap", type=float, default=DEFAULT_MAX_GAP,
                     help=f"max allowed |right_junction - left_junction| in bp (default {DEFAULT_MAX_GAP})")
    ap.add_argument("--out-dir", required=True, help="output directory for augmented zips")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if args.zip:
        zip_paths = [args.zip]
    else:
        with open(args.zip_list) as fh:
            zip_paths = [line.strip() for line in fh if line.strip()]

    for zip_path in zip_paths:
        try:
            process_zip(zip_path, args.max_gap, args.out_dir)
        except Exception as e:
            print(f"ERROR processing {zip_path}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
