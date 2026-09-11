#!/usr/bin/env python3
"""Add pairing columns to gene_disruption TSVs.

Reads are mapped against a MASKED reference genome (the IS sequence itself is
not present at the insertion site), so a true left/right overhang pair should
land back-to-back at the junction, with only a small gap (or overlap, e.g.
from a target-site duplication) between them.

TWO-TIER PAIRING. Every IS-element table in a genome's zip is loaded together
and paired in one pass, because a single real insertion does not always get
reported under a single IS name: the ISfinder database is collapsed at 99%
identity but still holds many near-identical elements, so one insertion's left
flank can align to one DB entry and its right flank to another. (ISMapper shows
the same thing -- in ismapper/ERR9860257 the ISEal1 and ISEc16 tables both
report an identical insertion at 384559-384562.) Pairing therefore runs in two
tiers:

  tier 1 (same_is)  -- left and right clusters of the SAME is_element
  tier 2 (cross_is) -- whatever is still unpaired, against clusters of ANY
                       is_element

Tier 1 is resolved completely before any tier-2 pair is considered, so a
same-IS partner always wins over a cross-IS one even when the cross-IS partner
sits closer. `pairing_type` records which tier produced each pair.

Both tiers apply the same physical criteria: same (ref_genome, contig), same
hit_strand, and abs(right.junction_pos - left.junction_pos) < max_gap. The
SIGN of that difference is strand-dependent -- for plus-strand pairs the right
junction sits below the left, for minus-strand pairs above -- because
blast_clusters_to_ref.py takes `send` for a left overhang and `sstart` for a
right one, and minus-strand HSPs have sstart > send. abs() is what lets one
window cover both; `pairing_gap` reports the signed value.

Within a tier, the globally nearest pair is assigned first, then the next
nearest among clusters still unassigned, and so on. Matching is one-to-one at
the CLUSTER level, not the row level: a cluster that BLASTs to several
reference genomes is paired independently in each (the flanking locus need not
exist, or be unique, across all six assemblies), but its duplicate rows within
one genome -- e.g. the two gene_rank rows emitted when a junction lands in
overlapping CDSs -- are one unit and always receive the same pairing verdict.

Rows with no position data (no_hit rows) and any row whose side is neither
left nor right are skipped, and keep pairing_type = "none".

Columns added after cluster_id:
  pairing        -- partner as "{is_element}:{side}:{cluster_id}", or "none".
                    Fully qualified because a bare cluster_id is ambiguous
                    once cross-IS pairs are possible.
  pairing_type   -- "same_is", "cross_is", or "none"
  pairing_gap    -- signed right.junction_pos - left.junction_pos, in bp

DEFERRED -- repeated reference regions. ISMapper's Figure 2 handles an
insertion next to a multi-copy repeat by anchoring on the flank that passed its
depth cutoff and re-searching for a sub-threshold partner near it
(closestBed of the depth-FILTERED flank against the UNFILTERED opposite side --
see create_bed_files/check_unpaired_hits in the ISMapper source), marking every
hit rescued that way with a "?" confidence suffix. That failure mode barely
applies here: ISMapper dilutes depth because BWA scatters multi-mapping reads
across repeat copies, whereas this pipeline clusters overhangs by SEQUENCE
before placing them, so a cluster stays intact and full-depth however
repetitive the locus. What remains is placement ambiguity -- one healthy
cluster with several equally-good loci -- and the Figure 2 remedy does adapt to
it: let the confident partner pick which repeat copy the ambiguous flank goes
to, and flag that the pairing, not the alignment, made the choice. That is not
implemented here because it cannot be: blast_clusters_to_ref.py's
best_hit_per_genome() keeps a single best-bitscore hit per (cluster,
ref_genome), so the alternative repeat copies are already gone by the time this
script runs. See the first item in README.md's "To Dos for the Future".
"""
import argparse
import csv
import io
import os
import sys
import zipfile

DEFAULT_MAX_GAP = 50

# inserted after cluster_id, in this order
NEW_COLUMNS = ["pairing", "pairing_type", "pairing_gap"]


def parse_float(s):
    return float(s) if s else None


class Placement:
    """One cluster at one reference locus, plus every row that reports it."""

    __slots__ = ("is_element", "side", "cluster_id", "ref_genome", "contig",
                 "junction", "strand", "row_idxs")

    def __init__(self, key, row_idxs):
        (self.is_element, self.side, self.cluster_id, self.ref_genome,
         self.contig, self.junction, self.strand) = key
        self.row_idxs = row_idxs

    @property
    def cluster_key(self):
        """The unit that may be paired at most once per reference genome."""
        return (self.is_element, self.side, self.cluster_id)

    @property
    def label(self):
        return f"{self.is_element}:{self.side}:{self.cluster_id}"


def build_placements(rows):
    """Collapse rows into one Placement per (cluster, locus).

    Rows that cannot be paired at all -- no_hit rows, and anything whose side
    is not left/right -- are simply left out, so they keep the "none" defaults.
    """
    grouped = {}
    for idx, r in enumerate(rows):
        junction = parse_float(r.get("junction_pos"))
        if (not r.get("ref_genome") or not r.get("contig")
                or junction is None or not r.get("hit_strand")):
            continue
        if r.get("side") not in ("left", "right"):
            continue
        key = (r.get("is_element"), r.get("side"), r.get("cluster_id"),
               r.get("ref_genome"), r.get("contig"), junction,
               r.get("hit_strand"))
        grouped.setdefault(key, []).append(idx)
    return [Placement(k, idxs) for k, idxs in grouped.items()]


def format_gap(gap):
    return str(int(gap)) if float(gap).is_integer() else repr(gap)


def compute_pairing(rows, max_gap):
    """Mutate rows in place, setting the NEW_COLUMNS keys on every row dict."""
    for r in rows:
        r["pairing"] = "none"
        r["pairing_type"] = "none"
        r["pairing_gap"] = ""

    by_genome = {}
    for p in build_placements(rows):
        by_genome.setdefault(p.ref_genome, []).append(p)

    # each reference genome is paired independently -- a cluster's flanking
    # locus need not exist, or be unique, across all six assemblies
    for placements in by_genome.values():
        lefts = [p for p in placements if p.side == "left"]
        rights = [p for p in placements if p.side == "right"]

        candidates = []
        for left in lefts:
            for right in rights:
                if left.contig != right.contig or left.strand != right.strand:
                    continue
                gap = right.junction - left.junction
                if abs(gap) >= max_gap:
                    continue
                tier = 0 if left.is_element == right.is_element else 1
                candidates.append((tier, abs(gap), gap, left, right))

        # tier first, so every same-IS pair is settled before the first
        # cross-IS pair is even considered; nearest end next; then the labels,
        # purely so ties break the same way regardless of the order the zip
        # happened to yield its tables (with 37 IS4 pairs all at gap 11 in a
        # single sample, ties are the common case, not the rare one).
        candidates.sort(key=lambda c: (c[0], c[1], c[3].label, c[4].label))

        assigned = set()
        for tier, _, gap, left, right in candidates:
            if left.cluster_key in assigned or right.cluster_key in assigned:
                continue
            assigned.add(left.cluster_key)
            assigned.add(right.cluster_key)

            pairing_type = "same_is" if tier == 0 else "cross_is"
            gap_str = format_gap(gap)
            for partner, placement in ((right, left), (left, right)):
                for idx in placement.row_idxs:
                    rows[idx]["pairing"] = partner.label
                    rows[idx]["pairing_type"] = pairing_type
                    rows[idx]["pairing_gap"] = gap_str


def process_zip(zip_path, max_gap, out_dir):
    accession = os.path.basename(zip_path)
    out_path = os.path.join(out_dir, accession)

    tables = {}
    all_rows = []
    with zipfile.ZipFile(zip_path, "r") as zin:
        for name in [n for n in zin.namelist() if n.endswith(".tsv")]:
            reader = csv.DictReader(io.StringIO(zin.read(name).decode()),
                                    delimiter="\t")
            if not reader.fieldnames:
                continue
            fieldnames = [f for f in reader.fieldnames if f not in NEW_COLUMNS]
            at = fieldnames.index("cluster_id") + 1
            fieldnames[at:at] = NEW_COLUMNS
            rows = list(reader)
            tables[name] = (fieldnames, rows)
            # same dict objects as in `tables`, so compute_pairing's in-place
            # edits are visible when the tables are written back out
            all_rows.extend(rows)

    # every IS element at once -- cross-IS pairing spans the whole zip
    compute_pairing(all_rows, max_gap)

    out_buf = io.BytesIO()
    with zipfile.ZipFile(out_buf, "w", zipfile.ZIP_DEFLATED) as zout:
        for name, (fieldnames, rows) in tables.items():
            out_io = io.StringIO()
            writer = csv.DictWriter(out_io, fieldnames=fieldnames,
                                    delimiter="\t", lineterminator="\n")
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
