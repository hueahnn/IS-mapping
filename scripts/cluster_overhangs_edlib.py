# purpose: cluster overhang reads within each (is_element, side) group using a
# position-anchored edit distance instead of CD-HIT's percent-identity/alignment
# approach -- replaces both the cd-hit-est loop and combine_cdhit_clusters.py's
# combine step in one script, producing {sample}.reads.tsv directly.
#
# ANCHORING: overhangs.py writes side=="left" overhangs as the read's LEADING
# soft-clip (seq[:left_clip]) -- the junction-adjacent base sits at the END of the
# string; index 0 is the distal end, which drifts read-to-read purely because
# soft-clip length varies with read length. side=="right" overhangs are the
# opposite (seq[len(seq)-right_clip:]) -- the junction-adjacent base is already at
# index 0. For two overhangs to be comparable over a shared genomic window
# regardless of length, comparison must anchor at the junction end for BOTH
# sides -- so left-side sequences are reversed before comparison (and reversed
# back for output: representative_sequence always stores the original forward
# orientation, matching combine_cdhit_clusters.py's existing convention).
#
# ALGORITHM: within each (is_element, side) group, sort sequences by
# (comparison-orientation) length descending. Take the longest unclustered
# sequence as a new cluster's representative/seed. Compare every other
# unclustered sequence to it with edlib's SHW ("prefix") alignment mode --
# query=candidate (always <= seed's length, by sort order), target=seed --
# which computes edit distance anchored at the shared start (the junction end,
# post-reversal) with a free gap at the END of target only, so the seed's extra
# tail beyond the candidate's length is never penalized. This is the "anchored,
# not length-penalized" comparison the CD-HIT replacement is meant to provide --
# unlike global Levenshtein (or the tool `starcode`), a 200bp vs 20bp pair
# doesn't cost >=180 just from the length gap. Any candidate within
# max_edit_frac * len(candidate) edits joins the cluster and is removed from the
# pool; repeat with the next-longest remaining sequence until the pool is empty.
#
# MERGE PASS: greedy clustering is path-dependent -- a true single junction can
# fragment into more than one cluster depending on which read happened to seed
# first. After the primary pass, the SAME algorithm re-runs over just the
# resulting cluster representatives (far fewer sequences) with a looser
# max-edit-fraction (--merge_edit_frac), merging clusters whose representatives
# still fall within that looser threshold. --min_cluster_size is applied AFTER
# the merge pass, so merging gets a chance to rescue clusters that would
# otherwise look like noise (mirrors combine_cdhit_clusters.py's singleton-drop
# convention, generalized to a configurable minimum).
#
# usage:
#   python cluster_overhangs_edlib.py /path/to/{sample}.manifest.tsv \
#       --sample_id SRR12705059 --output_dir out/ \
#       --max_edit_frac 0.10 --merge_edit_frac 0.15 --min_cluster_size 2

import argparse
import csv
import re
from pathlib import Path

import edlib
import pandas as pd

# same shape as combine_cdhit_clusters.py's HEADER_RE, but tolerant of the
# trailing "|mapq:{int}" field the current overhangs.py appends (that regex
# is anchored with $ right after pident and would fail to match real headers)
HEADER_RE = re.compile(
    r"^(?P<sample>[^|]+)\|(?P<read_id>[^|]+)\|(?P<is_element>[^|]+)"
    r"\|pos:(?P<pos>[^|]+)\|side:(?P<side>[^|]+)\|strand:(?P<strand>[+-])"
    r"\|pident:(?P<pident>[\d.]+)(?:\|mapq:\d+)?$"
)


def parse_header(header: str) -> dict:
    """Split one overhangs.py-style fasta header into its component fields."""
    m = HEADER_RE.match(header)
    if not m:
        raise ValueError(f"unrecognized fasta header: {header!r}")
    fields = m.groupdict()
    fields["pident"] = float(fields["pident"])
    return fields


def read_fasta(fasta_path: Path) -> list[dict]:
    """Parse one manifest-listed overhang FASTA into a list of per-read dicts
    (header fields + "seq" in original, forward orientation)."""
    records = []
    header = None
    seq_chunks: list[str] = []

    def flush():
        if header is None:
            return
        fields = parse_header(header)
        fields["seq"] = "".join(seq_chunks)
        records.append(fields)

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        flush()
    return records


def edit_threshold(query_len: int, max_edit_frac: float) -> int:
    """Max allowed edits for a query of this length -- proportional rather than
    a fixed count, since overhangs span ~5bp-220bp+ in this pipeline (overhangs.py's
    min_clip_len floor is 5bp, no ceiling) and a fixed count would be nonsensically
    strict on long overhangs / meaninglessly loose on short ones."""
    return max(1, round(query_len * max_edit_frac))


def greedy_cluster(records: list[dict], max_edit_frac: float) -> list[list[dict]]:
    """
    records: each must carry a "cmp_seq" key -- the comparison-orientation
    sequence (already reversed for side=="left" by the caller).

    Returns a list of clusters, each a list of member record dicts with the
    cluster's representative/seed first (clusters[i][0]).
    """
    pool = sorted(records, key=lambda r: (-len(r["cmp_seq"]), r["read_id"]))
    clusters = []
    while pool:
        seed = pool.pop(0)  # longest remaining unclustered = new representative
        members = [seed]
        still = []
        for cand in pool:  # cand["cmp_seq"] is always <= len(seed's), by sort order
            k = edit_threshold(len(cand["cmp_seq"]), max_edit_frac)
            result = edlib.align(cand["cmp_seq"], seed["cmp_seq"],
                                  mode="SHW", task="distance", k=k)
            (members if result["editDistance"] != -1 else still).append(cand)
        clusters.append(members)
        pool = still
    return clusters


def merge_clusters(clusters: list[list[dict]], merge_edit_frac: float) -> list[list[dict]]:
    """Second pass: re-cluster just the representatives (clusters[i][0]) with a
    looser threshold, to rescue path-dependent fragmentation from the greedy
    primary pass -- representatives are far fewer than raw reads, so this stays
    cheap even though it's the same O(n^2)-shaped algorithm."""
    if len(clusters) <= 1:
        return clusters

    reps = [c[0] for c in clusters]
    rep_clusters = greedy_cluster(reps, merge_edit_frac)

    rep_id_to_cluster_idx = {id(c[0]): i for i, c in enumerate(clusters)}

    merged = []
    for rep_group in rep_clusters:
        surviving_rep = rep_group[0]  # longest in this merged group (greedy_cluster's own sort)
        combined_members = []
        for rep_record in rep_group:
            orig_idx = rep_id_to_cluster_idx[id(rep_record)]
            combined_members.extend(clusters[orig_idx])
        combined_members.sort(key=lambda r: r is not surviving_rep)  # surviving rep first
        merged.append(combined_members)
    return merged


def cluster_sample(manifest_path: Path, output_dir: Path, sample_id: str,
                    max_edit_frac: float, merge_edit_frac: float,
                    min_cluster_size: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    reads_rows = []
    manifest_rows = []

    with open(manifest_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            is_element, side = row["is_element"], row["side"]
            fasta_path = row["fasta_path"]
            if not fasta_path:
                continue

            records = read_fasta(Path(fasta_path))
            for r in records:
                r["cmp_seq"] = r["seq"][::-1] if side == "left" else r["seq"]

            clusters = greedy_cluster(records, max_edit_frac)
            clusters = merge_clusters(clusters, merge_edit_frac)
            clusters = [c for c in clusters if len(c) >= min_cluster_size]

            for cluster_id, members in enumerate(clusters):
                seed = members[0]
                for m in members:
                    reads_rows.append({
                        "sample": sample_id,
                        "is_element": is_element,
                        "side": side,
                        "read_id": m["read_id"],
                        "pos": m["pos"],
                        "strand": m["strand"],
                        "pident": m["pident"],
                        "cluster_id": cluster_id,
                        "representative_read_id": seed["read_id"],
                        "representative_sequence": seed["seq"],  # forward orientation
                    })

            manifest_rows.append({
                "is_element": is_element, "side": side,
                "n_in": len(records), "n_clusters": len(clusters),
            })
            print(f"{is_element} {side}: {len(records)} reads -> {len(clusters)} clusters "
                  f"(min size {min_cluster_size})")

    reads_df = pd.DataFrame(reads_rows, columns=[
        "sample", "is_element", "side", "read_id", "pos", "strand", "pident",
        "cluster_id", "representative_read_id", "representative_sequence",
    ])
    reads_path = output_dir / f"{sample_id}.reads.tsv"
    reads_df.to_csv(reads_path, sep="\t", index=False)
    print(f"-> {reads_path} ({len(reads_df)} total reads)")

    manifest_df = pd.DataFrame(manifest_rows, columns=["is_element", "side", "n_in", "n_clusters"])
    manifest_out_path = output_dir / f"{sample_id}.cluster_manifest.tsv"
    manifest_df.to_csv(manifest_out_path, sep="\t", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest_path",
                    help="the overhang-extraction manifest, e.g. "
                         "overhangs/{sample}/{sample}.manifest.tsv (overhangs.py output)")
    p.add_argument("--sample_id", required=True)
    p.add_argument("--output_dir", required=True,
                    help="where to write {sample}.reads.tsv and {sample}.cluster_manifest.tsv")
    p.add_argument("--max_edit_frac", type=float, default=0.10,
                    help="max allowed edits as a fraction of the shorter (query) "
                         "sequence's length, for the primary clustering pass")
    p.add_argument("--merge_edit_frac", type=float, default=0.15,
                    help="looser max-edit-fraction used for the representative "
                         "merge pass (rescues path-dependent fragmentation from "
                         "the greedy primary pass)")
    p.add_argument("--min_cluster_size", type=int, default=2,
                    help="drop clusters with fewer than this many member reads, "
                         "applied AFTER the merge pass")
    args = p.parse_args()

    cluster_sample(
        Path(args.manifest_path), Path(args.output_dir), args.sample_id,
        args.max_edit_frac, args.merge_edit_frac, args.min_cluster_size,
    )
