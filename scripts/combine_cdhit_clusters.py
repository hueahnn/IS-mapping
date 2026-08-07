# purpose: combine one sample's CD-HIT centroid (.cdhit) and cluster-membership
# (.cdhit.clstr) files -- produced by the cd_hit rule in Snakefile -- into one TSV
# per IS-element/side pair, readable directly with pandas.read_csv(sep="\t").
# usage:
#   python combine_cdhit_clusters.py /n/scratch/.../cd-hit/DRR033895 --output_dir out/

# Each row is one clustered read (a member of a non-singleton cluster), reporting
# which cluster it fell into and that cluster's representative read + sequence.
# Singleton clusters (a read that didn't match anything else) are dropped entirely,
# since there's no cluster of interest to report for them.
#
# .cdhit fasta headers and .cdhit.clstr member-line headers share the format
# produced by scripts/overhangs.py's extract_overhangs_from_bam:
#   {sample}|{read_id}|{is_element}|pos:{p}|side:{left,right}|strand:{+,-}|pident:{pct}
# The .clstr file marks a cluster's representative sequence with a trailing "*";
# every other member line ends with "at cstart:cend:rstart:rend/strand/pident%" (its
# alignment to the representative). cd-hit-est is run with -d 0 in the cd_hit rule,
# so headers are never truncated in the .clstr file.

import argparse
import csv
import re
from pathlib import Path

import pandas as pd

HEADER_RE = re.compile(
    r"^(?P<sample>[^|]+)\|(?P<read_id>[^|]+)\|(?P<is_element>[^|]+)"
    r"\|pos:(?P<pos>[^|]+)\|side:(?P<side>[^|]+)\|strand:(?P<strand>[+-])"
    r"\|pident:(?P<pident>[\d.]+)$"
)

# a .clstr member line is either the cluster representative (ends "... *") or a
# member aligned to it (ends "... at cstart:cend:rstart:rend/strand/pident%")
MEMBER_RE = re.compile(
    r"^\d+\t\d+nt, >(?P<header>.+)\.\.\. (?:\*|at \d+:\d+:\d+:\d+/[+-]/[\d.]+%)$"
)


def parse_header(header: str) -> dict:
    """Split one overhangs.py-style fasta header into its component fields."""
    m = HEADER_RE.match(header)
    if not m:
        raise ValueError(f"unrecognized fasta header: {header!r}")
    fields = m.groupdict()
    fields["pident"] = float(fields["pident"])
    return fields


def load_sequences(cdhit_fasta_path: Path) -> dict:
    """Read a .cdhit fasta (cluster representatives only) into {header: sequence}."""
    sequences = {}
    header = None
    chunks = []
    with open(cdhit_fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            sequences[header] = "".join(chunks)
    return sequences


def parse_clusters(clstr_path: Path) -> list:
    """
    Parse a .cdhit.clstr file into a list of clusters, each
    {"members": [header, ...], "representative": header}, in cd-hit's own
    0-indexed cluster order (file order).
    """
    clusters = []
    current = None
    with open(clstr_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                if current is not None:
                    clusters.append(current)
                current = {"members": [], "representative": None}
                continue
            m = MEMBER_RE.match(line)
            if not m:
                raise ValueError(f"unrecognized .clstr line: {line!r}")
            header = m.group("header")
            current["members"].append(header)
            if line.endswith("*"):
                current["representative"] = header
        if current is not None:
            clusters.append(current)
    return clusters


def build_reads_table(cdhit_fasta_path: Path, clstr_path: Path) -> pd.DataFrame:
    """
    One row per read in a non-singleton cluster: the read's own parsed fields,
    which cluster it belongs to, and its cluster's representative read + sequence.
    Singleton clusters are dropped.
    """
    representative_sequences = load_sequences(cdhit_fasta_path)
    clusters = parse_clusters(clstr_path)

    rows = []
    for cluster_id, cluster in enumerate(clusters):
        if len(cluster["members"]) < 2:
            continue  # throw out singleton clusters

        rep_header = cluster["representative"]
        rep_fields = parse_header(rep_header)
        rep_sequence = representative_sequences[rep_header]

        for member_header in cluster["members"]:
            fields = parse_header(member_header)
            rows.append({
                "sample": fields["sample"],
                "is_element": fields["is_element"],
                "side": fields["side"],
                "read_id": fields["read_id"],
                "pos": fields["pos"],
                "strand": fields["strand"],
                "pident": fields["pident"],
                "cluster_id": cluster_id,
                "representative_read_id": rep_fields["read_id"],
                "representative_sequence": rep_sequence,
            })

    return pd.DataFrame(rows, columns=[
        "sample", "is_element", "side", "read_id", "pos", "strand", "pident",
        "cluster_id", "representative_read_id", "representative_sequence",
    ])


def process_sample(sample_dir: Path, output_dir: Path) -> None:
    """
    Read {sample}.cdhit_manifest.tsv in sample_dir and write one combined
    {sample}.reads.tsv covering every is_element/side pair (one row per read,
    same as before -- cluster_id is only unique within a given
    (is_element, side), so group by all three together, not cluster_id alone).
    """
    manifest_path = next(sample_dir.glob("*.cdhit_manifest.tsv"))
    sample = manifest_path.name.removesuffix(".cdhit_manifest.tsv")
    output_dir.mkdir(parents=True, exist_ok=True)

    dfs = []
    with open(manifest_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cdhit_fasta_path = Path(row["centroids_path"])
            clstr_path = Path(row["clstr_path"])

            df = build_reads_table(cdhit_fasta_path, clstr_path)
            dfs.append(df)
            print(f"{row['is_element']} {row['side']}: "
                  f"{len(df)} reads across {df['cluster_id'].nunique()} clusters")

    out_path = output_dir / f"{sample}.reads.tsv"
    combined = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(
        columns=["sample", "is_element", "side", "read_id", "pos", "strand", "pident",
                 "cluster_id", "representative_read_id", "representative_sequence"])
    combined.to_csv(out_path, sep="\t", index=False)
    print(f"-> {out_path} ({len(combined)} total reads)")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sample_dir",
                    help="one sample's cd-hit output dir, e.g. cd-hit/DRR033895 "
                         "(must contain {sample}.cdhit_manifest.tsv)")
    p.add_argument("--output_dir", default=None,
                    help="where to write {sample}__{is_element}__{side}.reads.tsv "
                         "(default: same as sample_dir)")
    args = p.parse_args()

    sample_dir = Path(args.sample_dir)
    output_dir = Path(args.output_dir) if args.output_dir else sample_dir
    process_sample(sample_dir, output_dir)
