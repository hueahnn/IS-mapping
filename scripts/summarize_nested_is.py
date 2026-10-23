# purpose: combine scripts/nested_is.py per-cluster tables (one per sample, or
# one combined table from a --zip_dir run) into two summaries:
#
#   class_summary.tsv        clusters and samples per class, and the fraction
#                            of all no_hit clusters in each
#   host_target_summary.tsv  one row per (host is_element, target, class) for
#                            the internal/terminal clusters -- the nested-
#                            insertion candidates. Start here: n_samples_paired
#                            (both flanks of a host copy pointing into one site
#                            of the target) is the strongest evidence. Targets
#                            hit by only one side (n_left or n_right ~0) at one
#                            site in many samples are often not nesting -- the
#                            host's database entry stops short and its flank
#                            runs on into a related element.
#
# usage:
#   python summarize_nested_is.py --table_list tables.txt --output_dir nested_is_summary
#   python summarize_nested_is.py --tables nested_is.clusters.tsv.gz --output_dir .

import argparse
import os
import sys

import numpy as np
import pandas as pd

CLASSES = ["internal", "terminal", "self", "distal_only", "no_is_hit", "too_short"]
USECOLS = ["sample", "is_element", "side", "class", "target", "target_db", "target_len",
           "junction_pos", "paired_in_target", "pident"]


def main():
    p = argparse.ArgumentParser(description="Summarize nested_is.py per-cluster tables.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--tables", nargs="+", help="nested_is.py output table(s)")
    src.add_argument("--table_list", help="file listing one table path per line")
    p.add_argument("--output_dir", required=True)
    p.add_argument("--site_tol", type=int, default=20,
                   help="junction positions this close count as one site (n_sites)")
    a = p.parse_args()

    if a.table_list:
        with open(a.table_list) as fh:
            tables = [line.strip() for line in fh if line.strip()]
    else:
        tables = a.tables

    parts = []
    for i, path in enumerate(tables):
        part = pd.read_csv(path, sep="\t", usecols=USECOLS,
                           dtype={"target": str, "target_db": str})
        if not part.empty:  # header-only: a sample with no no_hit clusters
            parts.append(part)
        if (i + 1) % 10000 == 0:
            print(f"  {i + 1}/{len(tables)} tables read", file=sys.stderr, flush=True)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=USECOLS)
    df["paired_in_target"] = df.paired_in_target.astype(str) == "True"
    print(f"{len(tables)} tables, {len(df)} clusters", file=sys.stderr)

    os.makedirs(a.output_dir, exist_ok=True)

    cls = (df.groupby("class").agg(n_clusters=("sample", "size"), n_samples=("sample", "nunique"))
           .reindex(CLASSES).fillna(0).astype(int))
    total = cls.n_clusters.sum()
    cls["frac_clusters"] = (cls.n_clusters / total).round(4) if total else 0.0
    cls.to_csv(os.path.join(a.output_dir, "class_summary.tsv"), sep="\t")

    def n_sites(pos):
        pos = np.sort(pos.to_numpy())
        return int(1 + (np.diff(pos) > a.site_tol).sum()) if len(pos) else 0

    cand = df[df["class"].isin(["internal", "terminal"])]
    paired_sample = cand["sample"].where(cand.paired_in_target)
    summ = (cand.assign(paired_sample=paired_sample,
                        is_left=cand.side == "left", is_right=cand.side == "right")
            .groupby(["is_element", "target", "target_db", "class"])
            .agg(n_clusters=("sample", "size"),
                 n_samples=("sample", "nunique"),
                 n_samples_paired=("paired_sample", "nunique"),
                 n_left=("is_left", "sum"),
                 n_right=("is_right", "sum"),
                 target_len=("target_len", "first"),
                 n_sites=("junction_pos", n_sites),
                 top_site=("junction_pos", lambda s: int(s.mode().iloc[0])),
                 median_pident=("pident", "median"))
            .reset_index().sort_values(["n_samples_paired", "n_samples"], ascending=False))
    summ["target_len"] = summ.target_len.astype(int)
    summ["median_pident"] = summ.median_pident.round(3)
    summ.to_csv(os.path.join(a.output_dir, "host_target_summary.tsv"), sep="\t", index=False)

    print(cls.to_string(), file=sys.stderr)
    print(f"{len(summ)} host->target pairs among internal/terminal clusters", file=sys.stderr)


if __name__ == "__main__":
    main()
