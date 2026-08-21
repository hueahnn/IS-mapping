# purpose: for a list of genome accessions, make one scatter plot per IS element
# seen in that set, comparing each genome's deduplicated left vs. right overhang
# cluster count (one point per genome). Variant of the aggregate-across-genomes
# dot plot in scripts/reads.ipynb (cells 39-43), broken out per genome
# instead of summed across all of them.
# usage:
#   python plot_left_right_overhangs_per_is.py atb/ecoli_atb_sra_accessions_random10.txt

import argparse
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CLUSTER_ROOT = Path("/n/scratch/users/h/hua575/atb_filtered/clusters")
DEFAULT_OUTPUT_DIR = Path("/home/hua575/baymlab/mapping/atb/downstream_analysis/left_right_overhang_plots/all")


def read_accessions(path: Path) -> list:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def collect_counts(accessions: list, cluster_root: Path) -> pd.DataFrame:
    """One row per (accession, is_element, side, n_clusters).

    n_clusters comes straight from each accession's {accession}.cluster_manifest.tsv
    -- cluster_overhangs_edlib.py already writes this column directly (no need
    to open per-group fasta files/archives and count representative records
    the way this used to work against CD-HIT's .cdhit centroid fastas, which
    the current clustering method doesn't produce)."""
    rows = []
    for accession in accessions:
        manifest_path = cluster_root / accession / f"{accession}.cluster_manifest.tsv"
        if not manifest_path.exists():
            print(f"WARNING: no cluster_manifest.tsv for {accession}, skipping")
            continue
        manifest = pd.read_csv(manifest_path, sep="\t")
        manifest["accession"] = accession
        rows.append(manifest[["accession", "is_element", "side", "n_clusters"]])

    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["accession", "is_element", "side", "n_clusters"])


def plot_is_element(is_element: str, plot_df: pd.DataFrame, n_total_accessions: int, output_dir: Path) -> None:
    """One dot plot for a single IS element: one point per genome, right vs. left."""
    max_val = max(plot_df["right"].max(), plot_df["left"].max())
    lims = [0.8, max_val * 1.2]

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(
        plot_df["right"], plot_df["left"],
        s=30, color="#2a78d6", alpha=0.75, edgecolor="white", linewidth=0.5, zorder=3,
    )

    # y = x reference line -- points above it have more left overhangs, below have more right
    ax.plot(lims, lims, color="black", linewidth=0.8, linestyle="--", zorder=1, label="left = right")

    # label every genome
    for _, row in plot_df.iterrows():
        ax.annotate(
            row["accession"], (row["right"], row["left"]),
            fontsize=6, xytext=(4, 4), textcoords="offset points",
        )

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Right overhangs")
    ax.set_ylabel("Left overhangs")
    ax.set_title(f"{is_element} — left vs. right overhangs per genome "
                 f"({len(plot_df)}/{n_total_accessions} genomes)")
    ax.legend(loc="upper left")

    plt.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{is_element}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"{is_element}: {len(plot_df)} genomes -> {out_path}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("accessions_file", type=Path,
                    help="file of genome accessions, one per line")
    p.add_argument("--cluster-root", type=Path, default=DEFAULT_CLUSTER_ROOT,
                    help=f"root dir containing {{accession}}/{{accession}}.cluster_manifest.tsv "
                         f"(default: {DEFAULT_CLUSTER_ROOT})")
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                    help=f"where to write {{is_element}}.png plots "
                         f"(default: {DEFAULT_OUTPUT_DIR})")
    args = p.parse_args()

    accessions = read_accessions(args.accessions_file)
    counts = collect_counts(accessions, args.cluster_root)

    per_genome = (
        counts
        .groupby(["accession", "is_element", "side"])["n_clusters"]
        .sum()
        .unstack("side", fill_value=0)
        .reset_index()
    )
    for side in ("left", "right"):
        if side not in per_genome.columns:
            per_genome[side] = 0

    print(f"{len(accessions)} genomes, {per_genome['is_element'].nunique()} IS elements")

    for is_element, plot_df in per_genome.groupby("is_element"):
        plot_is_element(is_element, plot_df, len(accessions), args.output_dir)


if __name__ == "__main__":
    main()
