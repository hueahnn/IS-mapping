# purpose: for a given bacterial species, reproduce the ATB + SRA coverage
# filtering pipeline from atb/../SRA_table.ipynb end-to-end:
#   1. atb query for the species -> {slug}_results.csv (one row per BioSample)
#   2. extract BioSample IDs -> {slug}_biosampleids.txt
#   3. filter the full NCBI SRA_Accessions.tab down to those BioSamples
#   4. filter to 10x-100x coverage (Bases in [min_bases, max_bases))
#   5. write the passing run accessions to {slug}_sra_accessions.txt
#   6. join back to the ATB table on run_accession -> {slug}_atb_filtered.csv
#      and {slug}_atb_sra_accessions.txt (the file to point a Snakemake
#      ACCESSIONS_PATH/BATCH_FILE at)
#
# usage:
#   python build_atb_species_table.py "Salmonella enterica" --slug salmonella

import argparse
import subprocess
from pathlib import Path

import pandas as pd
import polars as pl

DEFAULT_SRA_TABLE = "/n/scratch/users/h/hua575/SRA_Accessions.tab"
DEFAULT_ATB_DIR = "/home/hua575/baymlab/mapping/atb"


def slugify(species: str) -> str:
    return species.strip().lower().replace(" ", "_")


def query_atb(species: str, results_csv: Path) -> None:
    print(f"[1/6] atb query --species '{species}' -> {results_csv}")
    subprocess.run(
        [
            "atb", "query",
            "--species", species,
            "--hq-only",
            "--has-assembly",
            "--format", "csv",
            "-o", str(results_csv),
        ],
        check=True,
    )


def write_biosample_ids(results_csv: Path, biosample_ids_path: Path) -> list:
    print(f"[2/6] extracting BioSample IDs -> {biosample_ids_path}")
    atb_df = pd.read_csv(results_csv, low_memory=False)
    ids = atb_df["sample_accession"].dropna().astype(str).tolist()
    biosample_ids_path.write_text("\n".join(ids) + "\n")
    return ids


def filter_sra_table_by_biosample(sra_table_path: str, biosample_ids: list,
                                   filtered_parquet: Path) -> pl.DataFrame:
    print(f"[3/6] scanning {sra_table_path} for {len(biosample_ids)} BioSamples "
          f"-> {filtered_parquet}")
    sra_table = pl.scan_csv(
        sra_table_path,
        separator="\t",
        quote_char=None,  # raw NCBI TSV; stray '"' in free-text fields breaks quoted parsing
    ).select(["Accession", "Submission", "Bases", "Sample", "Type", "BioSample"]).filter(
        (pl.col("Bases") != "-") &
        (pl.col("BioSample") != "-") &
        (pl.col("BioSample").is_in(biosample_ids))
    ).collect()
    sra_table.write_parquet(filtered_parquet)
    return sra_table


def filter_by_coverage(sra_table: pl.DataFrame, min_bases: int, max_bases: int,
                        accessions_path: Path) -> set:
    print(f"[4/6] filtering to Bases in [{min_bases}, {max_bases}) "
          f"({sra_table.height} rows before)")
    sra_table = sra_table.with_columns(pl.col("Bases").cast(pl.Int64, strict=False))
    covered = sra_table.filter(
        (pl.col("Bases") >= min_bases) &
        (pl.col("Bases") < max_bases)
    )
    print(f"       {covered.height} rows pass coverage filter")

    print(f"[5/6] writing passing accessions -> {accessions_path}")
    entries = set(covered["Accession"].to_list())
    accessions_path.write_text("\n".join(sorted(entries)) + "\n")
    return entries


def join_to_atb(results_csv: Path, entries: set, atb_filtered_csv: Path,
                 atb_accessions_path: Path) -> None:
    print(f"[6/6] joining back to {results_csv} -> {atb_filtered_csv}")
    atb_df = pd.read_csv(results_csv, low_memory=False)
    matched = atb_df[atb_df["run_accession"].isin(entries)]
    matched.to_csv(atb_filtered_csv, index=False)
    atb_accessions_path.write_text(
        "\n".join(matched["run_accession"].astype(str)) + "\n"
    )
    print(f"       {len(matched)} samples -> {atb_accessions_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("species", help='exact species name, e.g. "Salmonella enterica"')
    p.add_argument("--slug", default=None,
                    help="short name used for output filenames/subdir "
                         "(default: species name lowercased with spaces -> underscores)")
    p.add_argument("--sra-table", default=DEFAULT_SRA_TABLE,
                    help=f"path to NCBI SRA_Accessions.tab (default: {DEFAULT_SRA_TABLE})")
    p.add_argument("--atb-dir", default=DEFAULT_ATB_DIR,
                    help=f"output directory root, gets a {{slug}}/ subdir (default: {DEFAULT_ATB_DIR})")
    p.add_argument("--min-bases", type=int, default=50_000_000,
                    help="minimum Bases, ~10x for a ~5Mb genome (default: 50000000)")
    p.add_argument("--max-bases", type=int, default=500_000_000,
                    help="maximum Bases, ~100x for a ~5Mb genome (default: 500000000)")
    args = p.parse_args()

    slug = args.slug or slugify(args.species)
    out_dir = Path(args.atb_dir) / slug
    out_dir.mkdir(parents=True, exist_ok=True)

    results_csv = out_dir / f"{slug}_results.csv"
    biosample_ids_path = out_dir / f"{slug}_biosampleids.txt"
    filtered_parquet = out_dir / f"{slug}_sra_filtered_results.parquet"
    sra_accessions_path = out_dir / f"{slug}_sra_accessions.txt"
    atb_filtered_csv = out_dir / f"{slug}_atb_filtered.csv"
    atb_accessions_path = out_dir / f"{slug}_atb_sra_accessions.txt"

    query_atb(args.species, results_csv)
    biosample_ids = write_biosample_ids(results_csv, biosample_ids_path)
    sra_table = filter_sra_table_by_biosample(args.sra_table, biosample_ids, filtered_parquet)
    entries = filter_by_coverage(sra_table, args.min_bases, args.max_bases, sra_accessions_path)
    join_to_atb(results_csv, entries, atb_filtered_csv, atb_accessions_path)


if __name__ == "__main__":
    main()
