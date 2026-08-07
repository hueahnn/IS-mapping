# purpose: remove IS elements from reference genomes (both from fasta and gbff
# annotation), where IS-element coordinates are found by BLASTing each
# genome's own sequence against the ISfinder reference database rather than
# by trusting Bakta's /product annotation text -- ISfinder is a curated,
# comprehensive IS collection, so this catches divergent/degenerate IS copies
# that never got annotated as a transposase CDS by Bakta in the first place.

import argparse
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

BLASTN_BIN = "/home/hua575/miniconda3/envs/blast/bin/blastn"
MAKEBLASTDB_BIN = "/home/hua575/miniconda3/envs/blast/bin/makeblastdb"

# curated, deduplicated (99% identity) ISfinder collection -- see
# ISfinder_database-master_2026/ for how IS.database.collapsed99.fa itself
# was built (fetch_one_is.py + a 99%-identity collapse pass)
ISFINDER_FASTA = Path(
    "/home/hua575/baymlab/mapping/ISfinder_database-master_2026/IS.database.collapsed99.fa"
)
ISFINDER_BLASTDB = ISFINDER_FASTA.with_suffix("")  # makeblastdb -out prefix (no .fa)

BASE_OUTPUT_DIR = Path("/home/hua575/baymlab/mapping/ref_genomes/ecoli")

# the 6 E. coli reference genomes this pipeline masks (must match Snakefile's
# REF_ACCESSIONS) -- each has an unmasked, Bakta-annotated genomic.gbff at
# BASE_OUTPUT_DIR/{acc}/genomic.gbff
REF_ACCESSIONS = [
    "GCA_000005845.2", "GCA_900096825.1", "GCA_000692435.1",
    "GCA_002473875.1", "GCA_000163235.1", "GCA_002966755.1",
]

BLAST_OUTFMT = "6 qseqid qstart qend pident length evalue"

# a BLAST hit against ISfinder must clear both of these to count as a real IS
# copy worth excising, rather than a short/spurious low-identity match
MIN_PIDENT = 90.0
MIN_LENGTH = 100
EVALUE = "1e-10"


def ensure_isfinder_blastdb(isfinder_fasta: Path, blastdb: Path) -> None:
    """Build the ISfinder nucleotide BLAST db from its FASTA if not already built."""
    if blastdb.with_suffix(".nin").exists():
        return
    print(f"Building ISfinder BLAST db at {blastdb}...")
    cmd = [
        MAKEBLASTDB_BIN, "-in", str(isfinder_fasta),
        "-dbtype", "nucl", "-out", str(blastdb),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"makeblastdb failed:\n{result.stderr}")


def find_is_regions_via_blast(
    genome_fasta_path: Path,
    blastdb: Path,
    min_pident: float = MIN_PIDENT,
    min_length: int = MIN_LENGTH,
    evalue: str = EVALUE,
    threads: int = 8,
) -> dict[str, list[tuple[int, int]]]:
    """
    BLAST a genome's own sequence (query) against the ISfinder db (subject)
    to find where IS elements sit on the genome.
    Returns dict of {contig_id: [(start, end), ...]} -- 0-indexed, half-open
    (Biopython convention), one entry per HSP that clears min_pident/min_length.
    """
    cmd = [
        BLASTN_BIN, "-task", "blastn",
        "-query", str(genome_fasta_path),
        "-db", str(blastdb),
        "-outfmt", BLAST_OUTFMT,
        "-num_threads", str(threads),
        "-evalue", evalue,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"blastn failed:\n{result.stderr}")

    is_coords: dict[str, list[tuple[int, int]]] = {}
    for line in result.stdout.splitlines():
        qseqid, qstart, qend, pident, length, _evalue = line.split("\t")
        pident = float(pident)
        length = int(length)
        if pident < min_pident or length < min_length:
            continue
        # qstart<=qend always in blastn tabular output regardless of subject
        # strand, so this is already query-forward orientation
        start0 = int(qstart) - 1
        end0 = int(qend)
        is_coords.setdefault(qseqid, []).append((start0, end0))

    for contig_id, regions in is_coords.items():
        print(f"  {contig_id}: {len(regions)} ISfinder BLAST hits found")

    return is_coords


def merge_overlapping(coords: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent IS regions before excision."""
    if not coords:
        return []
    sorted_coords = sorted(coords)
    merged = [sorted_coords[0]]
    for start, end in sorted_coords[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def remap_feature(feature: SeqFeature, excised: list[tuple[int, int]]) -> SeqFeature | None:
    """
    Remap a feature's coordinates after IS excision.
    Returns None if the feature falls entirely within an excised region.
    """
    # compute cumulative bases removed before each position
    def shift(pos: int) -> int | None:
        removed = 0
        for start, end in excised:
            if pos <= start:
                break
            elif pos >= end:
                removed += (end - start)
            else:
                return None  # position is inside an excised region
        return pos - removed

    new_start = shift(int(feature.location.start))
    new_end = shift(int(feature.location.end))

    if new_start is None or new_end is None:
        return None  # feature overlaps excised region — drop it

    new_feature = SeqFeature(
        FeatureLocation(new_start, new_end, strand=feature.location.strand),
        type=feature.type,
        qualifiers=feature.qualifiers,
    )
    return new_feature


def excise_is_elements(
    gbff_path: str,
    output_fasta_path: str,
    output_gbff_path: str,
    isfinder_db: Path = ISFINDER_BLASTDB,
    min_pident: float = MIN_PIDENT,
    min_length: int = MIN_LENGTH,
    evalue: str = EVALUE,
    threads: int = 8,
) -> None:
    """
    Remove IS elements (found by BLASTing the genome against ISfinder) from
    both sequence and annotation in a GBFF file.
    Writes cleaned FASTA and cleaned GBFF.
    """
    records = list(SeqIO.parse(gbff_path, "genbank"))

    # write the genome's own sequence out as a BLAST query, using the same
    # record.id as the gbff so ISfinder hits key back onto these records
    # unambiguously
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp_query_fasta = Path(tmp.name)
        SeqIO.write(records, tmp, "fasta")

    try:
        print("BLASTing genome against ISfinder database...")
        is_coords = find_is_regions_via_blast(
            tmp_query_fasta, isfinder_db,
            min_pident=min_pident, min_length=min_length,
            evalue=evalue, threads=threads,
        )
    finally:
        tmp_query_fasta.unlink()

    if not is_coords:
        print("WARNING: No IS elements found via ISfinder BLAST.")
        return

    cleaned_fasta_records = []
    cleaned_gbff_records = []

    for record in records:
        coords = is_coords.get(record.id)

        if coords is None:
            print(f"  No IS elements for {record.id} — keeping as-is")
            cleaned_fasta_records.append(record)
            cleaned_gbff_records.append(record)
            continue

        merged = merge_overlapping(coords)
        total_excised = sum(e - s for s, e in merged)
        print(f"  {record.id}: excising {len(merged)} regions "
              f"({total_excised} bp, {100*total_excised/len(record.seq):.2f}% of genome)")

        # --- excise sequence ---
        segments = []
        prev_end = 0
        for start, end in merged:
            if start > prev_end:
                segments.append(record.seq[prev_end:start])
            prev_end = end
        segments.append(record.seq[prev_end:])
        cleaned_seq = sum(segments, record.seq[:0])

        # --- remap annotations ---
        cleaned_features = []
        skipped = 0
        for feature in record.features:
            remapped = remap_feature(feature, merged)
            if remapped is not None:
                cleaned_features.append(remapped)
            else:
                skipped += 1

        print(f"    {skipped} features dropped (overlapping excised IS regions)")
        print(f"    {len(cleaned_features)} features retained")

        # --- build cleaned record ---
        cleaned_record = SeqRecord(
            cleaned_seq,
            id=record.id,
            name=record.name,
            description=f"{record.description} | IS elements excised ({total_excised} bp removed)",
            dbxrefs=record.dbxrefs,
            annotations=record.annotations,
            features=cleaned_features,
        )
        # update sequence length annotation
        cleaned_record.annotations["sequence_version"] = record.annotations.get("sequence_version", 1)

        cleaned_fasta_records.append(cleaned_record)
        cleaned_gbff_records.append(cleaned_record)

    SeqIO.write(cleaned_fasta_records, output_fasta_path, "fasta")
    print(f"\nWrote FASTA -> {output_fasta_path}")

    SeqIO.write(cleaned_gbff_records, output_gbff_path, "genbank")
    print(f"Wrote GBFF  -> {output_gbff_path}")


def mask_one(
    gbff_path: Path,
    sample: str,
    isfinder_db: Path,
    min_pident: float,
    min_length: int,
    evalue: str,
    threads: int,
) -> None:
    output_dir = BASE_OUTPUT_DIR / sample
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {sample} ===")
    excise_is_elements(
        gbff_path=str(gbff_path),
        output_fasta_path=str(output_dir / f"{sample}_no_IS_v2.fasta"),
        output_gbff_path=str(output_dir / f"{sample}_no_IS_v2.gbff"),
        isfinder_db=isfinder_db,
        min_pident=min_pident,
        min_length=min_length,
        evalue=evalue,
        threads=threads,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove IS elements (found via ISfinder BLAST) from a reference "
                     "genome's FASTA sequence + GBFF annotation."
    )
    parser.add_argument(
        "gbff_path", nargs="?",
        help="Input GBFF file containing the genome annotation. Not used with --all.",
    )
    parser.add_argument(
        "-s", "--sample",
        help="Sample name. Used as the output subdirectory (under "
             f"{BASE_OUTPUT_DIR}) and output filename prefix. "
             "Defaults to the input GBFF's stem.",
    )
    parser.add_argument(
        "--all", action="store_true",
        help=f"Process all {len(REF_ACCESSIONS)} reference genomes in "
             f"{BASE_OUTPUT_DIR} (reads {{acc}}/genomic.gbff for each of "
             "REF_ACCESSIONS) instead of a single --gbff_path.",
    )
    parser.add_argument("--isfinder_fasta", type=Path, default=ISFINDER_FASTA,
                         help="ISfinder FASTA to build/use as the BLAST db subject.")
    parser.add_argument("--isfinder_db", type=Path, default=ISFINDER_BLASTDB,
                         help="ISfinder BLAST db prefix (built from --isfinder_fasta if missing).")
    parser.add_argument("--min_pident", type=float, default=MIN_PIDENT,
                         help="Minimum %% identity for a BLAST hit to count as an IS copy.")
    parser.add_argument("--min_length", type=int, default=MIN_LENGTH,
                         help="Minimum alignment length (bp) for a BLAST hit to count as an IS copy.")
    parser.add_argument("--evalue", default=EVALUE, help="BLAST e-value threshold.")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    ensure_isfinder_blastdb(args.isfinder_fasta, args.isfinder_db)

    if args.all:
        for acc in REF_ACCESSIONS:
            gbff_path = BASE_OUTPUT_DIR / acc / "genomic.gbff"
            mask_one(
                gbff_path, acc, args.isfinder_db,
                args.min_pident, args.min_length, args.evalue, args.threads,
            )
        return

    if not args.gbff_path:
        parser.error("gbff_path is required unless --all is given")

    gbff_path = Path(args.gbff_path)
    sample = args.sample or gbff_path.stem
    mask_one(
        gbff_path, sample, args.isfinder_db,
        args.min_pident, args.min_length, args.evalue, args.threads,
    )


if __name__ == "__main__":
    main()
