# purpose: run the overhang pipeline on a per-IS-element level instead of
# per-accession: pool every accession's RAW overhangs (not the already
# per-accession-clustered reads.tsv reps -- see cluster_global_reps.py for
# that variant) for one or more specified IS elements across the entire
# overhangs/ archive, cd-hit-est cluster each (is_element, side) pool, then
# BLAST the resulting cluster representatives against the v2 (ISfinder-BLAST-
# masked) combined reference genome db and classify gene disruption at the
# junction-adjacent end -- same classification logic as
# blast_clusters_to_ref.py's --gbff_suffix no_IS_v2 path, reused directly so
# results stay comparable, just pooled across accessions instead of scoped to
# one. Finally, adds a "pairing" column linking each left-side cluster to the
# right-side cluster (if any) landing back-to-back at the same junction on the
# same reference contig/strand -- same pairing.add_pairing_column.py logic
# used on the per-accession gene_disruption zips, reused directly.
#
# overhangs/{accession}.tar.gz holds the actual per-(is_element, side) fasta
# files; overhangs/{accession}/{accession}.manifest.tsv (extracted alongside
# it) lists which (is_element, side) pairs that accession has and their
# n_seqs, without needing to open the archive -- used here to skip opening an
# accession's tar.gz entirely when it has none of the requested IS elements.
#
# usage:
#   /home/hua575/miniconda3/envs/bakta/bin/python cluster_and_align_by_is.py \
#       --is_elements IS1203 IS621 \
#       --threads 16
#
# scale note: ~148,743 accessions in overhangs/. Run via sbatch, not on a
# login node -- pooling is I/O-bound over that many small manifest files (plus
# opening a tar.gz for every accession that actually has a requested IS
# element), and can take a while. --skip_pooling lets you rerun just the
# cd-hit-est + BLAST steps (e.g. to try a different IS element list against
# already-pooled fastas, or retune thresholds) without repeating it.
#
# memory note: some IS elements collapse into far more distinct clusters than
# others -- confirmed empirically 2026-08-04: pooled across all accessions,
# IS1203 (~30M raw overhangs) collapsed to ~50K clusters, but IS621 (~49M raw
# overhangs) collapsed to ~940K, almost 20x less redundant. cd-hit-est itself
# handles that fine (auto-cycles its lookup table to fit available memory
# rather than OOMing), but gather_cluster_reps's per-cluster streaming here
# keeps memory proportional to cluster count, not raw read count, specifically
# so an IS element with unusually low redundancy doesn't repeat that.

import argparse
import csv
import subprocess
import sys
import tarfile
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from combine_cdhit_clusters import MEMBER_RE, load_sequences, parse_header
from blast_clusters_to_ref import (
    best_hit_per_genome,
    compute_junction_position,
    genes_at_point,
    load_cds_intervals,
    run_blast,
    write_query_fasta,
)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pairing"))
from add_pairing_column import DEFAULT_MAX_GAP, compute_pairing

CD_HIT_EST_BIN = "/home/hua575/miniconda3/envs/cd-hit/bin/cd-hit-est"

# same clustering parameters rule clip_and_cluster uses in
# /home/hua575/baymlab/mapping/Snakefile, so pooled-by-IS clustering stays
# directly comparable to the per-accession clustering it's pooling from
CDHIT_ARGS = ["-c", "0.9", "-n", "8", "-G", "0", "-aS", "0.8", "-d", "0", "-M", "0"]

SIDES = ["left", "right"]

DEFAULT_OVERHANGS_DIR = Path("/n/scratch/users/h/hua575/atb_filtered/overhangs")
DEFAULT_OUTPUT_DIR = Path("/n/scratch/users/h/hua575/atb_filtered/is_level")

REF_GENOMES_DIR = Path("/home/hua575/baymlab/mapping/ref_genomes/ecoli")
BLASTDB_V2 = "/home/hua575/baymlab/mapping/ref_genomes/blastdb/ecoli_masked_combined_v2"
REF_ACCESSIONS = [
    "GCA_000005845.2", "GCA_900096825.1", "GCA_000692435.1",
    "GCA_002473875.1", "GCA_000163235.1", "GCA_002966755.1",
]


# ---------------------------------------------------------------------------
# 1. pool raw overhangs for the requested IS elements across every accession
# ---------------------------------------------------------------------------

def pool_overhangs(overhangs_dir: Path, is_elements: list, output_dir: Path) -> dict:
    """
    Stream every accession's overhangs/{acc}.tar.gz into pooled FASTAs, one
    per (is_element, side), written to
    {output_dir}/{is_element}/{is_element}__{side}.overhangs.fasta.
    Sequence headers are copied through unchanged (still
    "{sample}|{read_id}|{is_element}|pos:...|side:...|strand:...|pident:..."),
    so combine_cdhit_clusters.build_reads_table parses pooled cd-hit output
    exactly like it does per-accession output.
    """
    pooled_paths = {
        (is_el, side): output_dir / is_el / f"{is_el}__{side}.overhangs.fasta"
        for is_el in is_elements for side in SIDES
    }
    for p in pooled_paths.values():
        p.parent.mkdir(parents=True, exist_ok=True)
    handles = {k: open(p, "wb") for k, p in pooled_paths.items()}
    n_seqs = {k: 0 for k in pooled_paths}
    is_elements_set = set(is_elements)

    archives = sorted(overhangs_dir.glob("*.tar.gz"))
    try:
        for i, archive in enumerate(archives, 1):
            sample = archive.name.removesuffix(".tar.gz")
            manifest_path = overhangs_dir / sample / f"{sample}.manifest.tsv"
            if not manifest_path.exists():
                print(f"WARNING: no manifest for {sample}, skipping", file=sys.stderr)
                continue

            wanted = []
            with open(manifest_path) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if row["is_element"] in is_elements_set and int(row["n_seqs"]) > 0:
                        wanted.append((row["is_element"], row["side"]))
            if not wanted:
                continue

            with tarfile.open(archive, "r:gz") as tf:
                index = {Path(m.name).name: m for m in tf.getmembers() if m.isfile()}
                for is_el, side in wanted:
                    member_name = f"{sample}__{is_el}__{side}.overhangs.fasta"
                    member = index.get(member_name)
                    if member is None:
                        print(f"WARNING: {member_name} listed in manifest but missing "
                              f"from {archive.name}", file=sys.stderr)
                        continue
                    with tf.extractfile(member) as mf:
                        data = mf.read()
                    handles[(is_el, side)].write(data)
                    n_seqs[(is_el, side)] += data.count(b">")

            if i % 5000 == 0:
                print(f"[{i}/{len(archives)}] accessions scanned "
                      f"({ {k: v for k, v in n_seqs.items()} })", file=sys.stderr, flush=True)
    finally:
        for fh in handles.values():
            fh.close()

    return pooled_paths, n_seqs


# ---------------------------------------------------------------------------
# 2. cd-hit-est cluster each pooled (is_element, side) fasta
# ---------------------------------------------------------------------------

def run_cdhit(fasta_path: Path, out_prefix: Path, threads: int) -> None:
    clstr_path = out_prefix.with_name(out_prefix.name + ".clstr")
    # incremental-build check (like make): lets a resumed run after a crash
    # further down the pipeline (e.g. the 2026-08-04 OOM below) skip redoing
    # cd-hit-est runs that already completed against the current pooled fasta,
    # rather than unconditionally re-clustering tens of millions of sequences
    if (out_prefix.exists() and clstr_path.exists()
            and out_prefix.stat().st_mtime >= fasta_path.stat().st_mtime):
        print(f"  {out_prefix} already up to date, skipping cd-hit-est", file=sys.stderr)
        return
    cmd = [CD_HIT_EST_BIN, "-i", str(fasta_path), "-o", str(out_prefix),
           *CDHIT_ARGS, "-T", str(threads)]
    print(f"  running: {' '.join(cmd)}", file=sys.stderr, flush=True)
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------------------
# 3. gather non-singleton cluster reps, streamed straight from the .clstr
#    file's member lines rather than through combine_cdhit_clusters.py's
#    build_reads_table (which builds one row per MEMBER READ via its
#    parse_clusters(), keeping every cluster's full member list resident at
#    once -- fine per-accession, but at this script's pooled-across-148k-
#    accessions scale that's tens of millions of rows, which OOM-killed a
#    32GB job building IS621's reads.tsv on 2026-08-04. Nothing downstream
#    of this (BLAST, pairing) ever needs per-read detail, only the
#    per-cluster aggregate, so this keeps only a running counter + a small
#    per-in-progress-cluster sample set in memory, discarded once each
#    cluster is flushed -- final memory is proportional to CLUSTER count
#    (hundreds of thousands here), not raw read count (tens of millions).
# ---------------------------------------------------------------------------

def gather_cluster_reps(cdhit_fasta_path: Path, clstr_path: Path) -> pd.DataFrame:
    representative_sequences = load_sequences(cdhit_fasta_path)
    columns = ["is_element", "side", "cluster_id", "n_seqs", "n_samples", "rep_seq"]
    rows = []

    cluster_id = -1
    member_count = 0
    sample_set = set()
    rep_header = None

    def flush():
        if rep_header is None or member_count < 2:
            return  # empty or singleton cluster -- drop, same as build_reads_table
        rep_fields = parse_header(rep_header)
        rows.append({
            "is_element": rep_fields["is_element"],
            "side": rep_fields["side"],
            "cluster_id": cluster_id,
            "n_seqs": member_count,
            "n_samples": len(sample_set),
            "rep_seq": representative_sequences[rep_header],
        })

    with open(clstr_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                flush()
                cluster_id += 1
                member_count = 0
                sample_set = set()
                rep_header = None
                continue
            m = MEMBER_RE.match(line)
            if not m:
                raise ValueError(f"unrecognized .clstr line: {line!r}")
            header = m.group("header")
            member_count += 1
            sample_set.add(header.split("|", 1)[0])  # "{sample}|..." -- cheaper than a full parse_header per member
            if line.endswith("*"):
                rep_header = header
        flush()

    df = pd.DataFrame(rows, columns=columns)
    df["cluster_key"] = (
        df["is_element"] + "__" + df["side"] + "__cluster" + df["cluster_id"].astype(str)
        if not df.empty else pd.Series(dtype=str)
    )
    return df


# ---------------------------------------------------------------------------
# 4. assemble final long-format table -- same shape/classification logic as
#    blast_clusters_to_ref.py's assemble_table, minus the per-accession
#    "sample" column (replaced by n_samples, carried through from step 3)
# ---------------------------------------------------------------------------

def assemble_pooled_table(clusters: pd.DataFrame, hits: pd.DataFrame, intervals: dict) -> pd.DataFrame:
    hits_by_query = {k: v for k, v in hits.groupby("qseqid")}

    rows = []
    for _, cl in clusters.iterrows():
        base = {
            "is_element": cl["is_element"],
            "side": cl["side"],
            "cluster_id": cl["cluster_id"],
            "rep_seq": cl["rep_seq"],
            "query_length": len(cl["rep_seq"]),
            "n_seqs": cl["n_seqs"],
            "n_samples": cl["n_samples"],
        }

        cluster_hits = hits_by_query.get(cl["cluster_key"])
        if cluster_hits is None or cluster_hits.empty:
            rows.append({
                **base, "ref_genome": None, "hit_type": "no_hit", "gene": None,
                "gene_rank": None, "contig": None, "hit_start": None, "hit_end": None,
                "hit_strand": None, "junction_pos": None, "pident": None, "evalue": None,
            })
            continue

        for _, h in cluster_hits.iterrows():
            start0, end0 = sorted((int(h["sstart"]) - 1, int(h["send"])))
            junction_pos0, reach_ok = compute_junction_position(h, cl["side"])

            row_common = {
                **base,
                "ref_genome": h["ref_genome"],
                "contig": h["contig"],
                "hit_start": start0 + 1,
                "hit_end": end0,
                "hit_strand": h["sstrand"],
                "junction_pos": junction_pos0 + 1,
                "pident": h["pident"],
                "evalue": h["evalue"],
            }

            genes = (
                genes_at_point(intervals, h["ref_genome"], h["contig"], junction_pos0)
                if reach_ok else []
            )

            if not genes:
                rows.append({**row_common, "hit_type": "intergenic", "gene": None, "gene_rank": None})
            else:
                for rank, gene in enumerate(genes, start=1):
                    rows.append({**row_common, "hit_type": "protein_coding", "gene": gene, "gene_rank": rank})

    column_order = [
        "is_element", "side", "cluster_id",
        "hit_type", "gene", "gene_rank",
        "rep_seq", "query_length", "n_seqs", "n_samples",
        "ref_genome", "contig", "junction_pos", "hit_start", "hit_end", "hit_strand",
        "pident", "evalue",
    ]
    return pd.DataFrame(rows)[column_order]


# ---------------------------------------------------------------------------
# 5. pair up left/right clusters landing at the same junction -- same
#    read-modify-rewrite as add_pairing_column.py's process_zip, just on a
#    plain TSV on disk instead of a TSV member inside a per-accession zip
# ---------------------------------------------------------------------------

def add_pairing(tsv_path: Path, max_gap: float) -> None:
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    if "pairing" in fieldnames:
        fieldnames.remove("pairing")
    fieldnames.insert(fieldnames.index("cluster_id") + 1, "pairing")

    compute_pairing(rows, max_gap)

    with open(tsv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--is_elements", nargs="+", default=["IS1203", "IS621"],
                    help="IS element name(s) to pool overhangs for, e.g. IS1203 IS621 "
                         "(must match the is_element names used in overhangs.py's manifests)")
    p.add_argument("--overhangs_dir", type=Path, default=DEFAULT_OVERHANGS_DIR,
                    help="per-accession overhangs dir (contains {acc}.tar.gz + "
                         "{acc}/{acc}.manifest.tsv)")
    p.add_argument("--output_dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                    help="where to write {is_element}/ subdirs with pooled fastas, "
                         "cd-hit output, clusters.tsv, and the final v2-alignment TSV")
    p.add_argument("--blastdb", default=BLASTDB_V2,
                    help="path prefix of the v2 masked combined ref genome BLAST db")
    p.add_argument("--ref_genomes_dir", type=Path, default=REF_GENOMES_DIR,
                    help="dir containing {accession}/{accession}_{gbff_suffix}.gbff")
    p.add_argument("--gbff_suffix", default="no_IS_v2",
                    help="masked-genome variant to read gbff annotations from -- must match "
                         "what --blastdb was built from (default: no_IS_v2)")
    p.add_argument("--accessions", nargs="+", default=REF_ACCESSIONS,
                    help="ref genome accessions making up --blastdb")
    p.add_argument("--max_gap", type=float, default=DEFAULT_MAX_GAP,
                    help="max allowed |right_junction - left_junction| in bp for pairing "
                         f"left/right clusters at the same insertion site (default {DEFAULT_MAX_GAP}, "
                         "same as pairing/add_pairing_column.py)")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--skip_pooling", action="store_true",
                    help="reuse existing pooled {is_element}/{is_element}__{side}.overhangs.fasta "
                         "in --output_dir and only (re-)run cd-hit-est + BLAST")
    args = p.parse_args()

    if args.skip_pooling:
        pooled_paths = {
            (is_el, side): args.output_dir / is_el / f"{is_el}__{side}.overhangs.fasta"
            for is_el in args.is_elements for side in SIDES
        }
        missing = [str(path) for path in pooled_paths.values() if not path.exists()]
        if missing:
            sys.exit(f"--skip_pooling but missing pooled fasta(s): {missing}")
    else:
        t0 = time.time()
        pooled_paths, n_seqs = pool_overhangs(args.overhangs_dir, args.is_elements, args.output_dir)
        print(f"pooling took {time.time() - t0:.0f}s: {n_seqs}", file=sys.stderr)

    print("Loading CDS intervals from v2-masked gbff files...")
    intervals = load_cds_intervals(args.ref_genomes_dir, args.accessions, gbff_suffix=args.gbff_suffix)

    for is_element in args.is_elements:
        print(f"=== {is_element} ===")
        is_dir = args.output_dir / is_element

        for side in SIDES:
            out_prefix = is_dir / f"{is_element}__{side}.cdhit"
            print(f"  cd-hit-est ({side})...")
            run_cdhit(pooled_paths[(is_element, side)], out_prefix, args.threads)

        clusters = pd.concat(
            [gather_cluster_reps(is_dir / f"{is_element}__{side}.cdhit",
                                  is_dir / f"{is_element}__{side}.cdhit.clstr")
             for side in SIDES],
            ignore_index=True,
        )
        clusters_tsv_path = is_dir / f"{is_element}.clusters.tsv"
        clusters.drop(columns=["cluster_key"]).to_csv(clusters_tsv_path, sep="\t", index=False)
        print(f"  {len(clusters)} non-singleton clusters -> {clusters_tsv_path}")

        if clusters.empty:
            print(f"  no non-singleton clusters for {is_element}, skipping BLAST")
            continue

        query_fasta = is_dir / f"{is_element}.clusters.query.fasta"
        write_query_fasta(clusters, query_fasta)

        print(f"  blastn against {args.blastdb} ...")
        blast_df = run_blast(query_fasta, args.blastdb, threads=args.threads)
        hits = best_hit_per_genome(blast_df)
        print(f"  {len(hits)}/{len(blast_df)} HSPs pass identity/coverage thresholds")

        final = assemble_pooled_table(clusters, hits, intervals)
        out_tsv = is_dir / f"{is_element}.v2_alignment.tsv"
        final.to_csv(out_tsv, sep="\t", index=False)

        print(f"  pairing left/right clusters (max_gap={args.max_gap})...")
        add_pairing(out_tsv, args.max_gap)
        print(f"  -> {out_tsv} ({len(final)} rows, {len(clusters)} clusters)")


if __name__ == "__main__":
    main()
