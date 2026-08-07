# purpose: BLAST cd-hit cluster representative sequences (overhang clusters from
# manual_pipeline/overhangs.py + scripts/combine_cdhit_clusters.py) against the
# masked (IS-excised) E. coli reference genomes, and determine whether the
# IS-element insertion each overhang came from actually DISRUPTED a gene.
#
# An overhang read has two ends, and only one is biologically meaningful for
# that question: the end immediately adjacent to the real IS-insertion
# junction. The other end just extends away into flanking host DNA and says
# nothing about where the transposon landed. So classification here is a
# single-point check at the junction-adjacent end (determined by `side`), not
# a whole-aligned-span overlap check -- a hit whose span merely brushes past a
# gene, without the junction-relevant end actually landing inside it, is not a
# disruption.
#
# Which end is junction-adjacent (traced in manual_pipeline/overhangs.py's
# extract_overhangs_from_bam, STAGE 3): for side=="left" the LAST base of the
# written overhang is junction-adjacent (its aligned/IS-anchored portion
# follows the clip in the read); for side=="right" the FIRST base is
# junction-adjacent (the aligned portion precedes the clip). This holds
# regardless of read mapping strand, since pysam already normalizes
# query_sequence/cigartuples to reference-forward order before that script
# ever slices anything, and cd-hit-est never revcomps a cluster's stored
# representative sequence either -- so this orientation survives untouched
# into the BLAST query here.
#
# In BLAST tabular output, qstart<=qend always, and sstart<->qstart /
# send<->qend positionally correspond regardless of subject strand (verified
# empirically with a tiny synthetic BLAST db before relying on it here). So:
# for side=="left" the junction-adjacent genomic position is `send`; for
# side=="right" it's `sstart`.
#
# input: one or more {sample}.cd-hit output directories, each already containing
# {sample}__{is_element}__{side}.reads.tsv files (combine_cdhit_clusters.py output).
# A "cluster" here is the unit (sample, is_element, side, cluster_id); singleton
# clusters never appear in reads.tsv (combine_cdhit_clusters.py drops them), so
# every cluster this script sees has >= 2 member reads.
#
# output: one long-format TSV, one row per (cluster, ref_genome) hit that passed
# the identity/coverage thresholds, plus one row per cluster with no passing hit
# in any of the 6 genomes (ref_genome/hit fields left empty, hit_type="no_hit").
# A cluster can therefore appear on multiple rows if its rep sequence hits more
# than one reference genome -- this is intentional (the flanking sequence isn't
# guaranteed to sit at the same locus, or exist at all, across the 6 assemblies).
# hit_type is "protein_coding" (junction point falls inside a CDS -- a plausible
# gene disruption), "intergenic" (junction point falls outside every CDS, OR the
# alignment didn't reliably reach the junction-relevant end -- conservatively
# treated the same as "not a disruption" rather than a separate uncertain
# state), or "no_hit". A junction point landing in 2+ overlapping/nested CDS
# produces up to 2 rows (gene_rank 1, 2), ranked by margin to the nearest CDS
# edge (deepest-inside-the-gene call first).
#
# usage:
#   python blast_clusters_to_ref.py \
#       --cdhit_dirs /n/scratch/.../cd-hit/SRR8275060 /n/scratch/.../cd-hit/SRR8456430 ... \
#       --blastdb /home/hua575/baymlab/mapping/ref_genomes/blastdb/ecoli_masked_combined \
#       --ref_genomes_dir /home/hua575/baymlab/mapping/ref_genomes/ecoli \
#       --output clusters_blasted.tsv

import argparse
import subprocess
from bisect import bisect_right
from pathlib import Path

import pandas as pd
from Bio import SeqIO

BLASTN_BIN = "/home/hua575/miniconda3/envs/blast/bin/blastn"

BLAST_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore qlen slen sstrand"
)
BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "sstrand",
]

# a hit must clear both of these to count as "the cluster landed here" rather
# than a spurious short/low-identity match
MIN_PIDENT = 95.0
MIN_QCOV = 0.90

# max bp an HSP may fall short of the overhang's true junction-adjacent end
# and still be trusted for gene-disruption/intergenic classification -- an
# HSP can clear MIN_QCOV overall while still stopping short specifically at
# the one end that matters (see compute_junction_position). 0 = the alignment
# must reach the exact junction-adjacent base, no slack, for maximum
# confidence that a "protein_coding" call is a real insertion event.
JUNCTION_TOLERANCE_BP = 0


# ---------------------------------------------------------------------------
# 1. gather cluster representatives from combine_cdhit_clusters.py's reads.tsv
# ---------------------------------------------------------------------------

def gather_cluster_reps(cdhit_dirs: list[Path]) -> pd.DataFrame:
    """
    One row per (sample, is_element, side, cluster_id): its representative
    sequence and how many member reads it has.
    """
    rows = []
    for cdhit_dir in cdhit_dirs:
        for reads_tsv in sorted(cdhit_dir.glob("*.reads.tsv")):
            df = pd.read_csv(reads_tsv, sep="\t", dtype={"pos": str})
            if df.empty:
                continue
            grouped = df.groupby(
                ["sample", "is_element", "side", "cluster_id"], as_index=False
            ).agg(
                n_seqs=("read_id", "count"),
                rep_seq=("representative_sequence", "first"),
            )
            rows.append(grouped)

    if not rows:
        return pd.DataFrame(
            columns=["sample", "is_element", "side", "cluster_id", "n_seqs", "rep_seq"]
        )

    clusters = pd.concat(rows, ignore_index=True)
    clusters["cluster_key"] = (
        clusters["sample"] + "__" + clusters["is_element"] + "__"
        + clusters["side"] + "__cluster" + clusters["cluster_id"].astype(str)
    )
    return clusters


def write_query_fasta(clusters: pd.DataFrame, fasta_path: Path) -> None:
    with open(fasta_path, "w") as fh:
        for _, row in clusters.iterrows():
            fh.write(f">{row['cluster_key']}\n{row['rep_seq']}\n")


# ---------------------------------------------------------------------------
# 2. CDS interval lookup per reference genome, built from the masked gbff
#    (masked, so coordinates match what we're BLASTing against)
# ---------------------------------------------------------------------------

def load_cds_intervals(ref_genomes_dir: Path, accessions: list[str], gbff_suffix: str = "no_IS") -> dict:
    """
    Returns {accession: {contig_id: (starts, ends, strands, genes, max_end)}}.
    All lists sorted by start. max_end[i] is the running max of ends[0..i],
    which lets genes_at_point's backward scan stop as soon as it's provably
    impossible for any earlier-starting CDS to still contain the query point.
    Gene label falls back product -> locus_tag when /gene is absent (bakta
    only sets /gene for named genes; most CDSs only carry /product).

    gbff_suffix: matches the masked-genome variant the --blastdb was built
    from, e.g. "no_IS" (masking.py's original annotation-based excision) or
    "no_IS_v2" (masking.py's ISfinder-BLAST-based excision) -- reads
    {ref_genomes_dir}/{acc}/{acc}_{gbff_suffix}.gbff.
    """
    intervals = {}
    for acc in accessions:
        gbff_path = ref_genomes_dir / acc / f"{acc}_{gbff_suffix}.gbff"
        per_contig = {}
        for record in SeqIO.parse(gbff_path, "genbank"):
            starts, ends, strands, genes = [], [], [], []
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                starts.append(int(feature.location.start))
                ends.append(int(feature.location.end))
                strands.append("+" if feature.location.strand == 1 else "-")
                gene = feature.qualifiers.get("gene", [None])[0]
                if not gene:
                    gene = feature.qualifiers.get("product", [None])[0]
                if not gene:
                    gene = feature.qualifiers.get("locus_tag", ["unnamed_CDS"])[0]
                genes.append(gene)

            order = sorted(range(len(starts)), key=lambda i: starts[i])
            starts = [starts[i] for i in order]
            ends = [ends[i] for i in order]
            strands = [strands[i] for i in order]
            genes = [genes[i] for i in order]

            max_end = []
            running = float("-inf")
            for e in ends:
                running = max(running, e)
                max_end.append(running)

            per_contig[record.id] = (starts, ends, strands, genes, max_end)
        intervals[acc] = per_contig
    return intervals


def genes_at_point(intervals: dict, accession: str, contig_id: str, point0: int,
                    max_hits: int = 2) -> list[str]:
    """
    point0: 0-indexed genomic coordinate of the single junction-adjacent base
    (see compute_junction_position) -- a POINT-containment lookup, not a
    span-overlap one. Callers must not pass a (start, end) range here.

    Returns up to `max_hits` gene names for every CDS whose half-open
    [start, end) interval contains point0, ranked by margin = distance from
    point0 to the nearer edge of that CDS, descending (largest margin first
    = point sits most deeply inside that CDS -- least likely to be a
    boundary-annotation or alignment-slop artifact). Empty list if point0
    falls inside no CDS on this contig.

    Tie-break (deterministic, for equal margins -- e.g. perfectly nested or
    symmetric CDS pairs): longer CDS wins; if still tied, gene name ascending.
    """
    contig_intervals = intervals.get(accession, {}).get(contig_id)
    if not contig_intervals or not contig_intervals[0]:
        return []

    starts, ends, strands, genes, max_end = contig_intervals
    # candidates: any CDS whose start is before or at our point
    idx = bisect_right(starts, point0)
    candidates = []
    for i in range(idx - 1, -1, -1):
        if max_end[i] <= point0:
            # no CDS at or before i can possibly contain point0 anymore
            break
        if starts[i] <= point0 < ends[i]:
            margin = min(point0 - starts[i], ends[i] - 1 - point0)
            candidates.append((genes[i], ends[i] - starts[i], margin))

    # rank: margin desc, then CDS length desc, then gene name asc. Negating
    # the numeric fields (rather than sort(..., reverse=True)) keeps the
    # gene-name tie-break ascending -- reverse=True would flip that too.
    candidates.sort(key=lambda c: (-c[2], -c[1], c[0]))
    return [c[0] for c in candidates[:max_hits]]


def compute_junction_position(h: pd.Series, side: str,
                               tolerance_bp: int = JUNCTION_TOLERANCE_BP) -> tuple[int, bool]:
    """
    h: one row of the `hits` DataFrame (has qstart/qend/qlen/sstart/send from
    BLAST_OUTFMT). side: the cluster's "left" or "right".

    Returns (junction_pos0, reach_ok):
      junction_pos0 -- 0-indexed genomic coordinate, in the matched contig's
        coordinate space, of the single junction-adjacent query base (see
        module docstring for the side->end and qstart/qend<->sstart/send
        correspondence this relies on -- verified empirically, not just
        assumed from BLAST documentation).
      reach_ok -- True iff the alignment actually extends to within
        `tolerance_bp` of the overhang's true junction-adjacent end. Passing
        the overall MIN_QCOV filter does NOT guarantee this: an HSP can lose
        disproportionate coverage right at one end (e.g. a mismatch-dense
        stretch at the true junction gets clipped by blastn) while still
        covering >=90% of the query overall.
    """
    if side == "left":
        junction_pos0 = int(h["send"]) - 1
        reach_ok = (int(h["qlen"]) - int(h["qend"])) <= tolerance_bp
    elif side == "right":
        junction_pos0 = int(h["sstart"]) - 1
        reach_ok = (int(h["qstart"]) - 1) <= tolerance_bp
    else:
        raise ValueError(f"unexpected side value: {side!r}")
    return junction_pos0, reach_ok


# ---------------------------------------------------------------------------
# 3. run blastn, keep best hit per (cluster, genome)
# ---------------------------------------------------------------------------

def run_blast(query_fasta: Path, blastdb: str, threads: int = 8) -> pd.DataFrame:
    cmd = [
        BLASTN_BIN, "-task", "blastn",
        "-query", str(query_fasta),
        "-db", blastdb,
        "-outfmt", BLAST_OUTFMT,
        "-num_threads", str(threads),
        "-evalue", "1e-10",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"blastn failed:\n{result.stderr}")

    from io import StringIO
    if not result.stdout.strip():
        return pd.DataFrame(columns=BLAST_COLUMNS)
    df = pd.read_csv(StringIO(result.stdout), sep="\t", names=BLAST_COLUMNS)
    return df


def best_hit_per_genome(blast_df: pd.DataFrame) -> pd.DataFrame:
    """Split sseqid into (ref_genome, contig), filter by identity/coverage,
    then keep the single best (by bitscore) hit per (qseqid, ref_genome)."""
    if blast_df.empty:
        return blast_df.assign(ref_genome=[], contig=[])

    df = blast_df.copy()
    df[["ref_genome", "contig"]] = df["sseqid"].str.split("__", n=1, expand=True)
    df["qcov"] = df["length"] / df["qlen"]
    df = df[(df["pident"] >= MIN_PIDENT) & (df["qcov"] >= MIN_QCOV)]
    df = df.sort_values("bitscore", ascending=False)
    df = df.drop_duplicates(subset=["qseqid", "ref_genome"], keep="first")
    return df


# ---------------------------------------------------------------------------
# 4. assemble final long-format table
# ---------------------------------------------------------------------------

def assemble_table(clusters: pd.DataFrame, hits: pd.DataFrame, intervals: dict) -> pd.DataFrame:
    hits_by_query = {k: v for k, v in hits.groupby("qseqid")}

    rows = []
    for _, cl in clusters.iterrows():
        key = cl["cluster_key"]
        base = {
            "sample": cl["sample"],
            "is_element": cl["is_element"],
            "side": cl["side"],
            "cluster_id": cl["cluster_id"],
            "rep_seq": cl["rep_seq"],
            "query_length": len(cl["rep_seq"]),
            "n_seqs": cl["n_seqs"],
        }

        cluster_hits = hits_by_query.get(key)
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
                "hit_start": start0 + 1,  # whole-span context only, 1-based inclusive
                "hit_end": end0,
                "hit_strand": h["sstrand"],
                "junction_pos": junction_pos0 + 1,  # the coordinate actually used for classification
                "pident": h["pident"],
                "evalue": h["evalue"],
            }

            genes = (
                genes_at_point(intervals, h["ref_genome"], h["contig"], junction_pos0)
                if reach_ok else []
            )

            if not genes:
                # covers both: junction point outside every CDS, and the
                # alignment not reliably reaching the junction end at all --
                # conservatively lumped together rather than reported as if
                # a confirmed non-disruption
                rows.append({**row_common, "hit_type": "intergenic", "gene": None, "gene_rank": None})
            else:
                for rank, gene in enumerate(genes, start=1):
                    rows.append({**row_common, "hit_type": "protein_coding", "gene": gene, "gene_rank": rank})

    # explicit column order (not left to dict-insertion order, which differs
    # between the no_hit branch and the per-hit branches above) -- hit_type/
    # gene/gene_rank pulled up front right after cluster_id since that's the
    # first thing anyone reads a row for
    column_order = [
        "sample", "is_element", "side", "cluster_id",
        "hit_type", "gene", "gene_rank",
        "rep_seq", "query_length", "n_seqs",
        "ref_genome", "contig", "junction_pos", "hit_start", "hit_end", "hit_strand",
        "pident", "evalue",
    ]
    return pd.DataFrame(rows)[column_order]


# ---------------------------------------------------------------------------
# 5. write one TSV per (sample, is_element) -- one zip per QUERY sample (the
#    SRA accession the clusters came from), containing one "{is_element}.tsv"
#    per IS element found in that sample. Rows that hit a ref genome and rows
#    that hit nothing live in the SAME per-IS table: "ref_genome" says which
#    genome the hit landed in, and "hit_type" is "no_hit" (with ref_genome
#    empty) for clusters that matched none of the 6 -- no separate no-hit
#    bucket, since which genomes DIDN'T match is just as useful to see
#    side-by-side with the ones that did.
#
# Written straight into the zip via writestr so no loose file ever touches
# disk -- this pipeline has already blown /n/scratch's inode quota once from
# unbounded small-file accumulation (see Snakefile's clip_and_cluster rule
# comment), so anything that fans out per-sample/per-IS-element gets zipped
# rather than left as loose files.
# ---------------------------------------------------------------------------

def write_split_tables(final: pd.DataFrame, output_dir: Path) -> None:
    import zipfile
    from io import StringIO

    output_dir.mkdir(parents=True, exist_ok=True)

    n_files = 0
    for sample, sample_df in final.groupby("sample"):
        zip_path = output_dir / f"{sample}.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for is_element, grp in sample_df.groupby("is_element"):
                buf = StringIO()
                grp.to_csv(buf, sep="\t", index=False)
                zf.writestr(f"{is_element}.tsv", buf.getvalue())
                n_files += 1

    print(f"Wrote {n_files} per-IS-element tables into "
          f"{final['sample'].nunique()} per-sample zip archives under {output_dir}/")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cdhit_dirs", nargs="+", required=True,
                    help="one or more {sample}'s cd-hit output dirs (containing *.reads.tsv)")
    p.add_argument("--blastdb", required=True,
                    help="path prefix of the combined masked-genome BLAST db")
    p.add_argument("--ref_genomes_dir", required=True,
                    help="dir containing {accession}/{accession}_{gbff_suffix}.gbff for each ref genome")
    p.add_argument("--gbff_suffix", default="no_IS",
                    help="masked-genome variant to read gbff annotations from -- must match "
                         "what --blastdb was built from, e.g. 'no_IS' (default) or 'no_IS_v2'")
    p.add_argument("--accessions", nargs="+", required=True,
                    help="ref genome accessions making up --blastdb (must match its sseqid prefixes)")
    p.add_argument("--output_dir", required=True,
                    help="output directory; writes {output_dir}/{sample}.zip, each containing "
                         "one {is_element}.tsv (ref_genome/hit_type columns note which genome, "
                         "if any, each cluster hit)")
    p.add_argument("--tmp_dir", default=None,
                    help="dir for the intermediate BLAST query FASTA (deleted by the caller "
                         "afterward, e.g. a job's $TMPDIR). Defaults to --output_dir, which is "
                         "fine for one-off/ad hoc runs, but MUST be set to a per-job-unique dir "
                         "(not shared --output_dir) when multiple invocations of this script run "
                         "concurrently against the same --output_dir -- otherwise their query "
                         "FASTAs collide on the same filename.")
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args()

    cdhit_dirs = [Path(d) for d in args.cdhit_dirs]
    ref_genomes_dir = Path(args.ref_genomes_dir)
    output_dir = Path(args.output_dir)
    tmp_dir = Path(args.tmp_dir) if args.tmp_dir else output_dir

    print("Gathering cluster representatives...")
    clusters = gather_cluster_reps(cdhit_dirs)
    print(f"  {len(clusters)} clusters across {clusters['sample'].nunique() if len(clusters) else 0} samples")

    if clusters.empty:
        print("No clusters found -- nothing to write.")
        return

    tmp_dir.mkdir(parents=True, exist_ok=True)
    query_fasta = tmp_dir / "_clusters.query.fasta"
    write_query_fasta(clusters, query_fasta)

    print("Running blastn against masked ref genomes...")
    blast_df = run_blast(query_fasta, args.blastdb, threads=args.threads)
    print(f"  {len(blast_df)} raw HSPs")

    hits = best_hit_per_genome(blast_df)
    print(f"  {len(hits)} hits pass pident>={MIN_PIDENT} qcov>={MIN_QCOV}")

    print("Loading CDS intervals from masked gbff files...")
    intervals = load_cds_intervals(ref_genomes_dir, args.accessions, gbff_suffix=args.gbff_suffix)

    print("Classifying hits and assembling final table...")
    final = assemble_table(clusters, hits, intervals)
    write_split_tables(final, output_dir)


if __name__ == "__main__":
    main()
