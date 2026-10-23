# purpose: find IS elements inserted INTO other IS elements / transposons
# (nested insertions), from the overhang clusters blast_clusters_to_ref.py
# could not place in any reference genome.
#
# The reference genomes are IS-masked, so an overhang whose flanking sequence
# is itself IS/Tn DNA comes out of blast_clusters_to_ref as hit_type=="no_hit".
# This script takes those no_hit clusters, BLASTs their representative
# sequences against the same IS/Tn databases the reads were aligned to, and
# asks whether the flank of host element X lands inside another element Y.
#
# Orientation (same convention as blast_clusters_to_ref.py): rep_seq is in X's
# forward orientation; for side=="left" the LAST base is junction-adjacent,
# for side=="right" the FIRST base is. BLAST keeps qstart<->sstart /
# qend<->send positionally paired on either strand, so the insertion point
# inside Y is `send` for left and `sstart` for right (1-based Y coordinates).
#
# Each distinct rep_seq is BLASTed once (identical flanks recur across
# clusters and, in batch mode, across samples); hits are then joined back onto
# every cluster carrying that sequence. The class depends on the cluster's own
# side and is_element, not just its sequence.
#
# Per-cluster class (best-bitscore hit among hits that reach the junction end):
#   internal     junction lands inside Y: nested-insertion candidate
#   terminal     junction within --terminal_tol bp of either end of Y -- X
#                abuts Y's end (adjacent elements) rather than interrupting Y
#   self         best hit (or a tie for it) is X itself -- tandem copies, IS
#                circles/dimers, or X sequence let in by boundary_tolerance
#   distal_only  IS/Tn hits exist but none reaches within --junction_tol bp
#                of the junction end (host DNA between X and Y: X near Y, not
#                in it)
#   no_is_hit    no hit with pident >= --min_pident and length >= --min_aln
#   too_short    rep_seq shorter than --min_len, not BLASTed
#
# paired_in_target marks clusters whose opposite-side partner (same sample and
# host, other side) lands in the same Y on the same strand with junction
# positions within --pair_tol bp -- both flanks of one X copy pointing into
# one site of Y. This is the strongest evidence for a real nested insertion.
#
# input: one or more gene_disruption_paired/{sample}.zip files (add_pairing_column
#   output; --zips) or a folder of them (--zip_dir), the IS/Tn BLAST db and
#   references.tsv (reference, database, length) from the Snakefile's
#   reference_table rule.
# output: one gzipped TSV, one row per no_hit cluster (header only if none).
#   Summarize one or many of these with scripts/summarize_nested_is.py.
#
# usage:
#   per sample (the Snakefile's nested_is rule):
#     python nested_is.py --zips gene_disruption_paired/SRR123.zip \
#         --blastdb db_index/nested_is/IS_Tn --references db_index/references.tsv \
#         --output nested_is/SRR123.tsv.gz
#   every sample of an existing run at once (one BLAST, cross-sample dedup):
#     python nested_is.py --zip_dir /n/scratch/.../gene_disruption_paired ... \
#         --output nested_is.clusters.tsv.gz --threads 4

import argparse
import csv
import io
import os
import subprocess
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd

BLASTN_BIN = "/home/hua575/miniconda3/envs/blast/bin/blastn"

# columns carried over from the gene-disruption tables
KEEP = ["sample", "is_element", "side", "cluster_id", "pairing", "pairing_type",
        "rep_seq", "query_length", "n_seqs"]
BLAST_OUTFMT = ("qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
                "evalue bitscore qlen slen sstrand")
BLAST_COLS = ["seq_id", "target", "pident", "aln_len", "mismatch", "gapopen", "qstart",
              "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "sstrand"]
OUT_COLS = ["sample", "is_element", "side", "cluster_id", "n_seqs", "query_length",
            "pairing", "pairing_type", "class", "target", "target_db", "target_len",
            "junction_pos", "dist_to_target_end", "sstrand", "paired_in_target",
            "pident", "aln_len", "qcov", "bitscore", "evalue", "n_tied",
            "tied_targets", "rep_seq"]


# ---------------------------------------------------------------------------
# 1. pull the no_hit clusters out of the gene-disruption zips
# ---------------------------------------------------------------------------

def read_no_hit(zip_paths):
    rows, n_bad = [], 0
    for i, path in enumerate(zip_paths):
        try:
            with zipfile.ZipFile(path) as zf:
                for name in zf.namelist():
                    with zf.open(name) as fh:
                        for row in csv.DictReader(io.TextIOWrapper(fh), delimiter="\t"):
                            if row["hit_type"] == "no_hit":
                                rows.append([row[k] for k in KEEP])
        except zipfile.BadZipFile:
            n_bad += 1
            print(f"bad zip, skipped: {path}", file=sys.stderr)
        if (i + 1) % 10000 == 0:
            print(f"  {i + 1}/{len(zip_paths)} zips read, {len(rows)} no_hit clusters",
                  file=sys.stderr, flush=True)
    df = pd.DataFrame(rows, columns=KEEP)
    df["query_length"] = pd.to_numeric(df.query_length)
    print(f"{len(zip_paths)} zips ({n_bad} bad), {len(df)} no_hit clusters", file=sys.stderr)
    return df


# ---------------------------------------------------------------------------
# 2. BLAST each distinct sequence once
# ---------------------------------------------------------------------------

def run_blast(seqs, blastdb, blastn, threads, tmp_dir):
    """seqs: unique sequences. Returns (seq -> seq_id, hits DataFrame)."""
    seq_ids = {s: f"q{i:08d}" for i, s in enumerate(seqs)}
    if not seq_ids:
        return seq_ids, pd.DataFrame(columns=BLAST_COLS)
    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp:
        query = os.path.join(tmp, "queries.fa")
        with open(query, "w") as fh:
            for s, sid in seq_ids.items():
                fh.write(f">{sid}\n{s}\n")
        out = os.path.join(tmp, "hits.tsv")
        # megablast with word size 16 matched full blastn's sensitivity on
        # these flanks (766 vs 768 queries with a passing hit) at ~1/5 the time
        cmd = [blastn, "-task", "megablast", "-word_size", "16",
               "-query", query, "-db", blastdb, "-outfmt", f"6 {BLAST_OUTFMT}",
               "-evalue", "1e-5", "-max_target_seqs", "50", "-max_hsps", "1",
               "-num_threads", str(threads), "-out", out]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"blastn failed:\n{result.stderr}")
        hits = pd.read_csv(out, sep="\t", names=BLAST_COLS)
    print(f"{len(seq_ids)} unique query sequences, {len(hits)} raw hits", file=sys.stderr)
    return seq_ids, hits


# ---------------------------------------------------------------------------
# 3. classify every cluster from its sequence's hits
# ---------------------------------------------------------------------------

def classify(clusters, seq_ids, hits, references, a):
    clusters = clusters.copy()
    clusters["seq_id"] = clusters.rep_seq.map(seq_ids)
    hits = hits[(hits.pident >= a.min_pident) & (hits.aln_len >= a.min_aln)]
    any_hit = set(hits.seq_id)

    # which end is the junction depends on side, so expand hits per (seq, side)
    seq_side = clusters.loc[clusters.seq_id.notna(), ["seq_id", "side"]].drop_duplicates()
    h = hits.merge(seq_side, on="seq_id")
    left = h.side == "left"
    h = h[np.where(left, h.qlen - h.qend, h.qstart - 1) <= a.junction_tol].copy()
    h["junction_pos"] = np.where(h.side == "left", h.send, h.sstart)

    # best reaching hit per (seq, side); equal-bitscore ties are kept to test
    # for self and to report
    h["best"] = h.groupby(["seq_id", "side"]).bitscore.transform("max")
    tied = h[h.bitscore == h.best]
    tie_lists = (tied.groupby(["seq_id", "side"]).target
                 .agg(lambda s: ",".join(sorted(set(s)))).rename("tied_targets"))
    best = (tied.sort_values(["seq_id", "side", "pident", "target"],
                             ascending=[True, True, False, True])
            .drop_duplicates(["seq_id", "side"])
            .merge(tie_lists, left_on=["seq_id", "side"], right_index=True))
    best = best[["seq_id", "side", "target", "tied_targets", "pident", "aln_len", "slen",
                 "sstrand", "bitscore", "evalue", "junction_pos"]]

    df = clusters.merge(best, on=["seq_id", "side"], how="left")
    df["n_tied"] = df.tied_targets.str.count(",") + 1
    df["qcov"] = (df.aln_len / df.query_length).round(3)
    df["dist_to_target_end"] = np.minimum(df.junction_pos - 1, df.slen - df.junction_pos)
    is_self = pd.Series([isinstance(t, str) and e in t.split(",")
                         for e, t in zip(df.is_element, df.tied_targets)], index=df.index)
    df["class"] = np.select(
        [df.seq_id.isna(),
         df.target.isna() & ~df.seq_id.isin(any_hit),
         df.target.isna(),
         is_self,
         df.dist_to_target_end <= a.terminal_tol],
        ["too_short", "no_is_hit", "distal_only", "self", "terminal"],
        default="internal")
    df = df.merge(references.rename(columns={"reference": "target", "database": "target_db"})
                  [["target", "target_db"]], on="target", how="left")

    # left/right clusters of one host copy landing at one site of Y
    key = ["sample", "is_element", "target", "sstrand"]
    cand = df[df["class"].isin(["internal", "terminal"])]
    L = cand[cand.side == "left"][key + ["junction_pos", "cluster_id"]]
    R = cand[cand.side == "right"][key + ["junction_pos", "cluster_id"]]
    pairs = L.merge(R, on=key, suffixes=("_l", "_r"))
    pairs = pairs[(pairs.junction_pos_l - pairs.junction_pos_r).abs() <= a.pair_tol]
    paired = {(s, e, "left", c) for s, e, c in
              zip(pairs["sample"], pairs.is_element, pairs.cluster_id_l)} | \
             {(s, e, "right", c) for s, e, c in
              zip(pairs["sample"], pairs.is_element, pairs.cluster_id_r)}
    df["paired_in_target"] = [k in paired for k in
                              zip(df["sample"], df.is_element, df.side, df.cluster_id)]

    return df.rename(columns={"slen": "target_len"})[OUT_COLS]


def main():
    p = argparse.ArgumentParser(description="Classify no_hit overhang clusters against "
                                            "the IS/Tn databases (nested insertions).")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--zips", nargs="+", help="gene_disruption_paired/{sample}.zip file(s)")
    src.add_argument("--zip_dir", help="folder of {sample}.zip files (all are read)")
    p.add_argument("--blastdb", required=True, help="IS/Tn BLAST db prefix")
    p.add_argument("--references", required=True,
                   help="references.tsv (reference, database, length)")
    p.add_argument("--output", required=True, help="per-cluster table (.tsv.gz)")
    p.add_argument("--blastn", default=BLASTN_BIN)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--tmp_dir", default=None)
    p.add_argument("--min_len", type=int, default=30,
                   help="shorter rep_seqs are not BLASTed (class too_short)")
    p.add_argument("--min_pident", type=float, default=90.0)
    p.add_argument("--min_aln", type=int, default=30)
    p.add_argument("--junction_tol", type=int, default=5,
                   help="max bp between the alignment end and the overhang's junction end")
    p.add_argument("--terminal_tol", type=int, default=20,
                   help="junctions this close to either end of Y are class terminal")
    p.add_argument("--pair_tol", type=int, default=20,
                   help="max bp between left and right junction positions in Y")
    a = p.parse_args()

    zips = a.zips or sorted(os.path.join(a.zip_dir, f) for f in os.listdir(a.zip_dir)
                            if f.endswith(".zip"))
    references = pd.read_csv(a.references, sep="\t")

    clusters = read_no_hit(zips)
    if clusters.empty:
        # every cluster was placed in a reference genome (or the sample had
        # none): header-only table, so the sample reads as "nothing nested"
        pd.DataFrame(columns=OUT_COLS).to_csv(a.output, sep="\t", index=False)
        return
    seqs =clusters.rep_seq[clusters.rep_seq.str.len() >= a.min_len].unique()
    seq_ids, hits = run_blast(seqs, a.blastdb, a.blastn, a.threads, a.tmp_dir)
    df = classify(clusters, seq_ids, hits, references, a)
    df.to_csv(a.output, sep="\t", index=False)

    counts = df["class"].value_counts()
    print("clusters per class: " + ", ".join(f"{k}={v}" for k, v in counts.items()),
          file=sys.stderr)


if __name__ == "__main__":
    main()
