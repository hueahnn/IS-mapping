# purpose: fold the ambiguous-MAPQ overhangs that clustering skipped back into the
# clusters that cluster_overhangs_edlib.py already built, so each cluster's read
# count ("cluster coverage") reflects all of its supporting evidence rather than
# only the confidently-mapped subset.
#
# WHY THESE READS EXIST: IS elements carry terminal inverted repeats, so a read
# from one IR aligns about equally well to the other and bwa hands it a low MAPQ.
# overhangs.py stage 1 stamps those reads ZA=1 and stage 3 routes their clipped
# sequences into a per-IS "__ambiguous" FASTA instead of the "__left"/"__right"
# one they would otherwise land in. They are real junction reads sidelined by a
# mapping artifact, not noise -- but their side call is exactly the thing that is
# unreliable, which is why they are held back until representatives exist to test
# them against.
#
# ANCHORING / ORIENTATION: cluster_overhangs_edlib.py compares sequences in a
# junction-first orientation, reversing side=="left" overhangs to get there (see
# its module docstring). A left overhang and a right overhang are different
# genomic sequences -- the two flanks of one insertion -- so they never belong in
# the same cluster. But if bwa put a right-IR read on the left IR, it mapped in
# the opposite orientation, and that read's true right-side overhang is
# revcomp(seq). So each ambiguous overhang is tested under TWO hypotheses:
#
#   A. the side call was right   -> compare cmp_orient(seq, called_side)
#                                   against that side's representatives
#   B. bwa picked the wrong IR   -> compare cmp_orient(revcomp(seq), other_side)
#                                   against the OTHER side's representatives
#
# Each representative is therefore compared under exactly one hypothesis. Both
# comparison sequences are built by running cmp_orient() -- the same rule
# cluster_overhangs_edlib.py applies -- over the flipped overhang, rather than by
# hand-deriving the composed transform. (reverse(revcomp(x)) collapses to
# complement(x), which is correct but reads like a bug.)
#
# SHW ARGUMENT ORDER: greedy_cluster() can pass query=candidate, target=seed
# unconditionally because its length-descending sort guarantees the candidate is
# never longer. Recruitment has no such guarantee -- a recruited overhang can be
# longer than the representative it matches. edlib's SHW ("prefix") mode puts the
# free gap at the end of TARGET only, so the shorter sequence must always be the
# query, or the longer one's tail gets charged as edits. anchored_distance()
# orders them by length for this reason.
#
# ASSIGNMENT: the representative with the lowest length-normalized edit distance
# wins, provided it is within --recruit_edit_frac. Ties across different clusters
# are left UNASSIGNED: these reads were ambiguous to begin with, and forcing a
# call would re-import the artifact this step exists to correct. Every unassigned
# read is counted by reason in the coverage table, so the drop is visible rather
# than silent.
#
# Representatives are never recomputed -- a recruited read joins an existing
# cluster and inherits its (is_element, side, cluster_id) key, so the sequences
# blast_clusters_to_ref.py queries with are unchanged.
#
# usage:
#   python recruit_ambiguous_overhangs.py \
#       --reads_tsv clusters/{sample}/{sample}.reads.tsv \
#       --overhang_dir <dir holding extracted *__ambiguous.overhangs.fasta> \
#       --sample_id SRR12705059 --output_dir out/ \
#       --recruit_edit_frac 0.10

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import edlib
import pandas as pd

# reuse the header/FASTA parsing and the proportional edit budget from the
# clustering step rather than restating them -- if that script's conventions
# change, recruitment must follow automatically or the two will silently disagree
# about what counts as the same cluster.
from cluster_overhangs_edlib import edit_threshold, parse_header, read_fasta

# the reads.tsv schema cluster_overhangs_edlib.py writes, in order
CORE_COLUMNS = [
    "sample", "is_element", "side", "read_id", "pos", "strand", "pident",
    "cluster_id", "representative_read_id", "representative_sequence",
]

# provenance columns this script appends. gather_cluster_reps() in
# blast_clusters_to_ref.py selects only the columns it names and
# add_pairing_column.py copies fieldnames through, so extra columns are inert
# downstream.
RECRUIT_COLUMNS = ["recruited", "called_side", "recruit_norm_dist", "recruit_mapq"]

COVERAGE_COLUMNS = [
    "sample", "row_type", "is_element", "side", "cluster_id", "reason",
    "n_core", "n_recruited", "n_total",
]

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

# float slop when deciding whether two normalized distances are the same -- these
# are small-integer ratios, so anything this close is an exact tie in practice
TIE_EPSILON = 1e-9


def revcomp(seq: str) -> str:
    """Reverse complement. overhangs.py's stage -1 filter guarantees pure ACGT,
    and str.translate passes anything else through unchanged."""
    return seq.translate(_COMPLEMENT)[::-1]


def cmp_orient(seq: str, side: str) -> str:
    """Junction-first comparison orientation -- mirrors cluster_overhangs_edlib.py's
    `cmp_seq = seq[::-1] if side == "left" else seq`. Kept as a named function so
    both hypotheses below go through the identical rule."""
    return seq[::-1] if side == "left" else seq


def other_side(side: str) -> str:
    return "right" if side == "left" else "left"


def anchored_distance(a: str, b: str, max_edit_frac: float) -> float | None:
    """Length-normalized anchored edit distance between two comparison-orientation
    sequences, or None if they are further apart than max_edit_frac allows.

    The shorter sequence is always the query: SHW mode leaves a free gap at the end
    of the target only, so the longer sequence's tail beyond the shorter one is not
    penalized -- the same "anchored, not length-penalized" comparison
    greedy_cluster() relies on. edlib's k bound doubles as the max_edit_frac gate:
    over budget comes back as editDistance == -1.

    Normalizing by the query length lets representatives of different lengths
    compete fairly for the same read.
    """
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    if not q:
        return None
    k = edit_threshold(len(q), max_edit_frac)
    result = edlib.align(q, t, mode="SHW", task="distance", k=k)
    if result["editDistance"] == -1:
        return None
    return result["editDistance"] / len(q)


def load_representatives(reads_df: pd.DataFrame) -> dict[str, list[dict]]:
    """One representative per (is_element, side, cluster_id), grouped by IS element
    -- the same cluster identity gather_cluster_reps() in blast_clusters_to_ref.py
    uses. Sides other than left/right are skipped defensively; after the
    cluster_overhangs_edlib.py fix they cannot appear here."""
    reps_by_is: dict[str, list[dict]] = defaultdict(list)
    if reads_df.empty:
        return reps_by_is

    unique = reads_df.drop_duplicates(subset=["is_element", "side", "cluster_id"])
    for row in unique.itertuples(index=False):
        if row.side not in ("left", "right"):
            continue
        reps_by_is[row.is_element].append({
            "is_element": row.is_element,
            "side": row.side,
            "cluster_id": row.cluster_id,
            "representative_read_id": row.representative_read_id,
            "representative_sequence": row.representative_sequence,
            "cmp_seq": cmp_orient(row.representative_sequence, row.side),
        })
    return reps_by_is


def read_ambiguous_overhangs(overhang_dir: Path) -> list[dict]:
    """Parse every *__ambiguous.overhangs.fasta in overhang_dir.

    Globs the directory rather than reading fasta_path out of the overhang
    manifest: clip_and_cluster copies that manifest verbatim from the job's
    $TMPDIR (only the CLUSTER manifest gets sed-rewritten), so its paths point at
    a directory that no longer exists by the time this runs.
    """
    records = []
    for fasta_path in sorted(overhang_dir.glob("*__ambiguous.overhangs.fasta")):
        for rec in read_fasta(fasta_path):
            # parse_header captures mapq as a named group, None on pre-ZA headers
            rec["mapq"] = rec.get("mapq") or ""
            records.append(rec)
    return records


def assign_one(record: dict, reps: list[dict], recruit_edit_frac: float):
    """Best-matching representative for one ambiguous overhang, or (None, reason).

    Returns (rep, normalized_distance) on a clean win, or (None, "tie") /
    (None, "no_match").
    """
    seq, called = record["seq"], record["side"]
    if called not in ("left", "right"):
        return None, "bad_side"

    # hypothesis A keeps the called side; hypothesis B assumes bwa picked the
    # wrong IR, which flips both the side and the strand of the overhang
    cmp_by_side = {
        called: cmp_orient(seq, called),
        other_side(called): cmp_orient(revcomp(seq), other_side(called)),
    }

    best_rep, best_dist, tied = None, None, False
    for rep in reps:
        cand = cmp_by_side.get(rep["side"])
        if cand is None:
            continue
        dist = anchored_distance(cand, rep["cmp_seq"], recruit_edit_frac)
        if dist is None:
            continue
        if best_dist is None or dist < best_dist - TIE_EPSILON:
            best_rep, best_dist, tied = rep, dist, False
        elif abs(dist - best_dist) <= TIE_EPSILON:
            tied = True

    if best_rep is None:
        return None, "no_match"
    if tied:
        return None, "tie"
    return best_rep, best_dist


def recruit(reads_df: pd.DataFrame, ambiguous: list[dict], sample_id: str,
            recruit_edit_frac: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (augmented reads table, cluster coverage table)."""
    reps_by_is = load_representatives(reads_df)

    recruited_rows = []
    # (is_element, side, cluster_id) -> count, and (is_element, called_side, reason) -> count
    recruited_counts: Counter = Counter()
    unassigned_counts: Counter = Counter()

    for rec in ambiguous:
        reps = reps_by_is.get(rec["is_element"], [])
        rep, outcome = assign_one(rec, reps, recruit_edit_frac)

        if rep is None:
            unassigned_counts[(rec["is_element"], rec["side"], outcome)] += 1
            continue

        # the winning cluster's full (is_element, side, cluster_id) key is what
        # this read joins -- NOT the side parsed from its own header. Writing the
        # called side next to the winner's cluster_id would invent a cluster that
        # does not exist (e.g. "IS2/left/cluster 3" when cluster 3 only exists
        # under "right"), which gather_cluster_reps() would then treat as a real
        # cluster, BLAST a second time, and classify with the wrong junction end.
        # pos/strand below stay as recorded: on a side-flipped recruit they
        # describe the wrong-end alignment, so they are provenance only. Nothing
        # downstream reads them out of reads.tsv.
        recruited_rows.append({
            "sample": sample_id,
            "is_element": rep["is_element"],
            "side": rep["side"],
            "read_id": rec["read_id"],
            "pos": rec["pos"],
            "strand": rec["strand"],
            "pident": rec["pident"],
            "cluster_id": rep["cluster_id"],
            "representative_read_id": rep["representative_read_id"],
            "representative_sequence": rep["representative_sequence"],
            "recruited": 1,
            "called_side": rec["side"],
            "recruit_norm_dist": round(outcome, 6),
            "recruit_mapq": rec["mapq"],
        })
        recruited_counts[(rep["is_element"], rep["side"], rep["cluster_id"])] += 1

    core = reads_df.copy()
    if core.empty:
        core = pd.DataFrame(columns=CORE_COLUMNS)
    core["recruited"] = 0
    for col in ("called_side", "recruit_norm_dist", "recruit_mapq"):
        core[col] = ""

    # concatenating an empty frame is deprecated in pandas (it changes the result
    # dtypes), and "nothing recruited" is the normal case for any sample whose
    # overhangs were extracted before overhangs.py grew the ZA tag
    if recruited_rows:
        recruited_df = pd.DataFrame(recruited_rows, columns=CORE_COLUMNS + RECRUIT_COLUMNS)
        augmented = pd.concat([core, recruited_df], ignore_index=True)
    else:
        augmented = core
    augmented = augmented[CORE_COLUMNS + RECRUIT_COLUMNS]

    coverage = build_coverage(reads_df, recruited_counts, unassigned_counts, sample_id)
    return augmented, coverage


def build_coverage(reads_df: pd.DataFrame, recruited_counts: Counter,
                   unassigned_counts: Counter, sample_id: str) -> pd.DataFrame:
    """Per-cluster coverage (n_core + n_recruited = n_total), followed by rows
    accounting for every ambiguous read that was NOT recruited."""
    rows = []

    if not reads_df.empty:
        core_counts = (
            reads_df.groupby(["is_element", "side", "cluster_id"], as_index=False)
            .agg(n_core=("read_id", "count"))
        )
        for row in core_counts.itertuples(index=False):
            key = (row.is_element, row.side, row.cluster_id)
            n_recruited = recruited_counts.get(key, 0)
            rows.append({
                "sample": sample_id, "row_type": "cluster",
                "is_element": row.is_element, "side": row.side,
                "cluster_id": row.cluster_id, "reason": "",
                "n_core": row.n_core, "n_recruited": n_recruited,
                "n_total": row.n_core + n_recruited,
            })

    for (is_element, called_side, reason), n in sorted(unassigned_counts.items()):
        rows.append({
            "sample": sample_id, "row_type": "unassigned",
            "is_element": is_element, "side": called_side,
            "cluster_id": "", "reason": reason,
            "n_core": "", "n_recruited": n, "n_total": "",
        })

    return pd.DataFrame(rows, columns=COVERAGE_COLUMNS)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reads_tsv", required=True,
                    help="cluster_overhangs_edlib.py's {sample}.reads.tsv")
    p.add_argument("--overhang_dir", required=True,
                    help="directory holding the extracted "
                         "*__ambiguous.overhangs.fasta files (untar the "
                         "overhangs/{sample}.tar.gz archive into it)")
    p.add_argument("--sample_id", required=True)
    p.add_argument("--output_dir", required=True,
                    help="where to write {sample}.reads.tsv (augmented) and "
                         "{sample}.cluster_coverage.tsv")
    p.add_argument("--recruit_edit_frac", type=float, default=0.10,
                    help="max allowed edits as a fraction of the shorter (query) "
                         "sequence's length when matching an ambiguous overhang "
                         "to a cluster representative; defaults to the same value "
                         "as cluster_overhangs_edlib.py's primary pass")
    args = p.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    reads_df = pd.read_csv(args.reads_tsv, sep="\t", dtype={"pos": str})
    ambiguous = read_ambiguous_overhangs(Path(args.overhang_dir))
    print(f"{len(reads_df)} clustered reads in, {len(ambiguous)} ambiguous overhangs to recruit")

    augmented, coverage = recruit(reads_df, ambiguous, args.sample_id, args.recruit_edit_frac)

    reads_path = output_dir / f"{args.sample_id}.reads.tsv"
    augmented.to_csv(reads_path, sep="\t", index=False)

    coverage_path = output_dir / f"{args.sample_id}.cluster_coverage.tsv"
    coverage.to_csv(coverage_path, sep="\t", index=False)

    n_recruited = int((augmented["recruited"] == 1).sum())
    clusters = coverage[coverage["row_type"] == "cluster"]
    unassigned = coverage[coverage["row_type"] == "unassigned"]
    print(f"Recruited {n_recruited}/{len(ambiguous)} ambiguous overhangs into "
          f"{len(clusters)} clusters")
    if not unassigned.empty:
        by_reason = unassigned.groupby("reason")["n_recruited"].sum().to_dict()
        print(f"Left unassigned: {by_reason}")
    print(f"-> {reads_path} ({len(augmented)} total reads)")
    print(f"-> {coverage_path}")


if __name__ == "__main__":
    main()
