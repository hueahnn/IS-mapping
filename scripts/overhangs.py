# purpose: four-stage pipeline for soft-clipped (overhang) reads from IS-mapped BAM.
# usage:
#   python overhangs.py ERR9860257.bam output_dir/ --sample_id ERR9860257 \
#       --min_pident 90 --boundary_tolerance 1 --min_is_hits 10 --min_mapq 30 \
#       --min_is_coverage_depth 10 --min_is_coverage_fraction 0.9
#
# region: pipeline steps outline
# STAGE 0 (filter_by_is_coverage): cohort-level pre-filter run on the RAW input BAM,
# before any per-read filtering. For each IS element, computes per-position depth
# from all primary mapped reads and keeps only reads belonging to IS elements where
# >= min_is_coverage_fraction of reference positions have depth >= min_is_coverage_depth.
# This is a presence/sequencing-depth QC on the IS element as a whole -- it does NOT
# consider junction/boundary quality (that's what stages 1-4 are for). Writes an
# unmodified BAM (no tags) of reads mapped to qualifying IS elements only.
#
# STAGE 1 (filter_bam): applies the per-read quality filter and writes a SINGLE output
# BAM containing the full, unmodified original alignment records for reads that pass
# (same header/reference set as input — sortable/indexable with standard samtools
# tools). Each surviving read is stamped with custom tags recording the filtering
# decision, since these are needed by later stages and have no standard SAM field:
#   ZI:f  percent identity, (aligned_length - NM) / aligned_length * 100
#   ZL:i  1 if this read qualifies as a LEFT overhang, else 0
#   ZR:i  1 if this read qualifies as a RIGHT overhang, else 0
#   ZA:i  1 if this read's MAPQ is below min_mapq ("ambiguous mapping"), else 0
# (both ZL and ZR can be 1 for the same read if the IS element is shorter than the
# read's aligned span; a read can furnish both a left and a right overhang.)
# ZA does NOT gate whether the read is written here -- unlike every other check in
# this stage, a low-MAPQ read still passes filter_bam. It only changes which FASTA
# bucket stage 3 routes the read's overhang(s) into (see STAGE 3 below).
#
# filtering order (each stage gates what reaches the next):
#   -1. sequence content filter — the read's full query_sequence must contain only
#       A/C/G/T (case-insensitive); a read with an N or any other ambiguity code is
#       dropped outright, before any boundary/alignment checks, since an ambiguous
#       base makes the read unusable regardless of where it aligns.
#   0. minimum IS-reference alignment span — reference_end - reference_start must be
#      >= min_is_alignment_bases (default 20). This is how much of the IS ELEMENT
#      itself the alignment covers, not how much of the read is aligned; a read with
#      only a sliver of reference-anchored alignment makes reference_start/end
#      unreliable for the boundary check that follows, so this runs first.
#   1. positional index filter — read's alignment must sit near one end of the IS
#      reference (within --boundary_tolerance bp of position 0 or of the IS length).
#      This determines which side (left/right) is even eligible for that read; a read
#      aligned to the middle of the IS element cannot produce a real boundary overhang
#      even if it happens to be soft-clipped there.
#   1.5. dual-clip filter — reject reads soft-clipped on BOTH ends of the READ itself
#        (raw CIGAR, checked before gating by eligible side). A read clipped on both
#        ends isn't clean single-junction evidence, regardless of which end happened
#        to be near a boundary.
#   2. softclip presence — on the side made eligible by (1), is there actually a
#      soft-clip in the CIGAR?
#   3. overhang length filter — is the soft-clip >= min_clip_len (default 5)?
#   4. percent identity filter — (aligned_length - NM) / aligned_length * 100 >= min_pident
#
# STAGE 2 (filter_by_is_hit_count): separate cohort-level step, run on stage 1's
# output. Drops every read belonging to an IS element with fewer than min_is_hits
# surviving reads. Kept as its own function/output file (not folded into filter_bam)
# so it can be rerun with a different threshold, or skipped, independent of the
# per-read filtering — and so filter_bam can stay a simple single-pass function.
#
# STAGE 3 (extract_overhangs_from_bam): reads stage 2's output and, for each IS
# element, writes its left-side clipped subsequences into one FASTA, its right-side
# clipped subsequences into another, and -- for any read flagged ZA=1 (ambiguous
# MAPQ) -- its clipped subsequence(s) into a third, per-IS "ambiguous" FASTA instead
# of the left/right one it would otherwise land in. This keeps low-confidence
# mapping evidence available for inspection rather than discarding it outright.
# Reads the ZL/ZR/ZA tags from stage 1 rather than recomputing boundary/length/mapq
# eligibility — avoids having multiple functions independently re-derive the same
# decision from separately-passed parameters that could drift out of sync.
#
# NOTE on stage 0 vs stage 2: filter_by_is_coverage (stage 0) and
# filter_by_is_hit_count (stage 2) are NOT interchangeable despite both being
# cohort-level, per-IS-element filters. Stage 0 runs on the RAW BAM and asks "does
# this IS element have real sequencing depth at all, independent of junction
# quality." Stage 2 runs on STAGE-1 OUTPUT (post junction-quality filtering) and
# asks "does this IS element have enough surviving junction-quality reads, per side."
# Do not reorder these relative to filter_bam.
# endregion


import pysam
import shutil
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter


def percent_identity(read):
    """(aligned_length - NM) / aligned_length * 100. Returns None if NM tag is absent
    or aligned_length is 0."""
    if not read.has_tag('NM'):
        return None
    aligned_length = read.query_alignment_length
    if not aligned_length:
        return None
    nm = read.get_tag('NM')
    return 100.0 * (aligned_length - nm) / aligned_length


def is_pure_acgt(seq: str) -> bool:
    """True if seq contains only A/C/G/T (case-insensitive) -- False if it has an N
    or any other ambiguity code."""
    return set(seq.upper()) <= {"A", "C", "G", "T"}


def filter_by_is_coverage(bam_path: str, output_bam_path: str,
                           min_depth: int = 10, min_coverage_fraction: float = 0.9):
    """
    Stage 0: cohort-level pre-filter, run on the RAW input BAM before filter_bam.
    For each IS element, computes per-position depth across ALL primary mapped reads
    and keeps only IS elements where >= min_coverage_fraction of reference positions
    have depth >= min_depth. Writes a BAM containing only reads mapped to qualifying
    IS elements (unmodified records — no tags added, since this isn't a per-read
    quality decision).

    ASSUMPTION: depth is computed from all primary mapped reads regardless of
    whether they'd later pass filter_bam's junction filters (stages 0-4). This is
    "does this element have real sequencing depth at all" — a presence/absence QC —
    not "does it have depth among junction-quality reads." If you want the latter,
    this needs to run AFTER filter_bam instead of before it.

    Only M/=/X (aligned, reference-consuming) blocks count toward depth — a
    deletion spanned by a read doesn't count as that position being "covered".
    Secondary/supplementary alignments are excluded, consistent with filter_bam.

    Two passes over bam_path (depth first, then write) — same reason as
    filter_by_is_hit_count: total coverage per IS element isn't known until every
    read is counted.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        is_lengths = dict(zip(bam.header.references, bam.header.lengths))
        depth = {name: np.zeros(length, dtype=np.int32) for name, length in is_lengths.items()}

        # pass 1: accumulate depth per IS element
        reads = list(bam.fetch(until_eof=True))
        for read in reads:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_name is None:
                continue
            for start, end in read.get_blocks():  # M/=/X blocks only; deletions don't count as covered
                depth[read.reference_name][start:end] += 1

        coverage_frac = {}
        qualifying_is = set()
        for name, arr in depth.items():
            if len(arr) == 0:
                continue
            frac = float((arr >= min_depth).sum()) / len(arr)
            coverage_frac[name] = frac
            if frac >= min_coverage_fraction:
                qualifying_is.add(name)

        # pass 2: write only reads mapped to qualifying IS elements
        written = 0
        out_bam = pysam.AlignmentFile(output_bam_path, "wb", template=bam)
        for read in reads:
            if read.reference_name in qualifying_is:
                out_bam.write(read)
                written += 1
        out_bam.close()

    n_total = len(depth)
    n_pass = len(qualifying_is)
    print(f"IS coverage filter: {n_pass}/{n_total} IS elements have "
          f">= {min_coverage_fraction*100:.0f}% of positions at >= {min_depth}x depth")
    print(f"Wrote {written} records ({n_pass} distinct IS elements) to {output_bam_path}")
    print(f"{len(qualifying_is)} IS elements above coverage threshold: {qualifying_is}")

    return qualifying_is


def filter_bam(bam_path: str, filtered_bam_path: str,
               min_pident: float = 90.0, boundary_tolerance: int = 20,
               min_is_alignment_bases: int = 20, min_clip_len: int = 5,
               min_mapq: int = 30):
    """
    Stage 1: write a single BAM of reads passing all per-read filter stages. Each
    surviving read is tagged with ZI (pident), ZL (1/0 left-eligible), ZR (1/0
    right-eligible), ZA (1/0 ambiguous mapping). Does not consider cohort-level
    properties like per-IS-element hit count or reference coverage — see
    filter_by_is_coverage() (stage 0, run before this) and filter_by_is_hit_count()
    (stage 2, run after this).

    min_mapq: reads with MAPQ below this are NOT dropped -- they still pass every
    other check the same as any other read, but are stamped ZA=1 so stage 3 can
    route their overhang(s) to a separate "ambiguous" FASTA instead of left/right.
    """
    written = 0
    ambiguous = 0
    is_seen = set()
    dropped = {"non_acgt": 0, "short_is_alignment": 0, "not_at_boundary": 0,
               "dual_clip": 0, "no_softclip_on_eligible_side": 0, "short_clip": 0,
               "no_nm": 0, "low_pident": 0}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # reference (IS element) lengths, built once — avoids repeated header lookups
        # per read across ~20k+ reads per sample
        is_lengths = dict(zip(bam.header.references, bam.header.lengths))

        # template=bam copies the exact header (references, lengths, @PG, etc.) from
        # the input BAM, so output records stay valid alignments against the same
        # IS reference set.
        out_bam = pysam.AlignmentFile(filtered_bam_path, "wb", template=bam)

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.cigartuples is None:
                continue

            seq = read.query_sequence
            cigar = read.cigartuples

            if seq is None:
                continue

            # --- stage -1: sequence content filter ---
            # reject reads with an N or any other non-ACGT base anywhere in the full
            # read, before any boundary/alignment checks -- an ambiguous base makes
            # the read unusable regardless of where it aligns.
            if not is_pure_acgt(seq):
                dropped["non_acgt"] += 1
                continue

            # --- stage 0: minimum IS-reference alignment span ---
            # how many bp of the IS ELEMENT (not the read) this alignment actually
            # covers — reference_end - reference_start. A read with only a sliver of
            # reference-anchored alignment (e.g. 3bp before a long softclip) makes
            # reference_start/reference_end themselves unreliable to use for boundary
            # classification, so this is checked before the positional filter, not after.
            is_alignment_bases = read.reference_end - read.reference_start
            if is_alignment_bases < min_is_alignment_bases:
                dropped["short_is_alignment"] += 1
                continue

            # --- stage 1: positional index filter ---
            is_length = is_lengths.get(read.reference_name)
            if is_length is None:
                continue

            near_start = read.reference_start <= boundary_tolerance
            near_end = (is_length - read.reference_end) <= boundary_tolerance

            if not (near_start or near_end):
                dropped["not_at_boundary"] += 1
                continue

            # --- stage 1.5: dual-clip filter ---
            # raw (ungated) softclip lengths on BOTH ends of the read, regardless of
            # which side the positional filter made eligible. A read clipped on both
            # ends is ambiguous evidence — not a clean single-junction overhang — so
            # it's tossed outright rather than being allowed through on whichever
            # side happened to be near a boundary.
            raw_left_clip = cigar[0][1] if cigar[0][0] == 4 else 0
            raw_right_clip = cigar[-1][1] if cigar[-1][0] == 4 else 0

            if raw_left_clip > 0 and raw_right_clip > 0:
                dropped["dual_clip"] += 1
                continue

            # --- stage 2: softclip presence, gated to the side made eligible above ---
            left_clip = raw_left_clip if near_start else 0
            right_clip = raw_right_clip if near_end else 0

            if left_clip == 0 and right_clip == 0:
                dropped["no_softclip_on_eligible_side"] += 1
                continue

            # --- stage 3: overhang length filter ---
            left_clip_ok = left_clip >= min_clip_len
            right_clip_ok = right_clip >= min_clip_len
            if not (left_clip_ok or right_clip_ok):
                dropped["short_clip"] += 1
                continue

            # --- stage 4: percent identity filter ---
            pident = percent_identity(read)
            if pident is None:
                dropped["no_nm"] += 1
                continue
            if pident < min_pident:
                dropped["low_pident"] += 1
                continue

            # --- MAPQ check: tags but does NOT gate -- an ambiguously-mapped read
            # still passes filter_bam, it's just flagged for stage 3 to route into
            # the "ambiguous" FASTA instead of left/right.
            is_ambiguous = read.mapping_quality < min_mapq

            read.set_tag("ZI", float(pident))
            read.set_tag("ZL", int(left_clip_ok))
            read.set_tag("ZR", int(right_clip_ok))
            read.set_tag("ZA", int(is_ambiguous))
            out_bam.write(read)
            is_seen.add(read.reference_name)
            written += 1
            if is_ambiguous:
                ambiguous += 1

        out_bam.close()

    print(f"Wrote {written} filtered records ({len(is_seen)} distinct IS elements) "
          f"to {filtered_bam_path}, of which {ambiguous} flagged ambiguous (MAPQ < {min_mapq})")
    print(f"Dropped: {dropped['non_acgt']} (non-ACGT base in read), "
          f"{dropped['short_is_alignment']} (IS-reference alignment span < {min_is_alignment_bases}bp), "
          f"{dropped['not_at_boundary']} (not near IS boundary), "
          f"{dropped['dual_clip']} (clipped on both ends), "
          f"{dropped['no_softclip_on_eligible_side']} (no softclip on eligible side), "
          f"{dropped['short_clip']} (clip < {min_clip_len}), "
          f"{dropped['no_nm']} (no NM tag), "
          f"{dropped['low_pident']} (pident < {min_pident})")


def filter_by_is_hit_count(bam_path: str, output_bam_path: str, min_is_hits: int = 10):
    """
    Stage 2: separate cohort-level step, run on stage 1's output. For each IS
    element, LEFT and RIGHT surviving-hit counts are evaluated independently
    against min_is_hits — an IS element can keep its left overhangs while losing
    its right overhangs (or vice versa) if only one side clears the threshold.

    A read tagged both ZL=1 and ZR=1 (qualifies for both sides) is NOT dropped
    outright if only one side fails its threshold — its tag for the failing side is
    cleared (set to 0) and it is still written if at least one side survives. It is
    only dropped entirely if both sides fail (or it only ever had one side, and that
    side fails).

    Reads flagged ZA=1 (ambiguous MAPQ, from filter_bam) are excluded from the
    left/right hit counts computed below -- they are lower-confidence evidence and
    shouldn't count toward whether an IS element has "enough" confident junction
    support. They also bypass the resulting min_is_hits gate entirely and are always
    written through unchanged: they're a separate bucket handled at stage 3, not
    part of the confident-evidence pool this filter protects.

    Two passes over bam_path are unavoidable — total per-side hits per IS element
    aren't known until every read has been counted — but this is isolated to its own
    function/output file, so it can be rerun with a different threshold, or skipped,
    without touching the per-read filtering in filter_bam().
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        reads = list(bam.fetch(until_eof=True))  # pass 1: load + count

        left_counts = Counter(r.reference_name for r in reads
                               if r.get_tag("ZL") == 1 and r.get_tag("ZA") == 0)
        right_counts = Counter(r.reference_name for r in reads
                                if r.get_tag("ZR") == 1 and r.get_tag("ZA") == 0)

        left_dropped_is = {name for name, n in left_counts.items() if n < min_is_hits}
        right_dropped_is = {name for name, n in right_counts.items() if n < min_is_hits}

        written = 0
        is_seen = set()
        out_bam = pysam.AlignmentFile(output_bam_path, "wb", template=bam)
        for read in reads:  # pass 2: write survivors, clearing tags per-side as needed
            if read.get_tag("ZA") == 1:
                # ambiguous-mapping reads bypass the confident-hit-count gate above
                out_bam.write(read)
                is_seen.add(read.reference_name)
                written += 1
                continue

            zl = read.get_tag("ZL")
            zr = read.get_tag("ZR")

            if zl == 1 and read.reference_name in left_dropped_is:
                zl = 0
            if zr == 1 and read.reference_name in right_dropped_is:
                zr = 0

            if zl == 0 and zr == 0:
                continue  # neither side survives for this read's IS element

            read.set_tag("ZL", int(zl))
            read.set_tag("ZR", int(zr))
            out_bam.write(read)
            is_seen.add(read.reference_name)
            written += 1
        out_bam.close()

    print(f"Wrote {written} records ({len(is_seen)} distinct IS elements) to {output_bam_path}")
    if left_dropped_is:
        print(f"Excluded LEFT overhangs for {len(left_dropped_is)} IS elements with "
              f"< {min_is_hits} left hits: {sorted(left_dropped_is)}")
    if right_dropped_is:
        print(f"Excluded RIGHT overhangs for {len(right_dropped_is)} IS elements with "
              f"< {min_is_hits} right hits: {sorted(right_dropped_is)}")


def extract_overhangs_from_bam(filtered_bam_path: str, output_dir: str, sample_id: str):
    """
    Stage 3: read the filtered BAM and, per IS element, write left overhangs to one
    FASTA, right overhangs to another, and ambiguous-MAPQ overhangs (ZA=1, from
    either side) to a third. Uses the ZL/ZR/ZI/ZA tags stamped by filter_bam — does
    not recompute filtering decisions.

    A read flagged ZA=1 still has its clipped sequence(s) extracted exactly as it
    would if confident (per ZL/ZR), it just lands in the "ambiguous" bucket instead
    of "left"/"right" for that IS element. The header's side:left/right field still
    records which end of the read was clipped -- only the output bucket changes.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    is_overhangs = defaultdict(list)  # (is_name, bucket) -> [(header, seq), ...]
    written = {"left": 0, "right": 0, "ambiguous": 0}

    with pysam.AlignmentFile(filtered_bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not (read.has_tag("ZL") and read.has_tag("ZR") and read.has_tag("ZI")
                    and read.has_tag("ZA")):
                # not a record produced by filter_bam — skip rather than guess
                continue

            seq = read.query_sequence
            cigar = read.cigartuples
            if seq is None or cigar is None:
                continue

            strand = "-" if read.is_reverse else "+"
            pident = read.get_tag("ZI")
            ambiguous = read.get_tag("ZA") == 1

            if read.get_tag("ZL") == 1:
                left_clip = cigar[0][1]
                clipped_seq = seq[:left_clip]
                bucket = "ambiguous" if ambiguous else "left"
                header = (f">{sample_id}|{read.query_name}|{read.reference_name}"
                          f"|pos:{read.reference_start}|side:left|strand:{strand}"
                          f"|pident:{pident:.2f}|mapq:{read.mapping_quality}")
                is_overhangs[(read.reference_name, bucket)].append((header, clipped_seq))
                written[bucket] += 1

            if read.get_tag("ZR") == 1:
                right_clip = cigar[-1][1]
                clipped_seq = seq[len(seq) - right_clip:]
                bucket = "ambiguous" if ambiguous else "right"
                header = (f">{sample_id}|{read.query_name}|{read.reference_name}"
                          f"|pos:{read.reference_end}|side:right|strand:{strand}"
                          f"|pident:{pident:.2f}|mapq:{read.mapping_quality}")
                is_overhangs[(read.reference_name, bucket)].append((header, clipped_seq))
                written[bucket] += 1

    manifest_path = output_dir / f"{sample_id}.manifest.tsv"
    with open(manifest_path, "w") as mf:
        mf.write("is_element\tside\tn_seqs\tfasta_path\n")
        for (is_name, side), records in is_overhangs.items():
            safe_name = is_name.replace("/", "_")
            out_path = output_dir / f"{sample_id}__{safe_name}__{side}.overhangs.fasta"
            with open(out_path, "w") as fout:
                for header, seq in records:
                    fout.write(f"{header}\n{seq}\n")
            mf.write(f"{is_name}\t{side}\t{len(records)}\t{out_path}\n")

    n_is = len({is_name for is_name, side in is_overhangs.keys()})
    print(f"Wrote {written['left']} left, {written['right']} right, and "
          f"{written['ambiguous']} ambiguous overhang sequences across {n_is} IS "
          f"elements to {output_dir}")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("bam_path")
    p.add_argument("output_dir")
    p.add_argument("--sample_id", required=True)
    p.add_argument("--min_clip_len", type=int, default=5,
                    help="minimum overhang (softclip) length in bp")
    p.add_argument("--min_pident", type=float, default=90.0,
                    help="minimum (aligned_length - NM) / aligned_length * 100")
    p.add_argument("--min_mapq", type=int, default=30,
                    help="reads below this MAPQ are not dropped, but are routed to "
                         "the 'ambiguous' FASTA instead of left/right at stage 3")
    p.add_argument("--boundary_tolerance", type=int, default=20,
                    help="max bp a read's alignment can sit from position 0 or the IS "
                         "reference's full length to still count as 'at the boundary'")
    p.add_argument("--min_is_alignment_bases", type=int, default=20,
                    help="minimum bp of the IS reference (not the read) this "
                         "alignment must span (reference_end - reference_start)")
    p.add_argument("--min_is_hits", type=int, default=10,
                    help="minimum surviving reads an IS element must have to be kept "
                         "in the output at all (evaluated independently per side)")
    p.add_argument("--min_is_coverage_depth", type=int, default=10,
                    help="minimum per-position depth for the IS reference coverage "
                         "pre-filter (stage 0)")
    p.add_argument("--min_is_coverage_fraction", type=float, default=0.9,
                    help="minimum fraction of IS reference positions that must meet "
                         "min_is_coverage_depth for that IS element to be kept "
                         "(stage 0)")
    p.add_argument("--filter_only", action="store_true",
                    help="run stage 1 only (per-read filter), skip IS-hit-count "
                         "filtering and FASTA extraction")
    p.add_argument("--skip_is_hit_filter", action="store_true",
                    help="run stages 0, 1 and 3, skip the per-IS-element hit-count "
                         "filter (stage 2)")
    p.add_argument("--skip_is_coverage_filter", action="store_true",
                    help="skip the IS reference coverage pre-filter (stage 0) and "
                         "run filter_bam directly on the raw input BAM")
    args = p.parse_args()

    output_dir = Path(args.output_dir)
    if output_dir.exists():
        print(f"Clearing existing output directory: {output_dir}")
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # stage 0: IS reference coverage pre-filter (raw BAM in, separate file, separate
    # function — same rationale as stage 2: rerunnable/skippable independent of the
    # per-read filtering in filter_bam())
    if args.skip_is_coverage_filter:
        bam_for_filter_bam = args.bam_path
    else:
        coverage_filtered_bam_path = output_dir / f"{args.sample_id}.coverage_filtered.bam"
        filter_by_is_coverage(args.bam_path, str(coverage_filtered_bam_path),
                               min_depth=args.min_is_coverage_depth,
                               min_coverage_fraction=args.min_is_coverage_fraction)
        bam_for_filter_bam = str(coverage_filtered_bam_path)

    filtered_bam_path = output_dir / f"{args.sample_id}.filtered.bam"

    # stage 1: per-read filter (keyword args — avoids silent argument-order bugs
    # when filter_bam's signature changes)
    filter_bam(bam_for_filter_bam, str(filtered_bam_path),
               min_pident=args.min_pident,
               boundary_tolerance=args.boundary_tolerance,
               min_is_alignment_bases=args.min_is_alignment_bases,
               min_clip_len=args.min_clip_len,
               min_mapq=args.min_mapq)

    if args.filter_only:
        raise SystemExit(0)

    # stage 2: per-IS-element hit count filter (separate file, separate function)
    if args.skip_is_hit_filter:
        final_bam_path = filtered_bam_path
    else:
        final_bam_path = output_dir / f"{args.sample_id}.filtered.min_hits.bam"
        filter_by_is_hit_count(str(filtered_bam_path), str(final_bam_path), args.min_is_hits)

    # stage 3: FASTA extraction
    extract_overhangs_from_bam(str(final_bam_path), str(output_dir), args.sample_id)