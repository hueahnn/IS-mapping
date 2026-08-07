#!/usr/bin/env python3
"""Extract per-(genome, IS element) coverage statistics from BAM files aligned to ISFinder.

For every reference (IS element) with >=1 mapped read in a BAM, computes zero-filled
mean/variance of per-base depth across the full element length, and across the first
and last 200bp windows (window shrinks to min(200, length) for short elements, so the
windows may overlap for elements under 400bp).
"""
import argparse
import os
import subprocess
import sys

HEADER = [
    "genome_accession",
    "is_element",
    "is_length",
    "mean_coverage",
    "var_coverage",
    "mean_coverage_first200",
    "var_coverage_first200",
    "mean_coverage_last200",
    "var_coverage_last200",
]


def get_ref_lengths(bam_path):
    proc = subprocess.run(
        ["samtools", "view", "-H", bam_path],
        capture_output=True, text=True, check=True,
    )
    lengths = {}
    for line in proc.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        sn = ln = None
        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                sn = field[3:]
            elif field.startswith("LN:"):
                ln = int(field[3:])
        if sn is not None and ln is not None:
            lengths[sn] = ln
    return lengths


def process_bam(bam_path):
    """Return a list of (accession, ref, length, mean_full, var_full,
    mean_left200, var_left200, mean_right200, var_right200) for every
    reference with >=1 mapped read, using zero-filled full-length and
    edge-window statistics computed via running sums (no per-base array)."""
    accession = os.path.basename(bam_path)
    if accession.endswith(".bam"):
        accession = accession[:-4]

    ref_lengths = get_ref_lengths(bam_path)
    # ref -> [sum_full, sumsq_full, sum_left, sumsq_left, sum_right, sumsq_right]
    stats = {}

    proc = subprocess.Popen(
        ["samtools", "depth", bam_path],
        stdout=subprocess.PIPE, text=True,
    )
    for line in proc.stdout:
        ref, pos_s, depth_s = line.rstrip("\n").split("\t")
        length = ref_lengths.get(ref)
        if length is None:
            continue
        pos = int(pos_s)
        depth = float(depth_s)

        s = stats.get(ref)
        if s is None:
            s = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            stats[ref] = s

        s[0] += depth
        s[1] += depth * depth

        win = min(200, length)
        if pos <= win:
            s[2] += depth
            s[3] += depth * depth
        if pos > length - win:
            s[4] += depth
            s[5] += depth * depth

    proc.stdout.close()
    retcode = proc.wait()
    if retcode != 0:
        raise RuntimeError(f"samtools depth failed (exit {retcode}) for {bam_path}")

    rows = []
    for ref, s in stats.items():
        length = ref_lengths[ref]
        win = min(200, length)

        mean_full = s[0] / length
        var_full = max(s[1] / length - mean_full ** 2, 0.0)

        mean_left = s[2] / win
        var_left = max(s[3] / win - mean_left ** 2, 0.0)

        mean_right = s[4] / win
        var_right = max(s[5] / win - mean_right ** 2, 0.0)

        rows.append((
            accession, ref, length,
            mean_full, var_full,
            mean_left, var_left,
            mean_right, var_right,
        ))
    return rows


def format_row(row):
    accession, ref, length, mean_full, var_full, mean_left, var_left, mean_right, var_right = row
    return "\t".join([
        accession, ref, str(length),
        f"{mean_full:.4f}", f"{var_full:.4f}",
        f"{mean_left:.4f}", f"{var_left:.4f}",
        f"{mean_right:.4f}", f"{var_right:.4f}",
    ])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--bam", help="single BAM file to process")
    group.add_argument("--bam-list", help="file with one BAM path per line")
    ap.add_argument("--out", required=True, help="output TSV path")
    args = ap.parse_args()

    if args.bam:
        bam_paths = [args.bam]
    else:
        with open(args.bam_list) as fh:
            bam_paths = [line.strip() for line in fh if line.strip()]

    with open(args.out, "w") as out_fh:
        out_fh.write("\t".join(HEADER) + "\n")
        for bam_path in bam_paths:
            try:
                rows = process_bam(bam_path)
            except Exception as e:
                print(f"ERROR processing {bam_path}: {e}", file=sys.stderr)
                continue
            for row in rows:
                out_fh.write(format_row(row) + "\n")


if __name__ == "__main__":
    main()
