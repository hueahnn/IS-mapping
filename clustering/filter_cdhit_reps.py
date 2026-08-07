#!/usr/bin/env python3
"""Filter a cd-hit .cdhit representative-sequence fasta down to only the
representatives whose cluster has more than --min-size members, using the
matching .clstr file to determine cluster sizes and representative IDs.
"""
import argparse
import sys


def get_qualifying_reps(clstr_path, min_size):
    qualifying = set()
    rep_id = None
    count = 0

    def finalize():
        if rep_id is not None and count > min_size:
            qualifying.add(rep_id)

    with open(clstr_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                finalize()
                rep_id = None
                count = 0
                continue
            count += 1
            # member line looks like: "0\t219nt, >SOME_ID... at 1:219:101:319/+/100.00%"
            # representative line ends with "*" instead of an "at ...%" clause
            start = line.find(">")
            end = line.find("...", start)
            seq_id = line[start + 1:end]
            if line.rstrip().endswith("*"):
                rep_id = seq_id
        finalize()

    return qualifying


def filter_fasta(fasta_path, qualifying, out_path):
    written = 0
    with open(fasta_path) as fin, open(out_path, "w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                header_id = line[1:].strip()
                keep = header_id in qualifying
                if keep:
                    written += 1
            if keep:
                fout.write(line)
    return written


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cdhit-fasta", required=True, help="cd-hit representative fasta (e.g. global_left.cdhit)")
    ap.add_argument("--clstr", required=True, help="matching .clstr file")
    ap.add_argument("--min-size", type=int, default=10, help="minimum cluster size (exclusive); default 10")
    ap.add_argument("--out", required=True, help="output filtered fasta path")
    args = ap.parse_args()

    qualifying = get_qualifying_reps(args.clstr, args.min_size)
    print(f"{len(qualifying)} clusters with >{args.min_size} members", file=sys.stderr)

    written = filter_fasta(args.cdhit_fasta, qualifying, args.out)
    print(f"wrote {written} representative sequences to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
