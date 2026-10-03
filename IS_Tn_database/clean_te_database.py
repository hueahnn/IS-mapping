#!/usr/bin/env python
# Write TE.database.clean.fa: TnCentral's TE.database.fa repaired and laid out
# like ISfinder's (">{name}" header, whole uppercase sequence on one line), so
# it can be indexed for alignment.
#
# TE.database.fa has 20 records whose previous sequence never got a trailing
# newline AND whose own header has none either, so a line reads
# "...seqA>In1223-KX784502seqB" (one line runs 4 records together). Read as-is,
# those 20 elements vanish into their neighbours and the ">" lands inside a
# sequence. They are split using the names in step1.TE.ID.list (the URLs
# TnCentral.step1.py scraped), longest match first, rather than guessing where
# an accession ends and the sequence starts. Headers that start their own line
# are kept as written, including "Tn7246-" (accession missing in
# TE.database.fa; step1.TE.ID.list has Tn7246-EU370913.1).
#
# usage:
#   python clean_te_database.py TE.database.fa step1.TE.ID.list TE.database.clean.fa
import sys


def main():
    raw, id_list, out = sys.argv[1:4]
    known = {l.strip().rstrip("/").rsplit("/", 1)[-1] for l in open(id_list) if l.strip()}

    lines = []
    for line in open(raw):
        line = line.rstrip("\n\r")
        i = line.find(">", 1)
        while i > 0:
            lines.append(line[:i])
            body = line[i + 1:]
            name = max((n for n in known if body.startswith(n)), key=len, default=None)
            if name is None:
                sys.exit(f"{raw}: can't split header from sequence: {body[:60]}")
            lines.append(">" + name)
            line = body[len(name):]
            i = line.find(">")
        lines.append(line)

    records, name, seq = [], None, []
    for line in lines:
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq)))

    names = [n for n, _ in records]
    if len(set(names)) != len(names):
        sys.exit("duplicate names after splitting")
    with open(out, "w") as fh:
        for n, s in records:
            if not s:
                print(f"WARNING: {n} has no sequence, skipping", file=sys.stderr)
                continue
            fh.write(f">{n}\n{s.upper()}\n")
    print(f"wrote {out} ({len(records)} sequences)")


if __name__ == "__main__":
    main()
