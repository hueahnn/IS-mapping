#!/usr/bin/env python
# Give ISfinder names to collapsed clusters that cd-hit-est represented with an
# ISOSDB sequence. ISOSDB names ("ISOSDB123") say nothing about the element,
# while ISfinder's do, so whenever a cluster holds at least one ISfinder member
# its representative is renamed after one of them. The representative SEQUENCE
# is left as cd-hit chose it (its longest member); only the header changes.
#
# With several ISfinder members, the name comes from the one with the highest
# identity to the representative, then the longest, then the first listed.
#
# The .clstr file is left as cd-hit wrote it (original ISOSDB names); the
# rename table records new_name -> original representative for going between
# the two.
#
# usage:
#   python rename_reps_to_isfinder.py IS_DB_merged.collapsed95.fa \
#       IS_DB_merged.collapsed95.fa.clstr IS_DB_merged.collapsed95.renamed.tsv
import re
import sys

MEMBER_RE = re.compile(r"\t(\d+)nt, >(\S+?)\.\.\. (?:(\*)|at (?:[+-]/)?([\d.]+)%)")


def is_isfinder(name):
    return not name.startswith("ISOSDB")


def parse_clstr(path):
    clusters, cur = [], None
    for line in open(path):
        if line.startswith(">Cluster"):
            cur = []
            clusters.append(cur)
            continue
        m = MEMBER_RE.search(line)
        if not m:
            sys.exit(f"unrecognized .clstr line: {line!r}")
        length, name, star, ident = m.groups()
        cur.append((name, int(length), 100.0 if star else float(ident), bool(star)))
    return clusters


def main():
    fasta, clstr, table = sys.argv[1:4]

    renames = {}
    with open(table, "w") as out:
        out.write("new_name\toriginal_representative\tn_members\tisfinder_members\n")
        for members in parse_clstr(clstr):
            rep = next(name for name, _, _, star in members if star)
            isf = [m for m in members if is_isfinder(m[0])]
            if is_isfinder(rep) or not isf:
                continue
            best = max(isf, key=lambda m: (m[2], m[1]))  # max keeps the first on ties
            renames[rep] = best[0]
            out.write(f"{best[0]}\t{rep}\t{len(members)}\t{','.join(m[0] for m in isf)}\n")

    lines = open(fasta).read().splitlines()
    with open(fasta, "w") as fh:
        for line in lines:
            if line.startswith(">"):
                line = ">" + renames.pop(line[1:].split()[0], line[1:])
            fh.write(line + "\n")
    if renames:
        sys.exit(f"{len(renames)} representatives from {clstr} not found in {fasta}")
    print(f"renamed {sum(1 for _ in open(table)) - 1} ISOSDB representatives to ISfinder names -> {table}")


if __name__ == "__main__":
    main()
