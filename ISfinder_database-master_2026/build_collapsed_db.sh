#!/bin/bash
# Collapse the raw ISfinder collection (IS.database.fa) to a given identity
# threshold with cd-hit-est, writing IS.database.collapsed{NN}.fa plus its
# .clstr membership file.
#
# WHY THIS EXISTS: the shipped IS.database.collapsed99.fa (5804 sequences) was
# built before anything recorded how, and it is NOT reproducible from
# IS.database.fa with this command -- cd-hit-est v4.8.1 at -c 0.99 yields 5794
# representatives, and the two representative sets differ by 73/63 names. The
# original tool, version or parameters are unknown. So: anything that compares
# thresholds must regenerate BOTH with this script, or it is comparing the
# thresholds AND an unrecorded methodology change at the same time.
#
# usage:
#   ./build_collapsed_db.sh 0.99
#   ./build_collapsed_db.sh 0.95
set -euo pipefail

CD_HIT_EST="${CD_HIT_EST:-/home/hua575/miniconda3/envs/cd-hit/bin/cd-hit-est}"
IN="${IN:-IS.database.fa}"
C="${1:?usage: $0 <identity, e.g. 0.95>}"

# -n 10 is cd-hit-est's recommended word size for the 0.95-1.0 band, so the
# same value is valid for every threshold compared so far; below 0.95 it must
# drop (8-9 for 0.90-0.95, 7 for 0.88-0.90, ...) or cd-hit-est will refuse.
# -d 0 keeps full IS names in the .clstr file, which the comparison depends on.
TAG=$(echo "$C" | sed 's/^0\.//')
OUT="IS.database.collapsed${TAG}.fa"

"$CD_HIT_EST" -i "$IN" -o "$OUT" -c "$C" -n 10 -d 0 -M 0 -T "${THREADS:-4}"

echo "wrote $OUT ($(grep -c '^>' "$OUT") representatives from $(grep -c '^>' "$IN") input sequences)"
