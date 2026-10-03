#!/bin/bash
# Build IS_DB_merged.fa: the ORIGINAL (uncollapsed) ISfinder collection plus
# ISOSDB V3, in ISfinder's layout (">{name}" header, whole sequence on one
# line), then collapse it at 95% identity with exactly the cd-hit-est settings
# ISfinder_database-master_2026/build_collapsed_db.sh used for
# IS.database.collapsed95.fa, so the two collapses are comparable.
#
# ISfinder records are copied through unchanged; ISOSDB's wrapped sequence
# lines are joined onto one line. Fails if a name appears in both sources.
#
# Where a collapsed cluster mixes ISfinder and ISOSDB sequences but cd-hit chose
# an ISOSDB representative, the representative is renamed after its ISfinder
# member (see rename_reps_to_isfinder.py).
#
# Also copies TnCentral's TE.database.fa here unchanged, writes a repaired
# TE.database.clean.fa from it, and builds minibwa indexes of the collapsed IS
# database and the clean TnCentral database for read alignment.
#
# usage (from this directory):
#   THREADS=4 ./build_merged_db.sh
set -euo pipefail

MAPPING=/home/hua575/baymlab/mapping
ISFINDER="$MAPPING/ISfinder_database-master_2026/IS.database.fa"
ISOSDB="$MAPPING/pseudoR/ISOSDB.V3.fna"
TNCENTRAL="$MAPPING/TnCentral_database/TE.database.fa"
CD_HIT_EST="${CD_HIT_EST:-/home/hua575/miniconda3/envs/cd-hit/bin/cd-hit-est}"
MINIBWA="${MINIBWA:-/home/hua575/miniconda3/envs/minibwa/bin/minibwa}"
OUT=IS_DB_merged.fa
COLLAPSED=IS_DB_merged.collapsed95.fa

{
    cat "$ISFINDER"
    awk '/^>/ { if (seq != "") print seq; print; seq = ""; next }
         { seq = seq $0 }
         END { if (seq != "") print seq }' "$ISOSDB"
} > "$OUT"

dups=$(grep '^>' "$OUT" | sort | uniq -d | head)
if [ -n "$dups" ]; then
    echo "duplicate names across sources:" >&2; echo "$dups" >&2; exit 1
fi
echo "wrote $OUT ($(grep -c '^>' "$OUT") sequences: $(grep -c '^>' "$ISFINDER") ISfinder + $(grep -c '^>' "$ISOSDB") ISOSDB)"

# same flags as build_collapsed_db.sh -- see there for why -n 10 / -d 0
"$CD_HIT_EST" -i "$OUT" -o "$COLLAPSED" -c 0.95 -n 10 -d 0 -M 0 -T "${THREADS:-4}"
echo "wrote $COLLAPSED ($(grep -c '^>' "$COLLAPSED") representatives from $(grep -c '^>' "$OUT") input sequences)"

# ISOSDB-represented clusters with an ISfinder member take the ISfinder name
python3 rename_reps_to_isfinder.py "$COLLAPSED" "$COLLAPSED.clstr" "${COLLAPSED%.fa}.renamed.tsv"

cp -p "$TNCENTRAL" .
echo "copied $(basename "$TNCENTRAL")"

# TE.database.fa has glued records -- see clean_te_database.py
python3 clean_te_database.py TE.database.fa "$MAPPING/TnCentral_database/step1.TE.ID.list" TE.database.clean.fa

# minibwa indexes for the Snakefile's bwa_isfinder rule ({fasta}.l2b/.mbw).
# Two separate indexes, not one combined one: TnCentral carries IS entries and
# transposons that contain IS copies, so in a combined index reads from those
# ISs would map equally well to both databases and come out MAPQ 0.
for fa in "$COLLAPSED" TE.database.clean.fa; do
    "$MINIBWA" index -t "${THREADS:-4}" "$fa"
done
