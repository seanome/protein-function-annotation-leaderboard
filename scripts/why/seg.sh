#!/bin/bash
# SEG low-complexity masking of the human and oyster one-protein-per-gene sets, for the
# seg_fraction column of labels.py. BLAST+ 2.16.0 segmasker, SEG defaults (window 12,
# locut 2.2, hicut 2.5), the same rule as the site's low_complexity axis. Seconds per proteome.
#
#   bash scripts/why/seg.sh ~/data/why-numbers/2026_03
# Writes <raw>/seg/{human,oyster}.segmasker.txt.gz and prints their sha256 for sources.json
# (keys human_segmasker and oyster_segmasker, release "BLAST+ 2.16.0 segmasker, defaults").
set -euo pipefail
RAW=${1:?raw folder, e.g. ~/data/why-numbers/2026_03}
mkdir -p "$RAW/seg"
segmasker -version-full | head -1 | tee "$RAW/seg/segmasker_version.txt"
grep -q "2.16.0" "$RAW/seg/segmasker_version.txt" || { echo "need BLAST+ 2.16.0" >&2; exit 1; }
for pair in human:UP000005640_9606 oyster:UP000005408_29159; do
    tag=${pair%%:*} stem=${pair#*:}
    gunzip -c "$RAW/$stem.fasta.gz" | segmasker -outfmt interval | gzip -n > "$RAW/seg/$tag.segmasker.txt.gz"
done
# The parser in labels.py reads '>header' lines then 'start - end' lines. Show the first
# record so a change in the format is seen before labels.py reads it.
gunzip -c "$RAW/seg/oyster.segmasker.txt.gz" | head -5
sha256sum "$RAW/seg/"*.segmasker.txt.gz
