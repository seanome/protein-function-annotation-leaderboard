#!/bin/bash
# Solvers unreadable: the laptop step once the tolerance atlas CSV (Prompt 1) exists.
# Web lookups only (UniProt and TimeTree); no search, no BUSCO, nothing heavy.
#   1. atlas species -> UniProt reference proteomes and TimeTree times (atlas_proteomes.csv)
#   2. random controls, seed 20261019, matched in 100-My bins (controls.csv, draw_record.json)
#   3. the run list (proteomes.tsv) and its resource prediction (resources_run.tsv)
# When a bin runs short, every candidate gets a TimeTree lookup: about 5 minutes for all 3,513
# eukaryotic candidates of 2026_03 (stand-in test, 2026-10-19).
#
#   bash scripts/nature_solved_it/unreadable/after_atlas.sh <atlas.csv> [species column, default "species"]
#
# Then commit the files it lists, push, and on Sherlock run submit_all.sh (see README.md).
set -euo pipefail

ATLAS=${1:?atlas CSV from Prompt 1}
COL=${2:-species}
REPO=$(cd "$(dirname "$0")/../../.." && pwd)
HERE=$REPO/scripts/nature_solved_it/unreadable
RES=$REPO/results/nature-solved-it/unreadable

python3 -I "$HERE/proteome_table.py" --atlas "$ATLAS" --species-column "$COL" --out "$RES/atlas_proteomes.csv"
python3 -I "$HERE/draw_controls.py" --atlas "$RES/atlas_proteomes.csv" --models "$RES/model_organisms.csv" --out-dir "$RES"
python3 -I "$HERE/predict_resources.py" --proteomes "$RES/proteomes.tsv" --sacct "$RES/sacct_47961379.txt" --out "$RES/resources_run.tsv"
echo "Files to commit:"
for f in atlas_proteomes.csv controls.csv draw_record.json proteomes.tsv resources_run.tsv timetree_cache.json reference_proteomes_README.txt.gz; do
    echo "  results/nature-solved-it/unreadable/$f"
done
