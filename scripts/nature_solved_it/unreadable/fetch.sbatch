#!/bin/bash
# Solvers unreadable, step 2 of 4: one array task per row of proteomes.tsv downloads that
# proteome's FASTA and its UniProt cross-references (fetch_inputs.py --row).
# Submitted by submit_all.sh (array size, RUN, REPO, PROTEOMES, logs). Not meant to be run alone.
#SBATCH --job-name=unread-fetch
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
set -euo pipefail

: "${RUN:?}" "${REPO:?}" "${PROTEOMES:?}"
/usr/bin/python3 -I "$REPO/scripts/nature_solved_it/unreadable/fetch_inputs.py" \
    --row "$SLURM_ARRAY_TASK_ID" --proteomes "$PROTEOMES" --out "$RUN/inputs"
