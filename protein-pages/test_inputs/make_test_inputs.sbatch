#!/bin/bash
# Make the 50-protein Botryllus test set for build.py (see make_test_inputs.py), then run
# hmmscan against Pfam on the 50 proteins and on the Swiss-Prot proteins they match at
# E <= 1e-3, with the leaderboard's hmmscan options.
#
#   cd <folder holding make_test_inputs.py> && sbatch make_test_inputs.sbatch $SCRATCH/protein-pages-test
#
# Slurm runs a copy of this file from its spool folder, so the Python script is found from
# the folder sbatch was called in ($SLURM_SUBMIT_DIR), not from this file's own path.
#
# Sizes: the script reads about 1 GB of gzipped HMMER tables twice and holds Swiss-Prot in
# memory (about 1 GB). On 2026-10-27 it picked 16,312 proteins for hmmscan (the 50 and the
# Swiss-Prot proteins they match at E <= 1e-3).
#SBATCH -p hns,normal
#SBATCH -A ayeletv
#SBATCH -J protein-pages-test
#SBATCH -c 16
#SBATCH --mem 16G
#SBATCH -t 2:00:00
set -euo pipefail
OUT=$1
HERE=$SLURM_SUBMIT_DIR
HMMER=$SCRATCH/apptainer-cache/quay.io-biocontainers-hmmer@sha256-7a2b317b8d2fd3650b4924a8482cddeb940d4a0746c6a1501ff03ac1b7439e0c.img
PFAM=$SCRATCH/data/pfam/Pfam-A.hmm
mkdir -p "$OUT/inputs"
python3 "$HERE/make_test_inputs.py" "$OUT/inputs"
cat "$OUT/inputs/proteins.fasta" "$OUT/hmmscan_targets.fasta" > "$OUT/hmmscan_queries.fasta"
apptainer exec -B /scratch "$HMMER" hmmscan --domtblout "$OUT/inputs/hmmscan.domtblout" \
    -o /dev/null --noali -E 10.0 --domE 10.0 --cpu "$SLURM_CPUS_PER_TASK" \
    "$PFAM" "$OUT/hmmscan_queries.fasta"
ls -la "$OUT/inputs"
