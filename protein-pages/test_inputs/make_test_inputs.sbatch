#!/bin/bash
# Make the 50-protein Botryllus test set for build.py (see make_test_inputs.py), then run
# hmmscan against Pfam on the 50 proteins and on the Swiss-Prot proteins they match at
# E <= 10, with the leaderboard's hmmscan options.
#
#   cd <folder holding make_test_inputs.py> && sbatch make_test_inputs.sbatch $SCRATCH/protein-pages-test
#
# An optional second argument is an earlier run's hmmscan.domtblout. Its rows are kept for the
# proteins it already covers and hmmscan runs only on the rest. hmmscan scores each protein on
# its own against the same Pfam file, so the kept rows are what a new run would write.
#
# Slurm runs a copy of this file from its spool folder, so the Python script is found from
# the folder sbatch was called in ($SLURM_SUBMIT_DIR), not from this file's own path.
#
# Sizes: the script reads about 1 GB of gzipped HMMER tables twice and holds Swiss-Prot in
# memory (about 1 GB). On 2026-10-27 it picked 16,312 proteins for hmmscan (the 50 and the
# Swiss-Prot proteins they match at E <= 1e-3; later runs added BHF and the targets up to E 10).
#SBATCH -p hns,normal
#SBATCH -A ayeletv
#SBATCH -J protein-pages-test
#SBATCH -c 16
#SBATCH --mem 16G
#SBATCH -t 2:00:00
set -euo pipefail
OUT=$1
PREV=${2:-}
HERE=$SLURM_SUBMIT_DIR
HMMER=$SCRATCH/apptainer-cache/quay.io-biocontainers-hmmer@sha256-7a2b317b8d2fd3650b4924a8482cddeb940d4a0746c6a1501ff03ac1b7439e0c.img
PFAM=$SCRATCH/data/pfam/Pfam-A.hmm
mkdir -p "$OUT/inputs"
python3 "$HERE/make_test_inputs.py" "$OUT/inputs"
cat "$OUT/inputs/proteins.fasta" "$OUT/hmmscan_targets.fasta" > "$OUT/hmmscan_queries.fasta"
if [ -n "$PREV" ]; then
    # names hmmscan was asked about before: column 4 of every row, plus the queries file of that run
    grep -v '^#' "$PREV" | awk '{print $4}' | sort -u > "$OUT/prev_names.txt"
    if [ -f "$(dirname "$PREV")/../hmmscan_queries.fasta" ]; then
        grep '^>' "$(dirname "$PREV")/../hmmscan_queries.fasta" | cut -c2- | awk '{print $1}' >> "$OUT/prev_names.txt"
        sort -u -o "$OUT/prev_names.txt" "$OUT/prev_names.txt"
    fi
    awk 'NR==FNR {seen[$1]=1; next} /^>/ {keep = !(substr($1,2) in seen)} keep' \
        "$OUT/prev_names.txt" "$OUT/hmmscan_queries.fasta" > "$OUT/hmmscan_new.fasta"
    grep -c '^>' "$OUT/hmmscan_new.fasta" || true
    awk 'NR==FNR {if (/^>/) want[substr($1,2)]=1; next} /^#/ || ($4 in want)' \
        "$OUT/hmmscan_queries.fasta" "$PREV" > "$OUT/hmmscan_kept.domtblout"
    apptainer exec -B /scratch "$HMMER" hmmscan --domtblout "$OUT/hmmscan_new.domtblout" \
        -o /dev/null --noali -E 10.0 --domE 10.0 --cpu "$SLURM_CPUS_PER_TASK" \
        "$PFAM" "$OUT/hmmscan_new.fasta"
    cat "$OUT/hmmscan_kept.domtblout" <(grep -v '^#' "$OUT/hmmscan_new.domtblout") > "$OUT/inputs/hmmscan.domtblout"
else
    apptainer exec -B /scratch "$HMMER" hmmscan --domtblout "$OUT/inputs/hmmscan.domtblout" \
        -o /dev/null --noali -E 10.0 --domE 10.0 --cpu "$SLURM_CPUS_PER_TASK" \
        "$PFAM" "$OUT/hmmscan_queries.fasta"
fi
ls -la "$OUT/inputs"
