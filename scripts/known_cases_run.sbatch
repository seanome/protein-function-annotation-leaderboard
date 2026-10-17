#!/bin/bash
# Search the 24 known-case query proteins and their 24 decoys against Swiss-Prot 2026_03 with
# the nine tools of Table 1, each at the setting Table 1 calls its best.
#
# Run on Sherlock, never on a laptop:
#   sbatch --export=ALL,KC_REPO=$PWD scripts/known_cases_run.sbatch
#
# Input: data/known_cases/queries.fasta, written by scripts/known_cases.py.
# Output: $OUT/<tool>.tsv plus $OUT/resources.tsv (core-hours and peak memory per step, read
# from sacct), which scripts/known_cases_results.py turns into the page's numbers.
#
# What each step was sized from (measured, not guessed). The numbers come from the Sherlock run
# of nextflow-runs/convergence-afp in seanome/2024-kmerseek-analysis, whose queries were 20
# proteins against this same Swiss-Prot release
# ($SCRATCH/2024-kmerseek-convergence-afp/nextflow-runs/convergence-afp/runs/convergence-afp/
#  trace.convergence-afp.txt, read 2026-10-17), except where another source is named:
#
#   step        measured there                      asked here      why
#   phmmer      16 cores, 15-37 s, 0.16 GB          16 cores 32 GB 2 h   48 queries, 2.4x
#   jackhmmer   16 cores, 24-57 s, 0.20 GB          16 cores 32 GB 6 h   3 rounds, long tail
#   mmseqs2     16 cores, 17 s, 2.4 GB              16 cores 32 GB 2 h
#   foldseek    16 cores, 13 s, 2.3 GB              16 cores 32 GB 2 h   AFDB Swiss-Prot db exists
#   prostt5     16 cores + 1 GPU, 42 s, 2.6 GB      16 cores 32 GB 2 h, 1 GPU
#   kmerseek    8 cores, 12-23 min, 45-52 GB        8 cores 64 GB 4 h    gbmr4 k20 and k22, scaled 1
#   hmmscan     Table 1: 4 cores, 0.39 core-h/1k    4 cores 8 GB 2 h     Pfam only, no Swiss-Prot
#   hhsearch    Table 1: 16 cores, 24.8 core-h/1k   16 cores 64 GB 48 h  Swiss-Prot is 14x one
#                                                                        proteome; see the note
#
# The hhsearch step is the one to watch: 24.8 core-hours per 1,000 queries against one 25,000
# protein proteome scales to roughly 0.7 core-hours per query against 575,748 Swiss-Prot
# entries, so 48 queries is about 34 core-hours, a little over 2 hours on 16 cores. The 48-hour
# ask is four times that, because no run has measured it.
#
#SBATCH --job-name=known-cases
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/known-cases-%j.out

set -uo pipefail

REPO=${KC_REPO:?set KC_REPO to the checkout of protein-function-annotation-leaderboard}
OUT=${KC_OUT:-$SCRATCH/known-cases}
CONV=$SCRATCH/2024-kmerseek-convergence-afp/nextflow-runs/convergence-afp
SPROT=${KC_SPROT:-$CONV/results/inputs/reference.fasta}       # Swiss-Prot 2026_03, 575,748 entries
AFDB=${KC_AFDB:-$CONV/results/db_cache/foldseek/afdb_swissprot}
QUERIES=$REPO/data/known_cases/queries.fasta
IMG=${KC_IMG:-$SCRATCH/apptainer_cache/kmerseek_0.4.0.sif}
mkdir -p "$OUT" logs

echo "queries: $(grep -c '^>' "$QUERIES")   target: $SPROT"

run() {  # run <name> <command...>, and record what it cost
  local name=$1; shift
  local t0=$SECONDS
  "$@" > "$OUT/$name.raw" 2> "$OUT/$name.err"
  local status=$?
  printf '%s\t%s\t%s\t%s\n' "$name" "$status" "$((SECONDS - t0))" "$SLURM_CPUS_PER_TASK" \
    >> "$OUT/resources.tsv"
  echo "$name: exit $status after $((SECONDS - t0)) s"
}

module load biology 2>/dev/null

# 1. phmmer and jackhmmer, HMMER 3.4, E-value cutoff 10, as Table 1 runs them.
run phmmer phmmer --domtblout "$OUT/phmmer.domtbl" -o /dev/null --noali -E 10.0 \
  --cpu "$SLURM_CPUS_PER_TASK" "$QUERIES" "$SPROT"
run jackhmmer jackhmmer -N 3 --domtblout "$OUT/jackhmmer.domtbl" -o /dev/null --noali -E 10.0 \
  --cpu "$SLURM_CPUS_PER_TASK" "$QUERIES" "$SPROT"

# 2. hmmscan against Pfam-A 38.2. It searches no relative, so Swiss-Prot does not enter.
run hmmscan hmmscan --domtblout "$OUT/hmmscan.domtbl" -o /dev/null --noali -E 10.0 --domE 10.0 \
  --cpu 4 "${KC_PFAM:?set KC_PFAM to Pfam-A.hmm}" "$QUERIES"

# 3. MMseqs2, one round and three rounds, -s 7, E-value cutoff 10.
mmseqs createdb "$QUERIES" "$OUT/qdb" && mmseqs createdb "$SPROT" "$OUT/tdb"
run mmseqs2_seqseq mmseqs search "$OUT/qdb" "$OUT/tdb" "$OUT/res1" "$OUT/tmp1" \
  --threads "$SLURM_CPUS_PER_TASK" -s 7 --max-seqs 1000 -e 10.0
mmseqs convertalis "$OUT/qdb" "$OUT/tdb" "$OUT/res1" "$OUT/mmseqs2_seqseq.tsv" \
  --format-output "query,target,qstart,qend,tstart,tend,bits,evalue"
run mmseqs2_iterative mmseqs search "$OUT/qdb" "$OUT/tdb" "$OUT/res3" "$OUT/tmp3" \
  --threads "$SLURM_CPUS_PER_TASK" -s 7 --num-iterations 3 --max-seqs 1000 -e 10.0
mmseqs convertalis "$OUT/qdb" "$OUT/tdb" "$OUT/res3" "$OUT/mmseqs2_iterative.tsv" \
  --format-output "query,target,qstart,qend,tstart,tend,bits,evalue"

# 4. Foldseek. Queries are the experimental structures named on the cards where the paper has
# one, and the AlphaFold Database model otherwise; data/known_cases/cases.json says which, and
# scripts/known_cases_structures.sh fetches them into $OUT/structures.
bash "$REPO/scripts/known_cases_structures.sh" "$REPO/data/known_cases/cases.json" "$OUT/structures"
foldseek createdb "$OUT/structures" "$OUT/qfs"
run foldseek foldseek search "$OUT/qfs" "$AFDB" "$OUT/fsres" "$OUT/fstmp" \
  --threads "$SLURM_CPUS_PER_TASK" -e 10.0 --max-seqs 1000
foldseek convertalis "$OUT/qfs" "$AFDB" "$OUT/fsres" "$OUT/foldseek.tsv" \
  --format-output "query,target,qstart,qend,tstart,tend,bits,evalue"

# 5. ProstT5 and Foldseek: the shape letters come from the sequence, so no structure is needed.
run prostt5_db foldseek createdb "$QUERIES" "$OUT/qprost" \
  --prostt5-model "${KC_PROSTT5:?set KC_PROSTT5 to the ProstT5 weights}" --gpu 1 --threads 8
run prostt5 foldseek search "$OUT/qprost" "$AFDB" "$OUT/prres" "$OUT/prtmp" \
  --threads "$SLURM_CPUS_PER_TASK" -e 10.0 --max-seqs 1000
foldseek convertalis "$OUT/qprost" "$AFDB" "$OUT/prres" "$OUT/prostt5.tsv" \
  --format-output "query,target,qstart,qend,tstart,tend,bits,evalue"

# 6. HH-suite: one profile per query protein, as Table 1 runs it (no UniRef30).
run hhblits ffindex_apply "$OUT/qhh_a3m.ffdata" "$OUT/qhh_a3m.ffindex" \
  -d "$OUT/hhres.ffdata" -i "$OUT/hhres.ffindex" \
  -- hhsearch -i stdin -d "${KC_HHDB:?set KC_HHDB to a Swiss-Prot HH-suite database}" \
     -blasttab /dev/stdout -cpu "$SLURM_CPUS_PER_TASK" -v 0

# 7. kmerseek at Table 1's best setting, gbmr4 k21, scaled 1.
run kmerseek apptainer exec --bind /scratch "$IMG" kmerseek search \
  --alphabet gbmr4 --ksize 21 --scaled 1 \
  --query "$QUERIES" --target "${KC_KMERSEEK_INDEX:?set KC_KMERSEEK_INDEX to a Swiss-Prot index}" \
  --threshold 0.0 --min-shared-kmers 2 --max-query-pvalue 0.05 --min-region-score 1.3

echo "done. Copy $OUT back and run: .venv/bin/python scripts/known_cases_results.py --run $OUT"
