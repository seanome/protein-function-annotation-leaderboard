#!/bin/bash
# Solvers unreadable, step 4 of 4: Slurm's record of every job in the run, then one table
# row per proteome (collect.py). Runs after the label array ends, whether or not every task
# succeeded (afterany), so failed proteomes are listed rather than lost.
# Submitted by submit_all.sh (RUN, REPO, PROTEOMES, the job ids). Not meant to be run alone.
#SBATCH --job-name=unread-collect
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
set -euo pipefail

: "${RUN:?}" "${REPO:?}" "${PROTEOMES:?}" "${SETUP_JOB:?}" "${FETCH_JOB:?}" "${LABEL_JOB:?}"
sacct -j "$SETUP_JOB,$FETCH_JOB,$LABEL_JOB" \
    --format=JobID%30,State,Partition,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,NodeList,Start,End -P \
    > "$RUN/sacct.txt"
/usr/bin/python3 -I "$REPO/scripts/nature_solved_it/unreadable/collect.py" --proteomes "$PROTEOMES" \
    --run "$RUN" --sacct "$RUN/sacct.txt" --fetch-job "$FETCH_JOB" --label-job "$LABEL_JOB"
