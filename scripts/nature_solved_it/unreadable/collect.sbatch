#!/bin/bash
# Solvers unreadable, step 4 of 4: Slurm's record of every job in the run, then one table
# row per proteome (collect.py). Runs only if setup succeeded (afterok) and after the label
# array ends, whether or not every task succeeded (afterany), so failed proteomes are listed
# rather than lost. collect.py also refuses unless $RUN/setup_ok names this run's setup job.
# Submitted by submit_all.sh (RUN, REPO, PROTEOMES, the job ids). Not meant to be run alone.
#SBATCH --job-name=unread-collect
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
set -euo pipefail

: "${RUN:?}" "${REPO:?}" "${PROTEOMES:?}" "${SETUP_JOB:?}" "${FETCH_JOB:?}" "${LABEL_JOB:?}"
# No setup, no result: stop before writing anything.
grep -qx "setup_job $SETUP_JOB" "$RUN/setup_ok" 2>/dev/null || {
    echo "refusing to collect: $RUN/setup_ok does not name setup job $SETUP_JOB" >&2
    exit 2
}
sacct -j "$SETUP_JOB,$FETCH_JOB,$LABEL_JOB" \
    --format=JobID%30,State,Partition,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,NodeList,Start,End -P \
    > "$RUN/sacct.txt"
/usr/bin/python3 -I "$REPO/scripts/nature_solved_it/unreadable/collect.py" --proteomes "$PROTEOMES" \
    --run "$RUN" --sacct "$RUN/sacct.txt" --setup-job "$SETUP_JOB" --fetch-job "$FETCH_JOB" --label-job "$LABEL_JOB"
