#!/bin/bash
# Solvers unreadable: submit the whole run on Sherlock as one chain, then walk away.
#
#   setup (1 job)  --afterok-->  fetch (array, one task per proteome)
#                  --aftercorr--> label (array; task N starts when fetch task N succeeded)
#   collect (1 job): afterok setup AND afterany label, so it runs if some label tasks failed
#                    but can never run after a failed setup (--kill-on-invalid-dep=yes)
#
# Run from a Sherlock login node, in the repository clone:
#   bash scripts/nature_solved_it/unreadable/submit_all.sh [RUN_DIR]
# RUN_DIR defaults to $SCRATCH/nature-solved-it/unreadable. The run list is
# results/nature-solved-it/unreadable/proteomes.tsv of this checkout (draw_controls.py).
# Job ids go to $RUN/job_ids.txt; logs to $RUN/logs/.
# Outputs of an earlier submission (unreadable_table.tsv, sacct.txt, setup_ok, version files,
# job_ids.txt, out/) are moved, not deleted, to $RUN/earlier-run-<UTC time>/ first, so they
# cannot be read as this run's result. Downloaded inputs are kept: each is checked by sha256.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../../.." && pwd)
HERE=$REPO/scripts/nature_solved_it/unreadable
RES=$REPO/results/nature-solved-it/unreadable
RUN=${1:-$SCRATCH/nature-solved-it/unreadable}
PROTEOMES=$RES/proteomes.tsv
MMSEQS_IMG=$SCRATCH/apptainer-cache/quay.io-biocontainers-mmseqs2@sha256-3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef.img

test -s "$PROTEOMES" || { echo "no run list at $PROTEOMES; run draw_controls.py first" >&2; exit 1; }
test -s "$MMSEQS_IMG" || { echo "MMseqs2 image missing: $MMSEQS_IMG" >&2; exit 1; }
N=$(( $(wc -l < "$PROTEOMES") - 1 ))
test "$N" -gt 0
LAST=$(( N - 1 ))
LABEL_TIME=$(/usr/bin/python3 -I "$HERE/predict_resources.py" --proteomes "$PROTEOMES" --sacct "$RES/sacct_47961379.txt" --time-ask label)
FETCH_TIME=$(/usr/bin/python3 -I "$HERE/predict_resources.py" --proteomes "$PROTEOMES" --sacct "$RES/sacct_47961379.txt" --time-ask fetch)
mkdir -p "$RUN/logs"
stale=()
for f in unreadable_table.tsv sacct.txt setup_ok busco_version.txt mmseqs_version.txt job_ids.txt out; do
    if [[ -e $RUN/$f ]]; then stale+=("$f"); fi
done
if (( ${#stale[@]} )); then
    OLD=$RUN/earlier-run-$(date -u +%Y%m%dT%H%M%SZ)
    mkdir "$OLD"
    for f in "${stale[@]}"; do mv "$RUN/$f" "$OLD/"; done
    echo "moved earlier outputs to $OLD: ${stale[*]}"
fi
# Sherlock's git is 1.8.3.1, which has no "git -C".
COMMIT=$(cd "$REPO" && git rev-parse HEAD 2>/dev/null || echo "not a git checkout")
ENV="ALL,RUN=$RUN,REPO=$REPO,PROTEOMES=$PROTEOMES"

SETUP=$(sbatch --parsable --chdir "$RUN" --export="$ENV" --output="$RUN/logs/%x_%j.log" "$HERE/setup.sbatch")
FETCH=$(sbatch --parsable --chdir "$RUN" --export="$ENV" --output="$RUN/logs/%x_%A_%a.log" \
    --array="0-$LAST%8" --time="$FETCH_TIME" \
    --dependency="afterok:$SETUP" --kill-on-invalid-dep=yes "$HERE/fetch.sbatch")
LABEL=$(sbatch --parsable --chdir "$RUN" --export="$ENV" --output="$RUN/logs/%x_%A_%a.log" \
    --array="0-$LAST%20" --time="$LABEL_TIME" \
    --dependency="aftercorr:$FETCH" --kill-on-invalid-dep=yes "$HERE/unreadable.sbatch")
COLLECT=$(sbatch --parsable --chdir "$RUN" --output="$RUN/logs/%x_%j.log" \
    --export="$ENV,SETUP_JOB=$SETUP,FETCH_JOB=$FETCH,LABEL_JOB=$LABEL" \
    --dependency="afterok:$SETUP,afterany:$LABEL" --kill-on-invalid-dep=yes "$HERE/collect.sbatch")

{
    echo "submitted_utc $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "commit $COMMIT"
    echo "proteomes $N"
    echo "setup $SETUP"
    echo "fetch $FETCH array 0-$LAST%8 time $FETCH_TIME"
    echo "label $LABEL array 0-$LAST%20 time $LABEL_TIME"
    echo "collect $COLLECT"
} | tee "$RUN/job_ids.txt"
echo "Result table when collect $COLLECT ends: $RUN/unreadable_table.tsv"
