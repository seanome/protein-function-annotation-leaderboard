#!/bin/bash
# Copy a finished run's small outputs from Sherlock into the repository (read-only on Sherlock).
# Copies the result table, Slurm's record, the setup marker and versions, the setup and collect
# logs, the earlier-run folders, the download manifests, and per proteome: summary.json (the
# counts every share is computed from), input checksums, tool versions and BUSCO's short
# summaries. The MMseqs2 tables (2 GB in all) stay on Sherlock.
#
#   bash scripts/nature_solved_it/unreadable/copy_results.sh <run date> <setup job> <collect job> [extra log job ids]
set -euo pipefail
DATE=${1:?run date, e.g. 2026-10-20}
SETUP=${2:?setup job id}
COLLECT=${3:?collect job id}
shift 3
REPO=$(cd "$(dirname "$0")/../../.." && pwd)
DEST=$REPO/results/nature-solved-it/unreadable/run_$DATE
mkdir -p "$DEST"
logs="logs/unread-setup_$SETUP.log logs/unread-collect_$COLLECT.log"
for j in "$@"; do logs="$logs logs/unread-setup_$j.log"; done
ssh -o BatchMode=yes sherlock "cd \$SCRATCH/nature-solved-it/unreadable && tar czf - \
    unreadable_table.tsv sacct.txt setup_ok busco_version.txt mmseqs_version.txt $logs \
    earlier-run-* job_ids*.txt \$(ls submissions.log 2>/dev/null) inputs/*manifest.json \
    out/*/summary.json out/*/inputs_sha256.txt out/*/busco_version.txt out/*/mmseqs_version.txt \
    out/*/busco_*/short_summary*.txt out/*/busco_*/short_summary*.json" | tar xzf - -C "$DEST"
echo "copied to $DEST: $(find "$DEST" -type f | wc -l | tr -d ' ') files, $(du -sh "$DEST" | cut -f1)"
