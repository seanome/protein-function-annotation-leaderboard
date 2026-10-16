#!/bin/bash
# Search one chunk of query proteins (array index = chunk number) against one index.
# Submitted by kmerseek_submit.sh after the index job, with --mem/--time from settings.tsv.
# Flags follow the dark-set pipeline's kmerseekSearch (main.nf at ec91da53): extension at
# the alphabet's own penalty, X-drop 4C, chains up to 30 residues apart and 10 off the
# diagonal, --ka-queries 0 so the E-value uses the fit stored in the index (a refused fit is
# recorded as "nofit", as in the pipeline), regions kept up to E = 10_000.
#SBATCH --job-name=why-kmerseek-search
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --output=%x_%A_%a.log
set -euo pipefail
: "${WHY:?}" "${KDIR:?}" "${ALPHABET:?}" "${KSIZE:?}" "${PENALTY:?}" "${XDROP:?}" "${QDIR:?}"
IMG=$SCRATCH/apptainer-cache/docker.io-olgabot-kmerseek-0.4.0-pr141-3c353d4.img
run() { apptainer exec --bind /scratch "$IMG" "$@"; }

CHUNK=$(printf 'chunk_%04d' "$SLURM_ARRAY_TASK_ID")
IDX=$KDIR/indexes/sprot_2026_03.${ALPHABET}.k${KSIZE}.lcfalse.kmerseek.rocksdb
OUTD=$KDIR/out
SLUG=$CHUNK.${ALPHABET}.k${KSIZE}.s1.lcfalse.extend-c${PENALTY}
mkdir -p "$OUTD"
[[ -f "$IDX/BUILD_DONE" ]] || { echo "no finished index at $IDX"; exit 1; }
if [[ -f "$OUTD/$SLUG.DONE" ]]; then
    echo "already searched: $SLUG (delete $OUTD/$SLUG.* to redo)"; exit 0
fi
cd "$OUTD"
set +e
run kmerseek search \
    --alphabet "$ALPHABET" --ksize "$KSIZE" \
    --query "$QDIR/$CHUNK.fasta" --target "$IDX" \
    --extend-mismatch-penalty "$PENALTY" --extend-xdrop "$XDROP" --ka-queries 0 \
    --chain-max-gap 30 --chain-max-shift 10 \
    --threshold 0.0 --min-shared-kmers 2 --max-query-pvalue 0.05 --min-region-score 1.3 \
    2> "$SLUG.log" \
  | awk -F, -v max=10000 'NR==1 { for (i=1;i<=NF;i++) if ($i=="region_evalue") c=i;
        if (!c) { print "region_evalue column missing" > "/dev/stderr"; exit 3 } print; next }
        ($c+0) <= max' \
  | run zstd -q -T2 -f -o "$SLUG.regions.csv.zst"
status=("${PIPESTATUS[@]}")
set -e
if [[ ${status[0]} -ne 0 ]] && grep -q "no Karlin-Altschul fit" "$SLUG.log"; then
    printf 'nofit\t%s\n' "$(grep -m1 'no Karlin-Altschul fit' "$SLUG.log" | tr '\t' ' ')" > "$SLUG.nofit"
elif [[ ${status[0]} -ne 0 || ${status[1]} -ne 0 || ${status[2]} -ne 0 ]]; then
    echo "search failed: exit codes ${status[*]}; see $OUTD/$SLUG.log"; exit 1
fi
date -u +%FT%TZ > "$SLUG.DONE"
echo "searched $SLUG"
