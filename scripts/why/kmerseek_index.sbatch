#!/bin/bash
# Build one kmerseek index of Swiss-Prot 2026_03 for one alphabet and k, with the
# Karlin-Altschul E-value fit stored in it. Submitted by kmerseek_submit.sh, which sets
# --mem and --time from results/why/kmerseek/settings.tsv and passes the row through
# --export (ALPHABET, KSIZE, PENALTY, XDROP, KA_QUERIES, WHY, KDIR).
# Flags follow the dark-set pipeline's kmerseekIndex (seanome/2024-kmerseek-analysis,
# nextflow-runs/invertebrate-dark-set/main.nf at ec91da53): scaled 1, no low-complexity mask.
#SBATCH --job-name=why-kmerseek-index
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --output=%x_%j.log
set -euo pipefail
: "${WHY:?}" "${KDIR:?}" "${ALPHABET:?}" "${KSIZE:?}" "${PENALTY:?}" "${XDROP:?}" "${KA_QUERIES:?}"
IMG=$SCRATCH/apptainer-cache/docker.io-olgabot-kmerseek-0.4.0-pr141-3c353d4.img
run() { apptainer exec --bind /scratch "$IMG" "$@"; }
stamp() { while IFS= read -r line; do printf '%(%Y-%m-%dT%H:%M:%S)T %s\n' -1 "$line"; done; }

REF=$KDIR/reference.fasta
IDX=$KDIR/indexes/sprot_2026_03.${ALPHABET}.k${KSIZE}.lcfalse.kmerseek.rocksdb
mkdir -p "$KDIR/indexes"

# One reference for every index: accession-only headers, as the pipeline writes them.
flock "$REF.lock" bash -c "
  [[ -s '$REF' ]] && exit 0
  zcat '${SPROT_GZ:-$WHY/inputs/uniprot_sprot.fasta.gz}' | awk '/^>/ { split(\$1, a, \"|\"); print \">\" a[2]; next } { print }' > '$REF.tmp'
  mv '$REF.tmp' '$REF'
"
if [[ -f "$IDX/BUILD_DONE" ]]; then
    echo "index already built: $IDX (delete it to rebuild)"; exit 0
fi
rm -rf "$IDX"   # a half-built index from a killed run
run kmerseek --version | stamp
echo start | stamp | tee "$IDX.log"
run kmerseek index \
    --alphabet "$ALPHABET" --ksize "$KSIZE" --scaled 1 \
    --input "$REF" --output "$IDX" \
    --extend-mismatch-penalty "$PENALTY" --extend-xdrop "$XDROP" \
    --ka-queries "$KA_QUERIES" --ka-survival-out "$IDX/ka_survival.C${PENALTY}.csv" \
    --kmer-stats-out "$IDX/spectrum.csv.gz" 2>&1 | stamp | tee -a "$IDX.log"
mv "$IDX.log" "$IDX/index.log"
# Searches open the index read-only; read-only files make a read-write open fail loudly.
find "$IDX" -type f -exec chmod a-w {} +
date -u +%FT%TZ > "$IDX/BUILD_DONE"
echo "built $IDX"
