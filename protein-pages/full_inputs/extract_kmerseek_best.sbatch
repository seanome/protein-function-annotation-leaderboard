#!/bin/bash
#SBATCH --job-name=pages-kmerseek-best
#SBATCH --time=02:00:00
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --output=%x-%j.out
# Each protein's best kmerseek region, over the 35 settings whose E-value fit worked, from the whole
# proteome's match files. See extract_kmerseek_best.py.
#
#   sbatch extract_kmerseek_best.sbatch
set -euo pipefail
K=$SCRATCH/botryllus-kmerseek-swissprot/results
W=$SCRATCH/protein-pages-final/kmerseek-best
command rm -rf "$W"
mkdir -p "$W/parts"
awk -F'\t' 'NR > 1 && $5 == "fitted" {print $1 ".k" $2}' "$K/summary/kmerseek_botryllus_swissprot.pairs.tsv" > "$W/settings.txt"
echo "$(wc -l < "$W/settings.txt") fitted settings"
for s in $(cat "$W/settings.txt"); do ls "$K/matches/$s"/*.matches.csv.zst; done > "$W/files.txt"
export W S=$SLURM_SUBMIT_DIR
xargs -P 16 -I{} bash -c 'set -euo pipefail; zstdcat "{}" | python3 "$S/extract_kmerseek_best.py" chunk > "$W/parts/$(basename {} .matches.csv.zst).tsv" 2>> "$W/errors.txt"' < "$W/files.txt"
sort "$W/errors.txt" | uniq -c | sort -rn | head
python3 extract_kmerseek_best.py merge "$K/summary/kmerseek_botryllus_swissprot.per_protein.tsv.gz" "$W"/parts/*.tsv | gzip > "$W/kmerseek_best_regions.tsv.gz"
zcat "$W/kmerseek_best_regions.tsv.gz" | wc -l
