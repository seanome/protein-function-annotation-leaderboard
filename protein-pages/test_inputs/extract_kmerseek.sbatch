#!/bin/bash
#SBATCH --job-name=pages-kmerseek
#SBATCH --time=02:00:00
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=%x-%j.out
# kmerseek.tsv for the 51 test proteins: every region with E-value at most 10 from the kmerseek
# run of the whole Botryllus proteome against all of Swiss-Prot 2026_03
# ($SCRATCH/botryllus-kmerseek-swissprot, 38 alphabet and k-mer size settings). Settings whose
# Karlin-Altschul fit was refused report a tail probability, not an E-value, and are left out
# (results/summary/kmerseek_botryllus_swissprot.pairs.tsv, column fit).
#
#   sbatch test_inputs/extract_kmerseek.sbatch
set -euo pipefail
K=$SCRATCH/botryllus-kmerseek-swissprot/results
IN=$SCRATCH/protein-pages-final/inputs
W=$SCRATCH/protein-pages-final/kmerseek-extract
command rm -rf "$W"
mkdir -p "$W/parts"
# The run renamed headers with a comma (query/query_renaming.tsv); the 51 have none, checked below
grep '^>' "$IN/proteins.fasta" | sed 's/^>//; s/ .*//' | sort -u > "$W/ids.txt"
if grep -qxFf "$W/ids.txt" <(cut -f1 "$K/query/query_renaming.tsv"); then echo "a test protein was renamed" >&2; exit 1; fi
awk -F'\t' 'NR > 1 && $5 == "fitted" {print $1 ".k" $2}' "$K/summary/kmerseek_botryllus_swissprot.pairs.tsv" > "$W/settings.txt"
echo "$(wc -l < "$W/ids.txt") proteins, $(wc -l < "$W/settings.txt") fitted settings"
for s in $(cat "$W/settings.txt"); do ls "$K/matches/$s"/*.matches.csv.zst; done > "$W/files.txt"
export W
# query_name is the first column and has no comma, so awk can pick the rows; python parses them
xargs -P 16 -I{} bash -c 'set -euo pipefail; zstdcat "{}" |
  awk -F, -v ids="$W/ids.txt" "BEGIN {while ((getline l < ids) > 0) want[l]=1} NR==1 || (\$1 in want)" > "$W/parts/$(basename {} .matches.csv.zst).csv"' \
  < "$W/files.txt"
python3 - "$W" "$IN/kmerseek.tsv" <<'EOF'
import csv, glob, sys
w, out = sys.argv[1], sys.argv[2]
KEEP = ["query_name", "target_name", "ksize", "scaled", "moltype", "region_start", "region_end", "region_subseq",
        "target_start", "target_end", "target_subseq", "moltype_seq", "region_length", "region_n_mismatches",
        "region_ka_bits", "region_evalue"]
n = 0
with open(out, "w") as fh:
    wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
    wr.writerow(KEEP)
    for p in sorted(glob.glob(w + "/parts/*.csv")):
        rd = csv.DictReader(open(p))
        for r in rd:
            if r["region_evalue"] and float(r["region_evalue"]) <= 10:
                wr.writerow([r[k] for k in KEEP])
                n += 1
print(n, "regions at E-value 10 or below written to", out)
EOF
