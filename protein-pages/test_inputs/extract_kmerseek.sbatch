#!/bin/bash
#SBATCH --job-name=pages-kmerseek
#SBATCH --time=04:00:00
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out
# Every kmerseek region with E-value at most 10 for the 51 test proteins, from the six settings
# the invertebrate dark-set run scored per protein (results-0.4/protein-min-evalue/botryllus).
# Replaces the first kmerseek.tsv, which came from dark-regions/ and so only held proteins the
# dark-set run found no other hit for, at one setting.
#
#   sbatch test_inputs/extract_kmerseek.sbatch
set -euo pipefail
R=$SCRATCH/invertebrate-dark-set/results-0.4
IN=$SCRATCH/protein-pages-final/inputs
W=$SCRATCH/protein-pages-final/kmerseek-extract
mkdir -p "$W/parts"
grep '^>' "$IN/proteins.fasta" | sed 's/^>//; s/ .*//' | sort -u > "$W/ids.txt"
SETTINGS="hp_lehninger2.k17.s1.lcfalse.extend-c1.51 hp_lehninger2.k25.s1.lcfalse.extend-c1.51
hp_lehninger_c_nonpolar2.k19.s1.lcfalse.extend-c1.52 hp_pbotc_1st_ed2.k19.s1.lcfalse.extend-c1.59
hp_pbotc_1st_ed2.k21.s1.lcfalse.extend-c1.59 hp_pbotc_1st_ed2.k23.s1.lcfalse.extend-c1.59"
# The chunks that hold the 51, read from the per-protein table (column 1 source_file, 2 query_name)
zcat "$R"/protein-min-evalue/botryllus/*.hp_lehninger2.k17.s1.lcfalse.extend-c1.51.csv.gz |
  awk -F, 'NR==FNR {want[$1]=1; next} ($2 in want) {split($1, a, "."); print a[1]}' "$W/ids.txt" - | sort -u > "$W/chunks.txt"
echo "$(wc -l < "$W/ids.txt") proteins in $(wc -l < "$W/chunks.txt") chunks"
for s in $SETTINGS; do while read -r c; do echo "$c.$s"; done < "$W/chunks.txt"; done > "$W/files.txt"
# query_name is the first column and has no comma, so awk can pick the rows; python parses them
export R W
xargs -P 16 -I{} bash -c 'set -euo pipefail; f=$R/botryllus/kmerseek/{}.regions.csv.zst
  zstdcat "$f" | awk -F, -v ids="$W/ids.txt" "BEGIN {while ((getline l < ids) > 0) want[l]=1} NR==1 || (\$1 in want)" > "$W/parts/{}.csv"' \
  < "$W/files.txt"
python3 - "$W" "$IN/kmerseek.tsv" <<'EOF'
import csv, glob, sys
w, out = sys.argv[1], sys.argv[2]
n, head = 0, None
with open(out, "w") as fh:
    wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
    for p in sorted(glob.glob(w + "/parts/*.csv")):
        rd = csv.reader(open(p))
        h = next(rd, None)
        if h is None:
            continue
        if head is None:
            head = h
            wr.writerow(head)
        assert h == head, p
        e = head.index("region_evalue")
        for r in rd:
            if float(r[e]) <= 10:
                wr.writerow(r)
                n += 1
print(n, "regions at E-value 10 or below written to", out)
EOF
