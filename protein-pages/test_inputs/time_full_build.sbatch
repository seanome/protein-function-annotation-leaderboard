#!/bin/bash
# Time build.py on every FUN protein that has search output: 7,216 proteins, with all of
# their phmmer, jackhmmer, MMseqs2 three-round and kmerseek (hp_lehninger2, k = 17) rows, and
# all of Swiss-Prot 2026_03 as swissprot.fasta. hmmscan output exists only for the 50 test
# proteins, so it is copied from the test set.
#
#   cd <folder holding build.py> && sbatch test_inputs/time_full_build.sbatch $SCRATCH/protein-pages-full
#
# 8 CPUs, the number a laptop has. Memory: build.py keeps only hits at E <= 1e-3.
#SBATCH -p hns,normal
#SBATCH -A ayeletv
#SBATCH -J protein-pages-time
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -t 2:00:00
set -euo pipefail
OUT=$1
HERE=$SLURM_SUBMIT_DIR
S=$SCRATCH
IN=$OUT/inputs
mkdir -p "$IN"
cat $S/botryllus-mhc/qpos/chunks/chunk_*.fasta > "$IN/proteins.fasta"
zcat $S/botryllus-mhc/qpos/out/mmseqs2.chunk_*.tsv.gz \
  | awk -F'\t' -v OFS='\t' '{print $1,$2,"","",$3,$4,$5,$6,$8,$7,""}' > "$IN/mmseqs2_3rounds.m8"
zcat $S/botryllus-mhc/qpos/out/phmmer.chunk_*.domtblout.gz > "$IN/phmmer.domtblout"
zcat $S/botryllus-mhc/qpos/out/jackhmmer.chunk_*.domtblout.gz > "$IN/jackhmmer.domtblout"
zcat $S/invertebrate-dark-set/results-0.4/dark-regions/botryllus_dark_regions.csv.gz \
  | python3 -c '
import csv, sys
r = csv.reader(sys.stdin); w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
h = next(r); w.writerow(h); im, ik = h.index("moltype"), h.index("ksize")
for f in r:
    if f[im] == "hp_lehninger2" and f[ik] == "17":
        w.writerow([x.replace("\t", " ") for x in f])' > "$IN/kmerseek.tsv"
zcat $S/botryllus-mhc-swissprot/db/swissprot/uniprot_sprot.fasta.gz > "$IN/swissprot.fasta"
cp $S/protein-pages-test/inputs/hmmscan.domtblout "$IN/"
cp "$HERE/test_inputs/run_settings.tsv" "$HERE/test_inputs/kmerseek_alphabets.rs" "$IN/"
ls -la "$IN"
ml biology py-biopython/1.79_py39
/usr/bin/time -v python3 "$HERE/build.py" "$IN" --out "$OUT/site" --jobs "$SLURM_CPUS_PER_TASK"
du -sh "$OUT/site"
