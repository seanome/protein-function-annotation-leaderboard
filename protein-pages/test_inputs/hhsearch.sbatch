#!/bin/bash
# HH-suite on the test-set proteins, one array task per protein: hhblits builds the protein's
# alignment from UniRef30, then hhsearch compares that alignment with the PDB70 profiles.
#
#   cd <folder holding this file> && sbatch hhsearch.sbatch <test set folder>
#
# Needs hhsuite_install.sbatch to have finished (image and databases). Reads
# <test set folder>/inputs/proteins.fasta; task N takes its Nth sequence. Writes:
#   inputs/hhsearch/<ID>.hhr      hhsearch's output against PDB70, which build.py reads (track
#                                 "hhsearch": one .hhr per protein; the Query line is the FASTA ID)
#   hhblits/<ID>.a3m              the alignment hhblits built (kept so hhsearch can be rerun
#                                 without hhblits)
#   hhblits/<ID>.uniref30.hhr     hhblits's own hit list against UniRef30
# A task whose .hhr already exists does nothing, so a rerun only does the missing proteins. Each
# output is written in node-local space, copied in under a .part name and then renamed, so a
# killed task leaves no partial .a3m or .hhr.
#
# Options, as HHpred runs it: hhblits -n 3 -e 1e-3 (three rounds, a hit joins the alignment at
# E-value 1e-3 or below), then hhsearch with its defaults for what is listed: hits with Prob of
# at least 20 (-p 20), at most 500 in the hit list (-Z 500) and 500 alignments (-B 500). build.py
# draws hits with Prob at least 50, so nothing it draws is cut by these. -maxmem 20 lets the
# realignment step use up to 20 GB (default 3 GB) so long hits on the longest proteins are still
# realigned rather than skipped.
#
# Sizes: the 51 proteins are 54-2,160 aa (18,295 in all). hhblits holds the UniRef30 column-state
# file (UniRef30_2023_02_cs219.ffdata; hhsuite_install.sbatch prints its size) for the prefilter,
# plus dynamic-programming matrices that grow with query length x template length per thread,
# which -maxmem caps at 20 GB. Measured on the 51 (job 50673585, 2026-11-01): 8-29 min per task,
# peak memory up to 31.9 GB, so 48 GB per task (32 GB was nearly used up by task 30). On 8 CPUs, 3 rounds against UniRef30
# usually take minutes for a protein of a few hundred residues and up to about an hour for one of
# ~2,000 aa with a deep alignment; hhsearch against PDB70 adds a few minutes. These are estimates,
# not measured here (no HH-suite run on Sherlock yet): 4 h per task leaves room. If a task runs
# out of time, resubmit only that index with more time, e.g. --array=17 -t 12:00:00.
#SBATCH -p hns,normal
#SBATCH -A ayeletv
#SBATCH -J protein-pages-hhsearch
#SBATCH --array=1-51
#SBATCH -c 8
#SBATCH --mem 48G
#SBATCH -t 4:00:00
#SBATCH -o hhsearch-%A_%a.out
set -euo pipefail
OUT=${1:?usage: sbatch hhsearch.sbatch <test set folder>}
Q=$OUT/inputs/proteins.fasta
HHR_DIR=$OUT/inputs/hhsearch
A3M_DIR=$OUT/hhblits
IMG=$SCRATCH/apptainer-cache/quay.io-biocontainers-hhsuite@sha256-4bf9bb5229de18f522a94f4443c19fdcbb0f0cb0e6ea92f5390aa170bcb0a24f.img
DB=$SCRATCH/databases/hhsuite
UNIREF=$DB/UniRef30_2023_02/UniRef30_2023_02
PDB70=$DB/pdb70/pdb70
N=$SLURM_ARRAY_TASK_ID
CPU=$SLURM_CPUS_PER_TASK

n_seq=$(grep -c '^>' "$Q")
if [ "$N" -gt "$n_seq" ]; then
    echo "task $N: $Q has only $n_seq sequences" >&2
    exit 1
fi
for f in "$IMG" "$UNIREF"_cs219.ffdata "$PDB70"_hhm.ffdata; do
    [ -s "$f" ] || { echo "missing $f: run hhsuite_install.sbatch first" >&2; exit 1; }
done
ID=$(grep '^>' "$Q" | sed -n "${N}p" | sed 's/^>//; s/[[:space:]].*//')
mkdir -p "$HHR_DIR" "$A3M_DIR"
if [ -s "$HHR_DIR/$ID.hhr" ]; then
    echo "task $N: $HHR_DIR/$ID.hhr exists, skipping"
    exit 0
fi

TMP=$(mktemp -d "$L_SCRATCH/pp-hh.XXXXXX" 2>/dev/null || mktemp -d)
trap 'rm -rf "$TMP"' EXIT
awk -v id=">$ID" '/^>/{p=($1==id)} p' "$Q" | sed '/^>/!s/\*$//' > "$TMP/$ID.fasta"
echo "task $N: $ID, $(grep -v '^>' "$TMP/$ID.fasta" | tr -d '\n' | wc -c) aa, $(date)"
X="apptainer exec -B /scratch -B $TMP $IMG"

if [ ! -s "$A3M_DIR/$ID.a3m" ]; then
    $X hhblits -i "$TMP/$ID.fasta" -d "$UNIREF" -n 3 -e 1e-3 -cpu "$CPU" -maxmem 20 \
        -oa3m "$TMP/$ID.a3m" -o "$TMP/$ID.uniref30.hhr"
    mv "$TMP/$ID.uniref30.hhr" "$A3M_DIR/"
    mv "$TMP/$ID.a3m" "$A3M_DIR/$ID.a3m.part" && mv "$A3M_DIR/$ID.a3m.part" "$A3M_DIR/$ID.a3m"
    echo "hhblits done, $(grep -c '^>' "$A3M_DIR/$ID.a3m") sequences in the alignment, $(date)"
fi
$X hhsearch -i "$A3M_DIR/$ID.a3m" -d "$PDB70" -cpu "$CPU" -maxmem 20 -o "$TMP/$ID.hhr"
got=$(awk '/^Query /{print $2; exit}' "$TMP/$ID.hhr")
if [ "$got" != "$ID" ]; then
    echo "task $N: the .hhr names the query '$got', not '$ID'" >&2
    exit 1
fi
# copy across filesystems under a temporary name, then rename, so the .hhr appears whole or not at all
mv "$TMP/$ID.hhr" "$HHR_DIR/$ID.hhr.part" && mv "$HHR_DIR/$ID.hhr.part" "$HHR_DIR/$ID.hhr"
echo "hhsearch done: $(grep -cE '^ *[0-9]+ ' <(sed -n '/^ No Hit/,/^$/p' "$HHR_DIR/$ID.hhr")) hits listed, $(date)"
