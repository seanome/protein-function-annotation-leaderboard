#!/bin/bash
# Rerun phmmer, jackhmmer and MMseqs2 (three rounds) on the test-set proteins with the tools'
# own alignments kept, so the protein pages can show each tool's alignment instead of one
# recomputed with Biopython.
#
#   cd <folder holding this file> && sbatch rerun_with_alignments.sbatch <test set folder>
#
# The original runs (2024-kmerseek-analysis notebook 281, $SCRATCH/botryllus-mhc/
# 281_tool_search_query_positions.sbatch) wrote only the tables: HMMER with -o /dev/null
# --noali, MMseqs2 without alignment columns. This uses the same container images, database
# and settings (E <= 10 reporting, jackhmmer -N 3, mmseqs2 -s 7 --num-iterations 3), and
# writes into <test set folder>/inputs:
#   phmmer.domtblout, phmmer.out         per-domain table and HMMER's main output (--notextw)
#   jackhmmer.domtblout, jackhmmer.out   the same for jackhmmer, all three rounds
#   mmseqs2_3rounds.m8                   query,target,fident,alnlen,qstart,qend,tstart,tend,
#                                        evalue,bits,theader,qaln,taln
#   hmmscan.out                          hmmscan's main output on the test proteins only (the
#                                        Swiss-Prot targets need only the table), same options
#                                        as make_test_inputs.sbatch
# The old tables are kept as *.from_notebook_281 so the two runs can be compared.
#
# Sizes: notebook 281's chunks of ~900 proteins took up to 1.2 h (phmmer), 3.8 h (jackhmmer)
# and 3 min (MMseqs2) on 16 CPUs; 51 proteins is about 1/18 of that, so 1.5 h with headroom.
#SBATCH -p hns,normal
#SBATCH -A ayeletv
#SBATCH -J protein-pages-aln
#SBATCH -c 16
#SBATCH --mem 32G
#SBATCH -t 1:30:00
set -euo pipefail
OUT=$1
IN=$OUT/inputs
REF=/scratch/users/olgabot/invertebrate-dark-set/reference/minus_Ascidiacea/reference.fasta
IMGS=$SCRATCH/apptainer-cache
HMMER=$IMGS/quay.io-biocontainers-hmmer@sha256-7a2b317b8d2fd3650b4924a8482cddeb940d4a0746c6a1501ff03ac1b7439e0c.img
MMSEQS=$IMGS/quay.io-biocontainers-mmseqs2@sha256-3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef.img
Q=$IN/proteins.fasta
TMP=$(mktemp -d "$L_SCRATCH/pp-aln.XXXXXX" 2>/dev/null || mktemp -d)
trap 'rm -rf "$TMP"' EXIT
for t in phmmer.domtblout jackhmmer.domtblout mmseqs2_3rounds.m8; do
    [ -f "$IN/$t.from_notebook_281" ] || mv "$IN/$t" "$IN/$t.from_notebook_281"
done

apptainer exec -B /scratch -B "$TMP" "$HMMER" phmmer --domtblout "$IN/phmmer.domtblout" --tblout /dev/null \
    -o "$IN/phmmer.out" --notextw -E 10 --cpu "$SLURM_CPUS_PER_TASK" "$Q" "$REF"
apptainer exec -B /scratch -B "$TMP" "$HMMER" jackhmmer -N 3 --domtblout "$IN/jackhmmer.domtblout" --tblout /dev/null \
    -o "$IN/jackhmmer.out" --notextw -E 10 --cpu "$SLURM_CPUS_PER_TASK" "$Q" "$REF"

apptainer exec -B /scratch "$HMMER" hmmscan --domtblout /dev/null -o "$IN/hmmscan.out" --notextw \
    -E 10.0 --domE 10.0 --cpu "$SLURM_CPUS_PER_TASK" "$SCRATCH/data/pfam/Pfam-A.hmm" "$Q"
cd "$TMP"
apptainer exec -B /scratch -B "$TMP" "$MMSEQS" mmseqs createdb "$Q" qdb
apptainer exec -B /scratch -B "$TMP" "$MMSEQS" mmseqs createdb "$REF" tdb
apptainer exec -B /scratch -B "$TMP" "$MMSEQS" mmseqs search qdb tdb res tmp -s 7 --num-iterations 3 -e 10 \
    --threads "$SLURM_CPUS_PER_TASK"
apptainer exec -B /scratch -B "$TMP" "$MMSEQS" mmseqs convertalis qdb tdb res "$IN/mmseqs2_3rounds.m8" \
    --format-output 'query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,theader,qaln,taln'

# Same hits as notebook 281's run? Compare (query, target, positions, E-value) row sets.
python3 - "$IN" <<'EOF'
import sys
d = sys.argv[1]
def dom(p):
    return {tuple(l.split()[i] for i in (0, 3, 9, 12, 15, 16, 17, 18)) for l in open(p) if not l.startswith("#")}
def m8_new(p):
    return {tuple(f[i] for i in (0, 1, 4, 5, 6, 7)) + ("%.3g" % float(f[8]),) for f in (l.rstrip("\n").split("\t") for l in open(p))}
def m8_old(p):
    return {tuple(f[i] for i in (0, 1, 4, 5, 6, 7)) + ("%.3g" % float(f[8]),) for f in (l.rstrip("\n").split("\t") for l in open(p))}
for name, a, b in (("phmmer", dom(d + "/phmmer.domtblout.from_notebook_281"), dom(d + "/phmmer.domtblout")),
                   ("jackhmmer", dom(d + "/jackhmmer.domtblout.from_notebook_281"), dom(d + "/jackhmmer.domtblout")),
                   ("mmseqs2_3rounds", m8_old(d + "/mmseqs2_3rounds.m8.from_notebook_281"), m8_new(d + "/mmseqs2_3rounds.m8"))):
    print("%s: notebook 281 %d rows, rerun %d, in both %d, only 281 %d, only rerun %d"
          % (name, len(a), len(b), len(a & b), len(a - b), len(b - a)), flush=True)
EOF
ls -la "$IN"
