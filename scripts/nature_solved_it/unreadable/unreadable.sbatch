#!/bin/bash
# Solvers unreadable, step 3 of 4: one array task per row of proteomes.tsv. Task N starts
# when fetch task N has finished (Slurm aftercorr). Each task runs, for one proteome:
#   1. the leaderboard's MMseqs2 search against Swiss-Prot (command copied verbatim from
#      scripts/why/sherlock_labels.sbatch, task 0),
#   2. BUSCO on the kingdom lineage, and on eukaryota_odb12 too (PLAN.md),
#   3. count_unreadable.py, which writes summary.json with the unlabelled share.
# Generalises scripts/why/sherlock_labels.sbatch from two fixed proteomes to a list.
# Submitted by submit_all.sh (array size, --time from predict_resources.py, RUN, REPO,
# PROTEOMES, logs). Not meant to be run alone.
#
# 16 CPUs and 32G are the asks that worked in job 47961379 (peak 4.24 GiB for MMseqs2 on
# 50,928 queries, 2.77 GiB for BUSCO); resources.md has the per-proteome prediction.
#SBATCH --job-name=unread-label
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
set -euo pipefail

: "${RUN:?}" "${REPO:?}" "${PROTEOMES:?}"
HERE=$REPO/scripts/nature_solved_it/unreadable
IMG_DIR=$SCRATCH/apptainer-cache
# The MMseqs2 image of scripts/why/sherlock_labels.sbatch (MMseqs2 18.8cc5c).
MMSEQS_IMG=$IMG_DIR/quay.io-biocontainers-mmseqs2@sha256-3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef.img
# This run's own BUSCO sandbox, built and checked by setup.sbatch (the shared one in
# $IMG_DIR lost its program loader on 2026-10-18; see setup.sbatch).
BUSCO_IMG=$RUN/images/busco-6.1.0--pyhdfd78af_2.sandbox
IN=$RUN/inputs
SPROT=$IN/uniprot_sprot.fasta.gz

# The row's fields, tab-separated: upid, taxid, kingdom lineage, all-eukaryote lineage, then
# the row as JSON for the summary.
IFS=$'\t' read -r UPID TAXID LINEAGE LINEAGE_ALL META < <(/usr/bin/python3 -I -c '
import csv, json, sys
r = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))[int(sys.argv[2])]
keep = ("group", "species", "upid", "proteome_taxid", "divergence_my", "busco_lineage", "busco_lineage_all")
print("\t".join([r["upid"], r["proteome_taxid"], r["busco_lineage"], r["busco_lineage_all"], json.dumps({k: r[k] for k in keep})]))
' "$PROTEOMES" "$SLURM_ARRAY_TASK_ID")
OUT=$RUN/out/$UPID
mkdir -p "$OUT"
cd "$OUT"
FASTA=$IN/${UPID}_${TAXID}.fasta.gz
sha256sum "$FASTA" "$SPROT" "$IN/${UPID}_uniprot.tsv.gz" > inputs_sha256.txt

# 1. MMseqs2, same command and settings as the leaderboard.
run() { apptainer exec --bind /scratch "$MMSEQS_IMG" "$@"; }
run mmseqs version > mmseqs_version.txt
command cp -f "$FASTA" queries.fasta.gz
tmp=$(mktemp -d "$OUT/mmseqs_tmp.XXXX")
    run mmseqs easy-search queries.fasta.gz "$SPROT" mmseqs_vs_sprot.tsv "$tmp" \
        -s 7.5 -e 1e-3 --max-seqs 1000 --threads "$SLURM_CPUS_PER_TASK" \
        --format-output query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen
gzip -f mmseqs_vs_sprot.tsv
command rm -rf "$tmp" queries.fasta.gz

# 2. BUSCO, same command as the leaderboard, once per lineage.
run() { apptainer exec --bind /scratch "$BUSCO_IMG" "$@"; }
run busco --version > busco_version.txt
gunzip -c "$FASTA" > proteins.fasta
busco_dirs=()
for LINEAGE_NOW in $(printf '%s\n' "$LINEAGE" "$LINEAGE_ALL" | awk '!seen[$0]++'); do
    fa=proteins tag=$LINEAGE_NOW
    # The leaderboard's names: "$fa.fasta", "busco_downloads_${tag}/lineages/$LINEAGE".
    ln -sfn "$RUN/busco_downloads" "busco_downloads_${tag}"
    LINEAGE=$LINEAGE_NOW
    run busco -i "$fa.fasta" -m proteins -l "busco_downloads_${tag}/lineages/$LINEAGE" --offline \
        -c "$SLURM_CPUS_PER_TASK" -o "busco_${tag}" -f --opt-out-run-stats
    busco_dirs+=(--busco "busco_${tag}")
done
command rm -f proteins.fasta

# 3. The shares.
/usr/bin/python3 -I "$HERE/count_unreadable.py" --fasta "$FASTA" --uniprot "$IN/${UPID}_uniprot.tsv.gz" \
    --mmseqs mmseqs_vs_sprot.tsv.gz "${busco_dirs[@]}" --meta "$META" --out summary.json
echo "done task $SLURM_ARRAY_TASK_ID $UPID"
