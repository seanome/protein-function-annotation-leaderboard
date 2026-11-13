#!/bin/bash
# Labels for step 4 and 6 of the "Why this project" numbers, run on Sherlock.
#   task 0: MMseqs2 search of human + oyster + the four proteins against Swiss-Prot 2026_03
#   task 1: BUSCO, human proteome, metazoa_odb12
#   task 2: BUSCO, oyster proteome, metazoa_odb12
# Inputs: $WHY/inputs/ (copied from the laptop; sha256 recorded in results/why/sources.json).
# Submit from the repo checkout on Sherlock:
#   sbatch --export=ALL,WHY=$SCRATCH/why-numbers/2026_03 scripts/why/sherlock_labels.sbatch
# Test on the four small proteins only (minutes), output in $WHY/out_test:
#   sbatch --export=ALL,WHY=$SCRATCH/why-numbers/2026_03,TEST=1 --time=00:30:00 scripts/why/sherlock_labels.sbatch
#SBATCH --job-name=why-labels
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.log
set -euo pipefail

: "${WHY:?set WHY to the run folder, e.g. \$SCRATCH/why-numbers/2026_03}"
IMG_DIR=$SCRATCH/apptainer-cache
# The MMseqs2 image the dark-set runs used, already converted on Sherlock.
MMSEQS_IMG=$IMG_DIR/quay.io-biocontainers-mmseqs2@sha256-3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef.img
BUSCO_URI=docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2
# A sandbox folder, not a .sif: building the .sif crashed in mksquashfs (exit 139) on
# Sherlock in the test job 47947485.
BUSCO_IMG=$IMG_DIR/quay.io-biocontainers-busco-6.1.0--pyhdfd78af_2.sandbox
LINEAGE=metazoa_odb12
export APPTAINER_CACHEDIR=$IMG_DIR/cache
IN=$WHY/inputs
OUT=$WHY/out
HUMAN=$IN/UP000005640_9606.fasta.gz
OYSTER=$IN/UP000005408_29159.fasta.gz
SPROT=$IN/uniprot_sprot.fasta.gz
if [[ ${TEST:-0} == 1 ]]; then
    OUT=$WHY/out_test HUMAN=$IN/four.fasta.gz OYSTER=$IN/four.fasta.gz SPROT=$IN/four.fasta.gz
fi
mkdir -p "$OUT"
cd "$OUT"
sha256sum "$IN"/*.fasta.gz > "inputs_sha256_task${SLURM_ARRAY_TASK_ID}.txt"

case "$SLURM_ARRAY_TASK_ID" in
0)
    run() { apptainer exec --bind /scratch "$MMSEQS_IMG" "$@"; }
    run mmseqs version > mmseqs_version.txt
    cat "$HUMAN" "$OYSTER" "$IN/four.fasta.gz" > queries.fasta.gz
    # On /scratch, which is bound into the container; node-local /lscratch is not.
    tmp=$(mktemp -d "$OUT/mmseqs_tmp.XXXX")
    # -s 7.5 is MMseqs2's most sensitive preset; E <= 1e-3 is its default; the 30% identity
    # cut is applied later on fident, so every hit is kept here, self-hits included.
    run mmseqs easy-search queries.fasta.gz "$SPROT" mmseqs_vs_sprot.tsv "$tmp" \
        -s 7.5 -e 1e-3 --max-seqs 1000 --threads "$SLURM_CPUS_PER_TASK" \
        --format-output query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen
    gzip -f mmseqs_vs_sprot.tsv
    rm -rf "$tmp"
    ;;
1|2)
    # Both BUSCO tasks need the image; the first to arrive pulls it, the other waits.
    flock "$BUSCO_IMG.lock" bash -c "[[ -d '$BUSCO_IMG/usr' ]] || apptainer build --sandbox '$BUSCO_IMG' '$BUSCO_URI'"
    run() { apptainer exec --bind /scratch "$BUSCO_IMG" "$@"; }
    if [[ $SLURM_ARRAY_TASK_ID == 1 ]]; then tag=human; fa_gz=$HUMAN; else tag=oyster; fa_gz=$OYSTER; fi
    fa=proteins_${tag}
    run busco --version > "busco_version_${tag}.txt"
    mkdir -p "busco_downloads_${tag}"
    run busco --download_path "busco_downloads_${tag}" --download "$LINEAGE" --opt-out-run-stats
    gunzip -c "$fa_gz" > "$fa.fasta"
    run busco -i "$fa.fasta" -m proteins -l "busco_downloads_${tag}/lineages/$LINEAGE" --offline \
        -c "$SLURM_CPUS_PER_TASK" -o "busco_${tag}" -f --opt-out-run-stats
    rm "$fa.fasta"
    ;;
esac
echo "done task $SLURM_ARRAY_TASK_ID"
