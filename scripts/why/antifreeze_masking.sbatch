#!/bin/bash
# Why does MMseqs2 find no hit for P24856 (antifreeze glycoprotein polyprotein, Notothenia
# neglecta, 790 residues, mostly the repeat AAT), not even its own Swiss-Prot entry?
# P24856 alone against Swiss-Prot 2026_03, eight ways:
#   MMseqs2 18.8cc5c (the image of job 47961379), -s 7.5 -e 1e-3 --max-seqs 1000, with
#     1 masking as job 47961379 ran it (defaults: --mask 1 --comp-bias-corr 1)
#     2 --mask 0            (no tantan low-complexity masking in the prefilter)
#     3 --comp-bias-corr 0  (no composition-bias correction)
#     4 --mask 0 --comp-bias-corr 0
#   phmmer and jackhmmer (HMMER 3.4), defaults, with and without --nobias (bias filter off).
# The flag names were read from `mmseqs easy-search -h` of the 18-8cc5c release binary and
# from `phmmer -h` / `jackhmmer -h` of HMMER 3.4; this job prints the same help from the
# images into help_*.txt so the PR can quote the images' own text.
#
#   sbatch --export=ALL,WHY=$SCRATCH/why-numbers/2026_03,HMMER_IMG=<HMMER 3.4 image the dark-set run uses> \
#       scripts/why/antifreeze_masking.sbatch
# Then, on the laptop, with the out_antifreeze/ folder copied to ~/data/why-numbers/2026_03/sherlock/:
#   python3 -I scripts/why/antifreeze.py --sherlock ~/data/why-numbers/2026_03/sherlock/out_antifreeze
#SBATCH --job-name=why-antifreeze
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.log
set -euo pipefail

: "${WHY:?set WHY to the run folder, e.g. \$SCRATCH/why-numbers/2026_03}"
: "${HMMER_IMG:?set HMMER_IMG to the HMMER 3.4 apptainer image}"
IMG_DIR=$SCRATCH/apptainer-cache
MMSEQS_IMG=$IMG_DIR/quay.io-biocontainers-mmseqs2@sha256-3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef.img
SPROT=$WHY/inputs/uniprot_sprot.fasta.gz
QUERY_ACC=P24856
OUT=$WHY/out_antifreeze
mkdir -p "$OUT"
cd "$OUT"
mm() { apptainer exec --bind /scratch "$MMSEQS_IMG" mmseqs "$@"; }
hm() { apptainer exec --bind /scratch "$HMMER_IMG" "$@"; }

sha256sum "$SPROT" > inputs_sha256.txt
echo "$MMSEQS_IMG" > images.txt
echo "$HMMER_IMG" >> images.txt
mm version > mmseqs_version.txt
mm easy-search -h > help_mmseqs_easy_search.txt 2>&1 || true
hm phmmer -h > help_phmmer.txt 2>&1 || true
hm jackhmmer -h > help_jackhmmer.txt 2>&1 || true
grep -E -- "--mask |--comp-bias-corr " help_mmseqs_easy_search.txt
grep -E -- "--nobias" help_phmmer.txt help_jackhmmer.txt

# The query, from the same Swiss-Prot file it is searched against.
gunzip -c "$SPROT" | awk -v a="$QUERY_ACC" '/^>/{p = ($0 ~ "^>sp\\|" a "\\|")} p' > query.fasta
grep -c '^>' query.fasta | grep -qx 1
gunzip -c "$SPROT" > sprot.fasta

cols=query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen
declare -A MM=(
    [mmseqs_as_job_47961379]=""
    [mmseqs_mask0]="--mask 0"
    [mmseqs_compbias0]="--comp-bias-corr 0"
    [mmseqs_mask0_compbias0]="--mask 0 --comp-bias-corr 0"
)
for name in mmseqs_as_job_47961379 mmseqs_mask0 mmseqs_compbias0 mmseqs_mask0_compbias0; do
    tmp=$(mktemp -d "$OUT/tmp.XXXX")
    # shellcheck disable=SC2086
    mm easy-search query.fasta sprot.fasta "$name.tsv" "$tmp" -s 7.5 -e 1e-3 --max-seqs 1000 --threads 1 \
        --format-output "$cols" ${MM[$name]} > "$name.log" 2>&1
    echo "mmseqs easy-search query.fasta sprot.fasta $name.tsv tmp -s 7.5 -e 1e-3 --max-seqs 1000 --threads 1 --format-output $cols ${MM[$name]}" > "$name.cmd"
    rm -rf "$tmp"
done
for tool in phmmer jackhmmer; do
    for bias in "" "--nobias"; do
        name=${tool}${bias:+_nobias}
        hm "$tool" --cpu 1 --noali $bias --tblout "$name.tblout" -o "$name.out" query.fasta sprot.fasta
        echo "$tool --cpu 1 --noali $bias --tblout $name.tblout -o $name.out query.fasta sprot.fasta" > "$name.cmd"
    done
done
rm sprot.fasta
sha256sum ./*.tsv ./*.tblout > outputs_sha256.txt
echo done
