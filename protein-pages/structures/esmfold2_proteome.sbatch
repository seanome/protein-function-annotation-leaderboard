#!/bin/bash
# ESMFold2 on the proteome: one array task per chunk of esmfold2_proteome_chunks.tsv, fused kernels.
#
#   cd structures && sbatch --array=1-<chunks>%4 esmfold2_proteome.sbatch
#
# Input: $SCRATCH/protein-pages-esmfold2/proteome/run.fasta, written by esmfold2_timing_sample.py
# fasta. Output: chunk_NNN/ beside it, holding the .pdb and .cif files and esmfold2_summary.tsv. A
# rerun of a task folds only the proteins of that chunk not yet done, so a time-out loses at most
# the protein in progress.
#
# Sizes from the fused-kernel timing test (job 50911817_1): a chunk is at most about 3 hours of
# folding, loading the model takes about 5 minutes, and the 5-hour limit leaves room for proteins
# slower than predicted. The most CPU memory used was 41 GB of the 96 GB asked. The 64 test proteins
# took 48 minutes, so a longer job could still grow past that; check the first task with sacct.
#SBATCH -p hns,gpu
#SBATCH -A ayeletv
#SBATCH -J esmfold2-proteome
#SBATCH -G 1
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -c 8
#SBATCH --mem 96G
#SBATCH -t 5:00:00
#SBATCH -o %x-%A_%a.out
set -euo pipefail
D=$GROUP_HOME/olgabot/envs/esmfold2  # the container and venv esmfold2_install.sbatch made
SIF=$D/pytorch-2.11.0-cuda12.8-cudnn9-runtime.sif
ENV=$D/venv-esm3.4.2
export HF_HOME=$SCRATCH/cache/hf HF_HUB_OFFLINE=1
run() { apptainer exec --nv -B /home/groups,/scratch "$SIF" "$@"; }  # --nv passes the GPU in
RUN=$SCRATCH/protein-pages-esmfold2/proteome
IDS=$(awk -F'\t' -v c="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $4 == c {print $1}' "$SLURM_SUBMIT_DIR/esmfold2_proteome_chunks.tsv")
if [ -z "$IDS" ]; then  # esmfold2_fold.py with no IDs would fold the whole FASTA
  echo "no proteins in chunk $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
run "$ENV/bin/python" "$SLURM_SUBMIT_DIR/esmfold2_fold.py" --kernels fused "$RUN/run.fasta" "$RUN/$(printf 'chunk_%03d' "$SLURM_ARRAY_TASK_ID")" $IDS
