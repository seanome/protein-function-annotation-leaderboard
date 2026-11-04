#!/bin/bash
# ESMFold2 timing test: the 64 proteins of esmfold2_timing_ids.tsv, folded twice, as two array tasks
# on two H100s at once: task 0 with the package's default code, task 1 with its fused Triton
# kernels (--kernels fused). Same seed and settings, so the two should give the same structures;
# the summaries say whether they do, and how many seconds and GB each protein took.
#
#   cd structures && sbatch esmfold2_timing.sbatch $SCRATCH/protein-pages-esmfold2/timing/sample.fasta
#
# Output: $SCRATCH/protein-pages-esmfold2/timing/{default,fused}/esmfold2_summary.tsv, and the
# .pdb/.cif files beside it. A rerun folds only the proteins not yet done.
#
# Time: BHF (252 aa) took 20 s after about 6 minutes of loading. The sample has more long proteins
# than the proteome, on purpose: if time grows with length to the power 1.5, the first 61 take
# about 95 minutes and the three long ones (2,501, 3,501 and 4,988 aa) about an hour; to the power
# 2, about 3 and 4 hours. The long ones run last, so a time-out loses only them. Memory as in
# esmfold2_fold.sbatch (the model is cast on the CPU while loading).
#SBATCH -p hns,gpu
#SBATCH -A ayeletv
#SBATCH -J esmfold2-timing
#SBATCH -G 1
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -c 8
#SBATCH --mem 96G
#SBATCH -t 6:00:00
#SBATCH --array 0-1
#SBATCH -o %x-%A_%a.out
set -euo pipefail
D=$GROUP_HOME/olgabot/envs/esmfold2  # the container and venv esmfold2_install.sbatch made
SIF=$D/pytorch-2.11.0-cuda12.8-cudnn9-runtime.sif
ENV=$D/venv-esm3.4.2
export HF_HOME=$SCRATCH/cache/hf HF_HUB_OFFLINE=1
run() { apptainer exec --nv -B /home/groups,/scratch "$SIF" "$@"; }  # --nv passes the GPU in
FASTA=$1
OUT=$SCRATCH/protein-pages-esmfold2/timing
if [ "$SLURM_ARRAY_TASK_ID" = 1 ]; then KERNELS=(--kernels fused); OUT=$OUT/fused; else KERNELS=(); OUT=$OUT/default; fi
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
IDS=$(cut -f1 "$SLURM_SUBMIT_DIR/esmfold2_timing_ids.tsv" | tail -n +2)  # in the file's order: BHF first, long ones last
run "$ENV/bin/python" "$SLURM_SUBMIT_DIR/esmfold2_fold.py" ${KERNELS[@]+"${KERNELS[@]}"} "$FASTA" "$OUT" $IDS
