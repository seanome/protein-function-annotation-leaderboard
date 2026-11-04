#!/bin/bash
# ESMFold2 on one GPU: proteins of a FASTA, all of them or the IDs given (esmfold2_fold.py).
#
#   sbatch esmfold2_fold.sbatch <proteins.fasta> <out folder> [protein ID ...]
#
# Needs esmfold2_install.sbatch to have finished. The ESMC backbone runs in fp8, which needs an
# H100, hence the constraint. Loading reads ESMC's fp32 copy into CPU memory before casting it
# (esm 3.4.2, EsmFold2Model.from_pretrained), about 24 GB for 6 billion parameters, plus the
# folding trunk: 96 GB leaves room. The time is for a start on one protein (BHF, 252 aa): loading
# and one fold are not measured yet; the first run's log and esmfold2_summary.tsv will give the
# seconds per protein to size a run on the whole proteome.
#SBATCH -p hns,gpu
#SBATCH -A ayeletv
#SBATCH -J esmfold2
#SBATCH -G 1
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -c 8
#SBATCH --mem 96G
#SBATCH -t 1:00:00
#SBATCH -o %x-%j.out
set -euo pipefail
D=$GROUP_HOME/olgabot/envs/esmfold2  # the container and venv esmfold2_install.sbatch made
SIF=$D/pytorch-2.11.0-cuda12.8-cudnn9-runtime.sif
ENV=$D/venv-esm3.4.2
export HF_HOME=$SCRATCH/cache/hf HF_HUB_OFFLINE=1
run() { apptainer exec --nv -B /home/groups,/scratch "$SIF" "$@"; }  # --nv passes the GPU in
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
run "$ENV/bin/python" -c "import esm, torch; print('esm', esm.__version__, 'torch', torch.__version__, 'GPU', torch.cuda.is_available())"
run "$ENV/bin/python" "$SLURM_SUBMIT_DIR/esmfold2_fold.py" "$@"  # Slurm runs a copy of this file, so not $0
