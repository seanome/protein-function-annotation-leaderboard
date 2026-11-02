#!/bin/bash
# ESMFold2 (Biohub, MIT license; github.com/biohub/esm): install the `esm` package and fetch the
# weights, so esmfold2_fold.sbatch can run with no downloads.
#
#   sbatch esmfold2_install.sbatch
#
# `esm` 3.4.2 needs Python 3.12 or newer and PyTorch 2.11. PyTorch 2.11's pip builds need glibc
# 2.28, and Sherlock runs CentOS 7 (glibc 2.17), where pip offers torch up to 2.6 only (the first
# try, job 50705798, failed on that). So everything runs inside PyTorch's own container,
# docker://pytorch/pytorch:2.11.0-cuda12.8-cudnn9-runtime (Ubuntu, Python 3.12, torch 2.11.0),
# through Apptainer, and `esm` goes into a venv on top of the container's torch. The container
# file and the venv go in $GROUP_HOME (Sherlock keeps environments out of $SCRATCH, which is
# purged); the Apptainer, pip and Hugging Face caches go in $SCRATCH, since all can be fetched
# again. Weights, read 2026-11-01 from huggingface.co: biohub/ESMFold2 (26.7 GB) and
# biohub/ESMC-6B (25.4 GB), the language model it is built on. A rerun skips what is already there.
#SBATCH -p hns,normal
#SBATCH -A ayeletv
#SBATCH -J esmfold2-install
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH -t 2:00:00
#SBATCH -o %x-%j.out
set -euo pipefail
D=$GROUP_HOME/olgabot/envs/esmfold2
SIF=$D/pytorch-2.11.0-cuda12.8-cudnn9-runtime.sif
ENV=$D/venv-esm3.4.2
export APPTAINER_CACHEDIR=$SCRATCH/cache/apptainer APPTAINER_TMPDIR=$SCRATCH/cache/apptainer/tmp PIP_CACHE_DIR=$SCRATCH/cache/pip HF_HOME=$SCRATCH/cache/hf
mkdir -p "$D" "$APPTAINER_TMPDIR" "$PIP_CACHE_DIR" "$HF_HOME"
[ -s "$SIF" ] || apptainer pull "$SIF" docker://pytorch/pytorch:2.11.0-cuda12.8-cudnn9-runtime
run() { apptainer exec -B /home/groups,/scratch "$SIF" "$@"; }
# The image has no ensurepip (job 50726351 failed on it), so the venv is made without pip and uses
# the image's own pip, which installs into the venv it is run from.
[ -f "$ENV/pyvenv.cfg" ] || run python3 -m venv --without-pip --system-site-packages "$ENV"
run "$ENV/bin/python" -m pip install "esm==3.4.2"
run "$ENV/bin/python" - <<'PY'
import esm, torch, huggingface_hub
from esm.models.esmfold2 import ESMFold2InputBuilder, EsmFold2Model  # the names esmfold2_fold.py uses
print("esm", esm.__version__, "torch", torch.__version__, "CUDA", torch.version.cuda)
for repo in ("biohub/ESMFold2", "biohub/ESMC-6B"):
    print(repo, huggingface_hub.snapshot_download(repo_id=repo))
PY
