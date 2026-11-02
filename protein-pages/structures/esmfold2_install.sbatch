#!/bin/bash
# ESMFold2 (Biohub, MIT license; github.com/biohub/esm): install the `esm` package and fetch the
# weights, so esmfold2_fold.sbatch can run with no downloads.
#
#   sbatch esmfold2_install.sbatch
#
# `esm` 3.4.2 needs Python 3.12 or newer and PyTorch 2.11 (pip pulls the CUDA build, about
# 5.6 GB). The environment goes in $GROUP_HOME (Sherlock keeps environments out of $SCRATCH,
# which is purged); pip's cache and the Hugging Face cache go in $SCRATCH, since both can be
# fetched again. Weights, read 2026-11-01 from huggingface.co: biohub/ESMFold2 (26.7 GB) and
# biohub/ESMC-6B (25.4 GB), the language model it is built on. A rerun skips what is already there.
#SBATCH -p hns,normal
#SBATCH -A ayeletv
#SBATCH -J esmfold2-install
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH -t 2:00:00
#SBATCH -o %x-%j.out
set -euo pipefail
ml purge
ml python/3.12.1
ENV=$GROUP_HOME/olgabot/envs/esmfold2-esm3.4.2
export PIP_CACHE_DIR=$SCRATCH/cache/pip HF_HOME=$SCRATCH/cache/hf
mkdir -p "$(dirname "$ENV")" "$PIP_CACHE_DIR" "$HF_HOME"
[ -x "$ENV/bin/python" ] || python3 -m venv "$ENV"
"$ENV/bin/pip" install --upgrade pip
"$ENV/bin/pip" install "esm==3.4.2"
"$ENV/bin/python" - <<'PY'
import esm, torch, huggingface_hub
from esm.models.esmfold2 import ESMFold2InputBuilder, EsmFold2Model  # the names esmfold2_fold.py uses
print("esm", esm.__version__, "torch", torch.__version__, "CUDA", torch.version.cuda)
for repo in ("biohub/ESMFold2", "biohub/ESMC-6B"):
    print(repo, huggingface_hub.snapshot_download(repo_id=repo))
PY
