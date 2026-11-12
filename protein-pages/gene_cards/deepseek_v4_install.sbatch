#!/bin/bash
#SBATCH --job-name=cards-v4-install
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
# vLLM 0.31.0 (released 2026-10-05) as vLLM's own container image, then the DeepSeek-V4-Flash
# weights (160 GB) into the Hugging Face cache on $SCRATCH. Sherlock's py-vllm module is
# 0.19.1, which cannot run DeepSeek V4, and pip cannot install a newer vLLM on Sherlock:
# CentOS 7's glibc 2.17 is too old for recent PyTorch wheels (job 50528214 found only
# torch up to 2.6). The image (docker://vllm/vllm-openai:v0.31.0, 9 GB) goes in $GROUP_HOME.
#
#   sbatch gene_cards/deepseek_v4_install.sbatch
set -euo pipefail
SIF=$GROUP_HOME/olgabot/containers/vllm-openai_v0.31.0.sif
export APPTAINER_CACHEDIR=$SCRATCH/.apptainer-cache APPTAINER_TMPDIR=${L_SCRATCH:-$SCRATCH}/apptainer-tmp
mkdir -p "$(dirname "$SIF")" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
[ -f "$SIF" ] || apptainer pull "$SIF" docker://vllm/vllm-openai:v0.31.0
apptainer exec "$SIF" python3 -c "import inspect, vllm, vllm.model_executor.models.registry as r; assert 'DeepseekV4ForCausalLM' in inspect.getsource(r); print('vllm', vllm.__version__, 'has DeepseekV4ForCausalLM')"
export HF_HUB_CACHE=$SCRATCH/.huggingface
apptainer exec --bind /scratch --env HF_HUB_CACHE="$HF_HUB_CACHE" "$SIF" hf download deepseek-ai/DeepSeek-V4-Flash
du -sh "$HF_HUB_CACHE/models--deepseek-ai--DeepSeek-V4-Flash"
