#!/bin/bash
#SBATCH --job-name=cards-v4-install
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
# vLLM 0.31.0 (released 2026-10-05) in a Python environment in $GROUP_HOME, then the
# DeepSeek-V4-Flash weights (160 GB) into the Hugging Face cache on $SCRATCH. Sherlock's
# py-vllm module is 0.19.1, which cannot run DeepSeek V4 (DeepseekV4ForCausalLM is not in its
# model list); vLLM 0.31.0 has the model, its chat format and its reasoning parser.
#
#   sbatch gene_cards/deepseek_v4_install.sbatch
set -euo pipefail
ENV=$GROUP_HOME/olgabot/envs/vllm-0.31.0
ml python/3.12.1
if [ ! -x "$ENV/bin/vllm" ]; then
  mkdir -p "$(dirname "$ENV")"
  python3 -m venv "$ENV"
  "$ENV/bin/pip" install --quiet --upgrade pip
  "$ENV/bin/pip" install --quiet vllm==0.31.0
fi
"$ENV/bin/python" -c "import vllm, torch; print('vllm', vllm.__version__, 'torch', torch.__version__, 'cuda', torch.version.cuda)"
"$ENV/bin/python" -c "import inspect, vllm.model_executor.models.registry as r; assert 'DeepseekV4ForCausalLM' in inspect.getsource(r); print('DeepseekV4ForCausalLM is in the model list')"
export HF_HUB_CACHE=$SCRATCH/.huggingface
"$ENV/bin/hf" download deepseek-ai/DeepSeek-V4-Flash
du -sh "$HF_HUB_CACHE/models--deepseek-ai--DeepSeek-V4-Flash"
