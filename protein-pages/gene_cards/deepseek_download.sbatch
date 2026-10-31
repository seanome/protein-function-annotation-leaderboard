#!/bin/bash
#SBATCH --job-name=cards-deepseek-download
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
# Download an open-weights DeepSeek model from Hugging Face into $SCRATCH/hf_models, for
# deepseek_cards.sbatch. Weights can be downloaded again, so they may live on purged scratch.
#
#   sbatch gene_cards/deepseek_download.sbatch [model]
set -euo pipefail
MODEL=${1:-deepseek-ai/DeepSeek-R1-Distill-Qwen-32B}
export HF_HOME=$SCRATCH/hf_models
ml math py-vllm/0.19.1_py314
hf download "$MODEL" --include "*.json" "*.safetensors" "tokenizer*"
du -sh "$HF_HOME"
