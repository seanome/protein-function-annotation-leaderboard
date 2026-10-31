#!/bin/bash
#SBATCH --job-name=cards-deepseek-download
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
# Download an open-weights DeepSeek model from Hugging Face into the cache the py-vllm module
# sets ($HF_HUB_CACHE, $SCRATCH/.huggingface), for deepseek_cards.sbatch. 65.5 GB for the
# 32B model; 8.5 minutes on 2026-10-31 (job 50515504). Weights can be downloaded again, so they may live on purged scratch.
#
#   sbatch gene_cards/deepseek_download.sbatch [model]
set -euo pipefail
MODEL=${1:-deepseek-ai/DeepSeek-R1-Distill-Qwen-32B}
ml math py-vllm/0.19.1_py314
hf download "$MODEL" --include "*.json" "*.safetensors" "tokenizer*"
du -sh "$HF_HUB_CACHE"
