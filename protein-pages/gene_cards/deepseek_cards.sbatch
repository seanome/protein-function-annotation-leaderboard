#!/bin/bash
#SBATCH --job-name=cards-deepseek
#SBATCH --partition=hns,gpu
#SBATCH -G 2
#SBATCH -C GPU_MEM:48GB|GPU_MEM:80GB
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
# Gene cards from an open-weights DeepSeek model, served by vLLM on this GPU node, so the
# evidence sheets of the unpublished proteome never leave Sherlock. vLLM listens on 127.0.0.1
# only and stops when the cards are written.
#
#   sbatch gene_cards/deepseek_cards.sbatch <evidence folder> [model]
#
# <evidence folder> holds evidence.jsonl (make_cards.py); cards go to
# <evidence folder>/cards.<model name>.jsonl, and a rerun asks only for the missing proteins.
# The weights must already be in the Hugging Face cache the py-vllm module sets
# ($HF_HUB_CACHE, $SCRATCH/.huggingface); gene_cards/deepseek_download.sbatch puts them there.
#
# Test on 2026-10-31 (job 50515504 then 50515505, two L40S): model loaded in 287 s, 100 cards
# in 558 s in all, no failed call.
# DeepSeek-R1-Distill-Qwen-32B in bf16 is 65.5 GB: split over two GPUs of 48 GB (L40S) or
# 80 GB (H100), which leaves room for the model's working memory on 16 proteins at a time.
set -euo pipefail
EV=$(readlink -f "$1")
MODEL=${2:-deepseek-ai/DeepSeek-R1-Distill-Qwen-32B}
NAME=$(basename "$MODEL")
HERE=$(dirname "$(readlink -f "$0")")
[ -f "$HERE/run_cards.py" ] || HERE=$SLURM_SUBMIT_DIR/gene_cards
export HF_HUB_OFFLINE=1 VLLM_NO_USAGE_STATS=1 DO_NOT_TRACK=1
ml math py-vllm/0.19.1_py314
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
PORT=$((20000 + SLURM_JOB_ID % 20000))
vllm serve "$MODEL" --served-model-name "$NAME" --host 127.0.0.1 --port $PORT \
  --tensor-parallel-size 2 --max-model-len 32768 --gpu-memory-utilization 0.92 \
  --reasoning-parser deepseek_r1 > "vllm-$SLURM_JOB_ID.log" 2>&1 &
VLLM=$!
trap 'kill $VLLM 2>/dev/null || true' EXIT
# Wait up to 20 minutes for the model to load, checking every 15 s
for i in $(seq 80); do
  if curl -sf "http://127.0.0.1:$PORT/health" > /dev/null; then break; fi
  kill -0 $VLLM 2>/dev/null || { echo "vLLM stopped; see vllm-$SLURM_JOB_ID.log" >&2; tail -20 "vllm-$SLURM_JOB_ID.log" >&2; exit 1; }
  sleep 15
done
curl -sf "http://127.0.0.1:$PORT/health" > /dev/null || { echo "vLLM did not start in 20 minutes" >&2; exit 1; }
echo "model loaded after $SECONDS s"
python3 "$HERE/run_cards.py" "$EV" --server "http://127.0.0.1:$PORT" --model "$NAME" --batch 1 --parallel 16
echo "$(wc -l < "$EV/cards.$NAME.jsonl") cards in $EV/cards.$NAME.jsonl after $SECONDS s"
