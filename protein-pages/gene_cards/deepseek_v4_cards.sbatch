#!/bin/bash
#SBATCH --job-name=cards-v4
#SBATCH --partition=hns,gpu
#SBATCH -G 4
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
# Gene cards from DeepSeek-V4-Flash, served by vLLM 0.31.0 (its container image, under
# Apptainer) on four H100 GPUs on this node,
# so the evidence sheets of the unpublished proteome never leave Sherlock. vLLM listens on
# 127.0.0.1 only and stops when the cards are written. Needs gene_cards/deepseek_v4_install.sbatch
# first (the container image and the 160 GB of weights).
#
#   sbatch gene_cards/deepseek_v4_cards.sbatch <evidence folder>
#
# Cards go to <evidence folder>/cards.DeepSeek-V4-Flash.jsonl; a rerun asks only for the missing
# proteins. Reasoning effort "high" (DeepSeek's Think High mode) and temperature 1.0, top_p
# 1.0, the sampling DeepSeek recommends for V4.
set -euo pipefail
EV=$(readlink -f "$1")
MODEL=deepseek-ai/DeepSeek-V4-Flash
NAME=$(basename "$MODEL")
HERE=$(dirname "$(readlink -f "$0")")
[ -f "$HERE/run_cards.py" ] || HERE=$SLURM_SUBMIT_DIR/gene_cards
SIF=$GROUP_HOME/olgabot/containers/vllm-openai_v0.31.0.sif
export HF_HUB_CACHE=$SCRATCH/.huggingface HF_HUB_OFFLINE=1 VLLM_NO_USAGE_STATS=1 DO_NOT_TRACK=1
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
PORT=$((20000 + SLURM_JOB_ID % 20000))
apptainer exec --nv --bind /scratch --env HF_HUB_CACHE="$HF_HUB_CACHE",HF_HUB_OFFLINE=1 "$SIF" \
  vllm serve "$MODEL" --served-model-name "$NAME" --host 127.0.0.1 --port $PORT \
  --tensor-parallel-size 4 --max-model-len 32768 --gpu-memory-utilization 0.92 \
  --reasoning-parser deepseek_v4 > "vllm-$SLURM_JOB_ID.log" 2>&1 &
VLLM=$!
trap 'kill $VLLM 2>/dev/null || true' EXIT
# Wait up to 40 minutes for the model to load, checking every 15 s
for i in $(seq 160); do
  if curl -sf "http://127.0.0.1:$PORT/health" > /dev/null; then break; fi
  kill -0 $VLLM 2>/dev/null || { echo "vLLM stopped; see vllm-$SLURM_JOB_ID.log" >&2; tail -30 "vllm-$SLURM_JOB_ID.log" >&2; exit 1; }
  sleep 15
done
curl -sf "http://127.0.0.1:$PORT/health" > /dev/null || { echo "vLLM did not start in 40 minutes" >&2; exit 1; }
echo "model loaded after $SECONDS s"
apptainer exec --bind /scratch "$SIF" python3 "$HERE/run_cards.py" "$EV" --server "http://127.0.0.1:$PORT" --model "$NAME" \
  --temperature 1 --top-p 1 --reasoning-effort high --batch 1 --parallel 16
echo "$(wc -l < "$EV/cards.$NAME.jsonl") cards in $EV/cards.$NAME.jsonl after $SECONDS s"
