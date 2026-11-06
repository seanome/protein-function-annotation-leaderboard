#!/bin/bash
# ESMFold2 on one RTX 3090: the 34 test proteins of 201 to 950 aa, to see whether short proteins of
# the proteome run could use this GPU. The H100 jobs are not touched.
#
#   cd structures && sbatch esmfold2_gpu_compare.sbatch
#
# Input: timing/sample.fasta and timing/{fused,default}/, the H100 folds of the same proteins with
# the same seed and settings. Output: gpu_compare/rtx3090/, with the .pdb and .cif files,
# esmfold2_summary.tsv and esmfold2_gpu_compare.tsv (esmfold2_gpu_compare.py).
#
# Sizes: these proteins took 9.7 minutes on the H100, and loading takes about 5, so 2 hours allows
# for a GPU 10 times slower. Shortest first, so a time-out loses the longest. GPU memory on the H100
# was at most 24 GB up to 983 aa, the RTX 3090's limit, so this stops at 950 aa. CPU memory as in
# esmfold2_proteome.sbatch. This GPU has no fp8 instructions, so the package runs that step in
# another number format, which is why speed and structures need testing.
#SBATCH -p hns,gpu
#SBATCH -A ayeletv
#SBATCH -J esmfold2-rtx3090
#SBATCH -G 1
#SBATCH -C GPU_SKU:RTX_3090
#SBATCH -c 8
#SBATCH --mem 96G
#SBATCH -t 2:00:00
#SBATCH -o %x-%j.out
set -euo pipefail
D=$GROUP_HOME/olgabot/envs/esmfold2  # the container and venv esmfold2_install.sbatch made
SIF=$D/pytorch-2.11.0-cuda12.8-cudnn9-runtime.sif
ENV=$D/venv-esm3.4.2
export HF_HOME=$SCRATCH/cache/hf HF_HUB_OFFLINE=1
run() { apptainer exec --nv -B /home/groups,/scratch "$SIF" "$@"; }  # --nv passes the GPU in
BASE=$SCRATCH/protein-pages-esmfold2
OUT=$BASE/gpu_compare/rtx3090
IDS=$(awk -F'\t' 'NR > 1 && $2 >= 201 && $2 <= 950' "$SLURM_SUBMIT_DIR/esmfold2_timing_fused_summary.tsv" | sort -t$'\t' -k2,2n | cut -f1)
if [ -z "$IDS" ]; then  # esmfold2_fold.py with no IDs would fold the whole FASTA
  echo "no test proteins" >&2
  exit 1
fi
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
set +e  # compare whatever finished, then exit with the fold's status
run "$ENV/bin/python" "$SLURM_SUBMIT_DIR/esmfold2_fold.py" --kernels fused "$BASE/timing/sample.fasta" "$OUT" $IDS
status=$?
set -e
run "$ENV/bin/python" "$SLURM_SUBMIT_DIR/esmfold2_gpu_compare.py" "$OUT" "$BASE/timing/fused" "$BASE/timing/default"
exit "$status"
