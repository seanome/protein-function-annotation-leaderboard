#!/bin/bash
# Solvers unreadable, step 1 of 4: fix the UniProt release, fetch Swiss-Prot of that release,
# and download the BUSCO lineage datasets once for every later task.
# Submitted by submit_all.sh, which sets RUN, REPO and the log path. Not meant to be run alone.
#SBATCH --job-name=unread-setup
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
set -euo pipefail

: "${RUN:?}" "${REPO:?}"
HERE=$REPO/scripts/nature_solved_it/unreadable
IMG_DIR=$SCRATCH/apptainer-cache
# Same BUSCO image as scripts/why/sherlock_labels.sbatch (a sandbox folder; building a .sif
# crashed in mksquashfs on Sherlock, job 47947485).
BUSCO_URI=docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2
BUSCO_IMG=$IMG_DIR/quay.io-biocontainers-busco-6.1.0--pyhdfd78af_2.sandbox
export APPTAINER_CACHEDIR=$IMG_DIR/cache

mkdir -p "$RUN/inputs" "$RUN/busco_downloads"
/usr/bin/python3 -I "$HERE/fetch_inputs.py" --setup --out "$RUN/inputs"

flock "$BUSCO_IMG.lock" bash -c "[[ -d '$BUSCO_IMG/usr' ]] || apptainer build --sandbox '$BUSCO_IMG' '$BUSCO_URI'"
run() { apptainer exec --bind /scratch "$BUSCO_IMG" "$@"; }
run busco --version > "$RUN/busco_version.txt"
# PLAN.md, "BUSCO lineage, chosen by kingdom before any run".
for lineage in metazoa_odb12 fungi_odb12 viridiplantae_odb12 eukaryota_odb12; do
    run busco --download_path "$RUN/busco_downloads" --download "$lineage" --opt-out-run-stats
    test -d "$RUN/busco_downloads/lineages/$lineage"
done
echo "setup done"
