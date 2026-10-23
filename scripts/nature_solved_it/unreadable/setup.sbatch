#!/bin/bash
# Solvers unreadable, step 1 of 4: fix the UniProt release, fetch Swiss-Prot of that release,
# build this run's own BUSCO image, check both images run, and download the BUSCO lineage
# datasets once for every later task. The last line writes $RUN/setup_ok, which collect.py
# requires.
# Submitted by submit_all.sh, which sets RUN, REPO, SLURM ids and the log path. Not meant to be
# run alone.
#
# Why a BUSCO image of its own (2026-10-20): the first submission (setup job 48668645) failed
# with "exec /usr/local/bin/busco failed: a shared library is likely missing in the image".
# The shared sandbox $SCRATCH/apptainer-cache/quay.io-biocontainers-busco-6.1.0--pyhdfd78af_2.sandbox,
# which ran BUSCO in job 47961379 on 2026-10-15, had lost usr/lib64 and every empty folder
# (tmp, var/tmp, ...) on 2026-10-18 at 07:37-07:38 PDT, by a process not identified. Its
# binaries ask for /lib64/ld-linux-x86-64.so.2, which no longer exists, so nothing in the
# image can start. This job builds a fresh sandbox inside the run folder, from the same image
# reference, checks that `busco --version` works, and makes it read-only so that a clean-up
# cannot delete files from it.
#SBATCH --job-name=unread-setup
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
set -euo pipefail

: "${RUN:?}" "${REPO:?}"
HERE=$REPO/scripts/nature_solved_it/unreadable
# Same image reference as scripts/why/sherlock_labels.sbatch (a sandbox folder; building a
# .sif crashed in mksquashfs on Sherlock, job 47947485).
BUSCO_URI=docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2
BUSCO_IMG=$RUN/images/busco-6.1.0--pyhdfd78af_2.sandbox
MMSEQS_IMG=$SCRATCH/apptainer-cache/quay.io-biocontainers-mmseqs2@sha256-3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef.img
export APPTAINER_CACHEDIR=$SCRATCH/apptainer-cache/cache
# Unpacking a 2.4 GB image needs room; node-local /tmp may be small.
export APPTAINER_TMPDIR=$RUN/apptainer_tmp
command rm -f "$RUN/setup_ok"

mkdir -p "$RUN/inputs" "$RUN/busco_downloads" "$RUN/images" "$APPTAINER_TMPDIR"
/usr/bin/python3 -I "$HERE/fetch_inputs.py" --setup --out "$RUN/inputs"

# One setup at a time may touch the image (2026-10-20: two submissions 13 s apart both built
# it, and the second setup's mv onto the finished read-only image was refused).
exec 9> "$RUN/images/busco_build.lock"
flock -n 9 || { echo "refusing: another setup holds $RUN/images/busco_build.lock" >&2; exit 1; }
# A sandbox left by an earlier setup is reused only if BUSCO starts in it.
if [[ -d $BUSCO_IMG ]] && ! apptainer exec "$BUSCO_IMG" busco --version > /dev/null 2>&1; then
    echo "existing $BUSCO_IMG does not run; moving it aside" >&2
    chmod -R u+w "$BUSCO_IMG"
    mv "$BUSCO_IMG" "$BUSCO_IMG.broken.$(date -u +%Y%m%dT%H%M%SZ)"
fi
if [[ ! -d $BUSCO_IMG ]]; then
    if [[ -e $BUSCO_IMG.part ]]; then
        echo "refusing: $BUSCO_IMG.part exists (an unfinished or duplicate build); remove it and resubmit" >&2
        exit 1
    fi
    apptainer build --sandbox "$BUSCO_IMG.part" "$BUSCO_URI"
    if [[ -e $BUSCO_IMG ]]; then
        echo "refusing: $BUSCO_IMG appeared during the build; not replacing it" >&2
        exit 1
    fi
    mv "$BUSCO_IMG.part" "$BUSCO_IMG"
fi
run() { apptainer exec --bind /scratch "$BUSCO_IMG" "$@"; }
run busco --version > "$RUN/busco_version.txt"
grep -q "BUSCO 6.1.0" "$RUN/busco_version.txt"
apptainer exec --bind /scratch "$MMSEQS_IMG" mmseqs version > "$RUN/mmseqs_version.txt"
test -s "$RUN/mmseqs_version.txt"
# Read-only from here on: no later job, and no clean-up, can delete files from it.
chmod -R a-w "$BUSCO_IMG"

# PLAN.md, "BUSCO lineage, chosen by kingdom before any run".
for lineage in metazoa_odb12 fungi_odb12 viridiplantae_odb12 eukaryota_odb12; do
    run busco --download_path "$RUN/busco_downloads" --download "$lineage" --opt-out-run-stats
    test -d "$RUN/busco_downloads/lineages/$lineage"
done
command rm -rf "$APPTAINER_TMPDIR"
{
    echo "setup_job ${SLURM_JOB_ID:?}"
    echo "finished_utc $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "busco_image $BUSCO_IMG"
    echo "busco_version $(cat "$RUN/busco_version.txt")"
    echo "mmseqs_version $(cat "$RUN/mmseqs_version.txt")"
} > "$RUN/setup_ok"
echo "setup done"
