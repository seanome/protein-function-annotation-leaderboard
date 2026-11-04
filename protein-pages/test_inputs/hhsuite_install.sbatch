#!/bin/bash
# Get HH-suite and its two databases onto Sherlock, for hhsearch.sbatch.
#
#   cd <folder holding this file> && sbatch hhsuite_install.sbatch
#
# HH-suite is not a Sherlock module (`module spider hh-suite` finds nothing, 2026-11-01), so it
# runs from the BioContainers image hhsuite 3.3.0 (quay.io/biocontainers/hhsuite:3.3.0--h503566f_15),
# pulled by its digest into $SCRATCH/apptainer-cache next to the other images.
#
# Databases, into $SCRATCH/databases/hhsuite/:
#   UniRef30_2023_02/   UniRef30 2023_02 (Uniclust), for hhblits to build each protein's
#                       alignment; 70,555,922,504 bytes packed
#   pdb70/              PDB70 "from_mmcif_latest" (the 2022-03-13 build), the profiles hhsearch
#                       searches; 28,627,472,935 bytes packed
# Neither server folder publishes an md5 file (checked 2026-11-01), so the checks here are: the
# downloaded size equals the server's Content-Length, tar reads the whole archive without an
# error, every ffindex/ffdata file hhblits and hhsearch open is there and not empty, and, if the
# archive holds an md5sum file, `md5sum -c` on it. The tarballs are deleted only after all that.
# A rerun resumes a part-downloaded file (wget -c) and skips a database already unpacked.
#
# Sizes: ~99 GB to download. At 10-20 MB/s from the GWDG server that is 1.4-2.8 h, with the two
# downloads running side by side; unpacking (pigz + tar) to roughly 3x the packed size takes about
# 1 h more. 12 h leaves room for a server at ~5 MB/s; if it still runs out, resubmit and wget
# carries on where it stopped. Downloading and unpacking need little CPU or memory.
#SBATCH -p hns,normal
#SBATCH -A ayeletv
#SBATCH -J hhsuite-install
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -t 12:00:00
#SBATCH -o hhsuite-install-%j.out
set -euo pipefail
IMGS=$SCRATCH/apptainer-cache
HHSUITE=$IMGS/quay.io-biocontainers-hhsuite@sha256-4bf9bb5229de18f522a94f4443c19fdcbb0f0cb0e6ea92f5390aa170bcb0a24f.img
DB=$SCRATCH/databases/hhsuite
export APPTAINER_CACHEDIR=$IMGS/cache
mkdir -p "$DB"

if [ ! -s "$HHSUITE" ]; then
    apptainer pull "$HHSUITE" \
        docker://quay.io/biocontainers/hhsuite@sha256:4bf9bb5229de18f522a94f4443c19fdcbb0f0cb0e6ea92f5390aa170bcb0a24f
fi
apptainer exec "$HHSUITE" hhblits -h | head -1

# get NAME URL BYTES PREFIX: download, check, unpack into $DB/NAME, check the files named PREFIX_*
get() {
    local name=$1 url=$2 bytes=$3 prefix=$4
    local tgz=$DB/$(basename "$url") dir=$DB/$name
    if [ -f "$dir/.unpacked" ]; then
        echo "$name: already unpacked in $dir"
        return 0
    fi
    wget -c -nv --tries=20 --waitretry=60 -O "$tgz" "$url"
    local got
    got=$(stat -c %s "$tgz")
    if [ "$got" != "$bytes" ]; then
        echo "$name: downloaded $got bytes, the server said $bytes" >&2
        return 1
    fi
    echo "$name: downloaded $got bytes, unpacking"
    mkdir -p "$dir"
    pigz -dc -p 2 "$tgz" | tar -x -C "$dir"
    # some archives put the files in a subfolder: move them up so the prefix is $dir/$prefix
    if [ ! -e "$dir/${prefix}_cs219.ffindex" ]; then
        local sub
        sub=$(find "$dir" -mindepth 2 -maxdepth 3 -name "${prefix}_cs219.ffindex" -printf '%h\n' | head -1) || true
        if [ -n "$sub" ]; then mv "$sub"/* "$dir"/; fi
    fi
    local f
    for f in a3m hhm cs219; do
        for x in ffindex ffdata; do
            if [ ! -s "$dir/${prefix}_$f.$x" ]; then
                echo "$name: $dir/${prefix}_$f.$x is missing or empty" >&2
                return 1
            fi
        done
    done
    if [ -f "$dir/md5sum" ]; then
        (cd "$dir" && md5sum -c md5sum)
    fi
    ls -la "$dir"
    touch "$dir/.unpacked"
    rm -f "$tgz"
    echo "$name: done, tarball deleted"
}

get UniRef30_2023_02 https://wwwuser.gwdg.de/~compbiol/uniclust/2023_02/UniRef30_2023_02_hhsuite.tar.gz \
    70555922504 UniRef30_2023_02 > "$DB/UniRef30_2023_02.log" 2>&1 &
p1=$!
get pdb70 https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz \
    28627472935 pdb70 > "$DB/pdb70.log" 2>&1 &
p2=$!
s1=0; s2=0
wait $p1 || s1=$?
wait $p2 || s2=$?
cat "$DB/UniRef30_2023_02.log" "$DB/pdb70.log"
echo "UniRef30 exit $s1, PDB70 exit $s2"
[ "$s1" = 0 ] && [ "$s2" = 0 ]
du -sh "$DB"/*/
