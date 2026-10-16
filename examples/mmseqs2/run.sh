#!/usr/bin/env bash
# The worked example of CONTRIBUTING.md: MMseqs2 on the human_midi_plus query set against one
# relative, written as a submission folder (predictions.tsv.gz, timing.tsv, tool.yaml).
#
#   examples/mmseqs2/run.sh <relative> <relative.fasta> <out folder> [threads]
#
# <relative> is a label from data/query_sets/human_midi_plus/relatives.tsv and
# <relative.fasta> its proteome from the QfO tarball (CONTRIBUTING.md, step 1 shows how to get
# it). Needs mmseqs (18.8cc5c is the version on the leaderboard) and GNU or BSD time.
set -euo pipefail

relative=$1
target_fasta=$2
out=$3
threads=${4:-8}
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
qs=$repo/data/query_sets/human_midi_plus
setting=s7

# The relative's proteome must be the exact file the leaderboard's truth tables were made for.
want=$(awk -F'\t' -v r="$relative" 'NR > 1 && $1 == r {print $8}' "$qs/relatives.tsv")
[ -n "$want" ] || { echo "no relative '$relative' in $qs/relatives.tsv" >&2; exit 1; }
got=$(shasum -a 256 "$target_fasta" | cut -d' ' -f1)
[ "$got" = "$want" ] || { echo "$target_fasta sha256 $got, expected $want" >&2; exit 1; }

mkdir -p "$out"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# Peak memory and wall time of one command. GNU time (Linux) prints kB with -v; BSD time
# (macOS) prints bytes with -l.
timed() {
    local log=$1; shift
    if /usr/bin/time -v true >/dev/null 2>&1; then
        /usr/bin/time -v -o "$log" "$@"
        awk -F': ' '/Maximum resident set size/ {printf "%.3f\n", $2 / 1024 / 1024}' "$log" > "$log.gb"
    else
        /usr/bin/time -l "$@" 2> "$log"
        awk '/maximum resident set size/ {printf "%.3f\n", $1 / 1024 / 1024 / 1024}' "$log" > "$log.gb"
    fi
}
now() { python3 -c 'import time; print(time.time())'; }

# Database build: the relative's proteome as an MMseqs2 database (once per relative).
t0=$(now)
timed "$work/db.time" mmseqs createdb "$target_fasta" "$work/target" >/dev/null
mmseqs createdb "$qs/queries.fasta" "$work/query" >/dev/null
t1=$(now)

# Search: the same flags as the leaderboard's MMseqs2 row (tools/mmseqs2_seqseq on the site).
timed "$work/search.time" mmseqs search "$work/query" "$work/target" "$work/result" "$work/tmp" \
    --threads "$threads" -s 7 --max-seqs 1000 -e 10.0 >/dev/null
mmseqs convertalis "$work/query" "$work/target" "$work/result" "$work/hits.tsv" \
    --threads "$threads" --format-output "query,target,qstart,qend,tstart,tend,bits,evalue" >/dev/null
t2=$(now)

# predictions.tsv: one row per hit. query_start/query_end become the call; the target and
# its interval let the scorer copy the target's known annotation onto the call. The bit
# score ranks the calls. No label: a search call carries the target, not a family name.
version=$(mmseqs version)
{
    printf 'query_set\trelative\tquery\tcall_start\tcall_end\ttarget\ttarget_start\ttarget_end\tlabel\tscore\ttool\tversion\tsetting\n'
    awk -F'\t' -v OFS='\t' -v rel="$relative" -v v="$version" -v s="$setting" \
        '{print "human_midi_plus", rel, $1, $3, $4, $2, $5, $6, "", $7, "mmseqs2_seqseq", v, s}' "$work/hits.tsv"
} | gzip -c > "$out/predictions.tsv.gz"

hardware=$( (sysctl -n machdep.cpu.brand_string 2>/dev/null || awk -F': ' '/model name/ {print $2; exit}' /proc/cpuinfo) | tr -d '\t')
{
    printf 'query_set\trelative\tsetting\tstep\twall_s\tcores\tpeak_ram_gb\tgpu_model\thardware\n'
    printf 'human_midi_plus\t%s\t%s\tdatabase_build\t%.1f\t1\t%s\tnone\t%s\n' \
        "$relative" "$setting" "$(echo "$t1 - $t0" | bc)" "$(cat "$work/db.time.gb")" "$hardware"
    printf 'human_midi_plus\t%s\t%s\tsearch\t%.1f\t%s\t%s\tnone\t%s\n' \
        "$relative" "$setting" "$(echo "$t2 - $t1" | bc)" "$threads" "$(cat "$work/search.time.gb")" "$hardware"
} > "$out/timing.tsv"

sed "s/@VERSION@/$version/; s/@THREADS@/$threads/" "$here/tool.yaml" > "$out/tool.yaml"
echo "wrote $out: $(gzip -dc "$out/predictions.tsv.gz" | tail -n +2 | wc -l | tr -d ' ') calls"
cat "$out/timing.tsv"
