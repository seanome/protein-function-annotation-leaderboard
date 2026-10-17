#!/bin/bash
# Submit the kmerseek runs: per alphabet, one index job, then one search array (one task
# per query chunk) that starts when that index is built. Memory and time per alphabet come
# from settings.tsv (scripts/why/kmerseek_plan.py). A finished index or search is skipped,
# never overwritten, so running this twice costs nothing.
#   bash kmerseek_submit.sh                  every alphabet, full Swiss-Prot 2026_03
#   TEST=1 bash kmerseek_submit.sh           two alphabets, a 5_000-entry Swiss-Prot, the four proteins
set -euo pipefail
WHY=${WHY:-$SCRATCH/why-numbers/2026_03}
HERE=$WHY/kmerseek_scripts
SETTINGS=$HERE/settings.tsv
if [[ ${TEST:-0} == 1 ]]; then
    KDIR=$WHY/kmerseek_test; QDIR=$KDIR/queries
    mkdir -p "$QDIR"
    zcat "$WHY/inputs/four.fasta.gz" | awk '/^>/ { split($1, a, "|"); print ">" a[2]; next } { print }' > "$QDIR/chunk_0000.fasta"
    zcat "$WHY/inputs/uniprot_sprot.fasta.gz" | awk '/^>/ { n++ } n <= 5000' > "$KDIR/sprot_5000.fasta"
    zcat "$WHY/inputs/four.fasta.gz" >> "$KDIR/sprot_5000.fasta"
    gzip -f "$KDIR/sprot_5000.fasta"
    export SPROT_GZ=$KDIR/sprot_5000.fasta.gz
    ONLY=${TEST_ALPHABETS:-hp_pbotc_1st_ed2,protein20}
    NCHUNK=1
else
    KDIR=$WHY/kmerseek; QDIR=$WHY/kmerseek_queries; ONLY=${ONLY:-}
    NCHUNK=$(ls "$QDIR"/chunk_*.fasta | wc -l)
fi
mkdir -p "$KDIR"
cd "$KDIR"
tail -n +2 "$SETTINGS" | while IFS=$'\t' read -r alphabet ksize classes bits fmax penalty xdrop kaq nmeas peak hours smem stime imem itime cpus; do
    if [[ -n "$ONLY" && ",$ONLY," != *",$alphabet,"* ]]; then continue; fi
    if [[ ${TEST:-0} == 1 ]]; then imem=32; itime=2; smem=32; stime=1; fi
    vars=ALL,WHY=$WHY,KDIR=$KDIR,QDIR=$QDIR,ALPHABET=$alphabet,KSIZE=$ksize,PENALTY=$penalty,XDROP=$xdrop,KA_QUERIES=$kaq${SPROT_GZ:+,SPROT_GZ=$SPROT_GZ}
    ijob=$(sbatch --parsable --mem=${imem}G --time=${itime}:00:00 --export="$vars" "$HERE/kmerseek_index.sbatch")
    sjob=$(sbatch --parsable --dependency=afterok:$ijob --array=0-$((NCHUNK - 1)) --mem=${smem}G --time=${stime}:00:00 \
           --export="$vars" "$HERE/kmerseek_search.sbatch")
    printf '%s\tk%s\tindex %s (%s GB, %s h)\tsearch %s x%s (%s GB, %s h each)\n' \
        "$alphabet" "$ksize" "$ijob" "$imem" "$itime" "$sjob" "$NCHUNK" "$smem" "$stime"
done
