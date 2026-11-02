"""The proteins of the ESMFold2 timing test, and their FASTA.

    python esmfold2_timing_sample.py pick ../full_inputs/proteome_lengths.tsv.gz > esmfold2_timing_ids.tsv
    python esmfold2_timing_sample.py fasta botryllus_new_genome_protein.fasta esmfold2_timing_ids.tsv > sample.fasta

pick: BHF (FUN008084) first, to load the model and warm up before anything is timed; then 10
proteins drawn at random (seed 0) from each length range below, in a shuffled order so that no
range always runs first; then three long proteins (about 2,500, 3,500 and 5,000 aa), last, to
find where a fold stops fitting in the GPU's 80 GB. A sequence is drawn once even when several
isoforms share it. The ranges cover the proteome: half its proteins are under 259 aa, 90% under
762 aa and 99% under 1,890 aa.

fasta: the sequences of those IDs, from the proteome FASTA, under the IDs the Biomni tables use
(FUN015821, then FUN015821.2 for the second isoform with that header, as make_proteome_lengths.py
names them). The proteome is unpublished, so this FASTA is written on Sherlock and never committed.
"""
import gzip
import random
import sys

RANGES = [(1, 150), (151, 300), (301, 500), (501, 800), (801, 1200), (1201, 2000)]
PER_RANGE = 10
LONG = [2500, 3500, 5000]
WARMUP = "FUN008084"


def pick(lengths_tsv):
    first = {}  # md5 -> (ID, length), the first protein with that sequence
    with gzip.open(lengths_tsv, "rt") as fh:
        next(fh)
        for line in fh:
            pid, length, md5 = line.rstrip("\n").split("\t")
            first.setdefault(md5, (pid, int(length)))
    prots = sorted(v for v in first.values() if v[0] != WARMUP)
    rng = random.Random(0)
    drawn = []
    for lo, hi in RANGES:
        drawn += [(p, n, "%d-%d" % (lo, hi)) for p, n in rng.sample([x for x in prots if lo <= x[1] <= hi], PER_RANGE)]
    rng.shuffle(drawn)
    longs = [min(prots, key=lambda x: (abs(x[1] - t), x[0])) + ("long",) for t in LONG]
    print("protein_id\tlength\tlength_range")
    for p, n, r in [(WARMUP, dict(first.values())[WARMUP], "warm-up")] + drawn + longs:
        print("%s\t%d\t%s" % (p, n, r))


def fasta(proteome, ids_tsv):
    want = [l.split("\t")[0] for l in open(ids_tsv).read().splitlines()[1:]]
    seqs, seen, cur = {}, {}, None
    for line in open(proteome):
        line = line.strip()
        if line.startswith(">"):
            gene = line[1:].split()[0].split("_")[0]
            seen[gene] = seen.get(gene, 0) + 1
            cur = gene if seen[gene] == 1 else "%s.%d" % (gene, seen[gene])
            seqs[cur] = []
        elif line and cur:
            seqs[cur].append(line)
    missing = [p for p in want if p not in seqs]
    if missing:
        sys.exit("not in %s: %s" % (proteome, ", ".join(missing)))
    for p in want:
        print(">%s\n%s" % (p, "".join(seqs[p]).rstrip("*")))


if __name__ == "__main__":
    {"pick": pick, "fasta": fasta}[sys.argv[1]](*sys.argv[2:])
