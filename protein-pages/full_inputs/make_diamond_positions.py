"""Write diamond_positions.tsv: where each protein's best DIAMOND hit lies, by realignment.

    python make_diamond_positions.py botryllus_protein_annotations_v2.tsv.gz proteome_lengths.tsv.gz proteome.fasta uniprot_sprot.fasta.gz | gzip > diamond_positions.tsv.gz

Biomni's table names DIAMOND's best Swiss-Prot hit but keeps no positions. Each pair is aligned again
here (Biopython PairwiseAligner, local, BLOSUM62, gap open -11, extend -1), so the positions and
identity are this realignment's, not DIAMOND's. Positions are 1-based and inclusive. Only coordinates
and counts are written: the proteome is unpublished. A pair whose matrix would pass MAX_CELLS is
skipped, with blank positions, to stay within the laptop's memory.
"""
import csv
import gzip
import multiprocessing as mp
import sys
from pathlib import Path

from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import build_full  # noqa: E402

MAX_CELLS = 40_000_000
AL = None


def init():
    global AL
    AL = PairwiseAligner(mode="local", open_gap_score=-11, extend_gap_score=-1)
    AL.substitution_matrix = substitution_matrices.load("BLOSUM62")


def align(job):
    pid, acc, q, t = job
    if len(q) * len(t) > MAX_CELLS:
        return [pid, acc, "", "", "", "", "", ""]
    q2 = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in q)
    t2 = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in t)
    a = AL.align(q2, t2)[0]
    qb, tb = a.aligned
    same = sum(q2[i] == t2[j] for (q0, q1), (t0, t1) in zip(qb, tb) for i, j in zip(range(q0, q1), range(t0, t1)))
    cols = sum(q1 - q0 for q0, q1 in qb) + sum(qb[i + 1][0] - qb[i][1] + tb[i + 1][0] - tb[i][1] for i in range(len(qb) - 1))
    return [pid, acc, qb[0][0] + 1, qb[-1][1], tb[0][0] + 1, tb[-1][1], same, cols]


def main(ann, lengths, proteome, swissprot):
    rows = [r for r in csv.DictReader(gzip.open(ann, "rt"), delimiter="\t") if r["sp_acc"]]
    ls = {}
    with gzip.open(lengths, "rt") as fh:
        next(fh)
        for line in fh:
            p, n, m = line.rstrip("\n").split("\t")
            ls[p] = (int(n), m)
    seqs = build_full.read_proteome(proteome, ls)
    tseq = build_full.read_swissprot(swissprot, {r["sp_acc"] for r in rows})
    jobs = [(r["protein_id"], r["sp_acc"], seqs[r["protein_id"]], tseq[r["sp_acc"]]) for r in rows if r["protein_id"] in seqs and r["sp_acc"] in tseq]
    print("%d pairs of %d" % (len(jobs), len(rows)), file=sys.stderr)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["protein_id", "target", "qstart", "qend", "tstart", "tend", "identical", "columns"])
    with mp.get_context("spawn").Pool(4, init) as pool:
        for i, row in enumerate(pool.imap(align, jobs, chunksize=50)):
            w.writerow(row)
            if i % 2000 == 0:
                print(i, file=sys.stderr, flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:5])
