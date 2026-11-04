"""One chunk of kmerseek matches on stdin; each query's best region on stdout. A second mode merges
the chunks and checks them against the run's own per-protein summary.

    zstdcat chunk.matches.csv.zst | python3 extract_kmerseek_best.py chunk > part.tsv
    python3 extract_kmerseek_best.py merge per_protein.tsv.gz part1.tsv part2.tsv ... | gzip > best_regions.tsv.gz

Best = lowest region E-value, then highest bit score. Region and target coordinates stay as kmerseek
wrote them, 0-based and end-exclusive. The query's residues are left out, because the proteome is
unpublished: the output holds coordinates and the Swiss-Prot residues only, so it can be committed.
"""
import csv
import gzip
import math
import sys

COLS = ["query_name", "target_name", "moltype", "ksize", "region_start", "region_end", "target_start", "target_end",
        "target_subseq", "region_length", "region_n_mismatches", "region_ka_bits", "region_evalue", "region_n_chained"]
csv.field_size_limit(10 ** 8)


def key(r):
    return float(r["region_evalue"]), -float(r["region_ka_bits"] or 0)


def chunk():
    best, bad = {}, 0
    for r in csv.DictReader(sys.stdin):
        try:
            e = float(r["region_evalue"])
        except ValueError:
            continue
        if not math.isfinite(e) or e > 10:
            continue
        if r["region_n_chained"] == "1" and len(r["region_subseq"]) != len(r["target_subseq"]):
            bad += 1
            continue
        q = r["query_name"]
        if q not in best or key(r) < key(best[q]):
            best[q] = r
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    for r in best.values():
        w.writerow([r[c] for c in COLS])
    print("%d one-piece regions with unequal strings" % bad, file=sys.stderr)


def merge(summary, parts):
    best = {}
    for p in parts:
        for f in csv.reader(open(p), delimiter="\t"):
            r = dict(zip(COLS, f))
            q = r["query_name"]
            if q not in best or key(r) < key(best[q]):
                best[q] = r
    same = diff = absent = 0
    for s in csv.DictReader(gzip.open(summary, "rt"), delimiter="\t"):
        r = best.get(s["query_name"])
        if not s["best_evalue"] or float(s["best_evalue"]) > 10:
            continue
        if r is None:
            absent += 1
        elif r["target_name"] == s["best_target"] and math.isclose(float(r["region_evalue"]), float(s["best_evalue"]), rel_tol=1e-4):
            same += 1
        else:
            diff += 1
    print("agrees with the summary: %d; differs: %d; missing: %d" % (same, diff, absent), file=sys.stderr)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(COLS)
    for q in sorted(best):
        w.writerow([best[q][c] for c in COLS])


if __name__ == "__main__":
    chunk() if sys.argv[1] == "chunk" else merge(sys.argv[2], sys.argv[3:])
