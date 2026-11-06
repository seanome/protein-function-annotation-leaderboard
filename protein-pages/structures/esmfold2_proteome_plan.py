"""Split the proteome into chunks for the ESMFold2 run, one chunk per Slurm array task.

    python esmfold2_proteome_plan.py esmfold2_timing_fused_summary.tsv ../full_inputs/proteome_lengths.tsv.gz esmfold2_proteome_chunks.tsv ../figures/esmfold2_proteome_plan

Takes each different sequence once, under the first ID that has it, and keeps lengths from 201 aa
(Biomni folded the shorter ones) to 2,070 aa (where the memory curve of the timing test reaches
the H100's 80 GB). They run shortest first, so most proteins are done early and the longest ones,
where the time and the fit are predictions, come last. A chunk ends once its predicted time
reaches 3 hours; the time per protein comes from the fused-kernel test (esmfold2_timing_estimate.py).
Writes protein_id, length, predicted_seconds, chunk. It holds IDs and lengths, no residues. The
figure shows how many sequences are done after each predicted GPU hour.
"""
import gzip
import sys

import numpy as np

from esmfold2_timing_estimate import pf, predict_seconds, read_timing

MIN_AA = 201
MAX_AA = 2070
CHUNK_S = 3 * 3600


def main(summary, lengths_tsv, out, figure):
    done, _ = read_timing(summary)
    n, secs, _ = done.T
    first = {}
    with gzip.open(lengths_tsv, "rt") as fh:
        next(fh)
        for line in fh:
            pid, length, md5 = line.rstrip("\n").split("\t")
            first.setdefault(md5, (pid, int(length)))
    prots = sorted(((length, pid) for pid, length in first.values() if MIN_AA <= length <= MAX_AA))
    est, _ = predict_seconds(n, secs, np.array([p[0] for p in prots]))
    chunk, used = 1, 0.0
    rows = []
    for (length, pid), s in zip(prots, est):
        if used >= CHUNK_S:
            chunk, used = chunk + 1, 0.0
        used += s
        rows.append((pid, length, s, chunk))
    with open(out, "w") as fh:
        fh.write("protein_id\tlength\tpredicted_seconds\tchunk\n")
        for pid, length, s, c in rows:
            fh.write("%s\t%d\t%.1f\t%d\n" % (pid, length, s, c))
    hours = np.bincount([r[3] for r in rows], [r[2] / 3600 for r in rows])[1:]
    print("%d sequences in %d chunks, %.0f predicted GPU hours; longest chunk %.2f h" % (len(rows), chunk, hours.sum(), hours.max()))
    print("chunk\tsequences\tshortest_aa\tlongest_aa\tpredicted_hours")
    for c in range(1, chunk + 1):
        got = [r for r in rows if r[3] == c]
        print("%d\t%d\t%d\t%d\t%.2f" % (c, len(got), got[0][1], got[-1][1], hours[c - 1]))
    gpu_hours = np.cumsum([r[2] for r in rows]) / 3600
    tested = np.array([r[1] for r in rows]) <= n.max()
    pf.use_style()
    fig, ax = pf.figure(pf.ONE_COLUMN_MM, 60)
    ax.plot(gpu_hours[tested], np.arange(1, len(rows) + 1)[tested], color=pf.OKABE_ITO["blue"], label="up to %d aa, the longest protein tested" % n.max())
    ax.plot(gpu_hours[~tested], np.arange(1, len(rows) + 1)[~tested], color=pf.OKABE_ITO["vermillion"], label="%d to %d aa, longer than any protein tested" % (n.max() + 1, MAX_AA))
    ax.set_xlabel("GPU hours into the run (predicted)")
    ax.set_ylabel("Sequences folded")
    ax.set_xlim(0)
    ax.set_ylim(0)
    ax.legend(frameon=False, loc="lower right")
    pf.save(fig, figure)


if __name__ == "__main__":
    main(*sys.argv[1:5])
