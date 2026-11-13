"""How long ESMFold2 would take on the whole proteome, read from the timing test.

    python esmfold2_timing_estimate.py esmfold2_timing_default_summary.tsv ../full_inputs/proteome_lengths.tsv.gz ../figures/esmfold2_timing

Seconds and GPU memory for a protein of a given length are read off the line between the test
proteins nearest in length (seconds on log-log axes), so nothing is assumed about the shape of the
curve. The longest length that fits on the 80 GB H100 comes from a quadratic through memory
against length. A sequence shared by several isoforms is counted once. Biomni already folded
proteins of 200 aa or less, so those are counted apart.
"""
import csv
import gzip
import sys

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent / "gene_cards"))
import pubfig as pf  # noqa: E402

GPU_GB = 80
BIOMNI_MAX_AA = 200
LOAD_S = 360  # about 6 minutes to load the model, once per job


def main(summary, lengths_tsv, out):
    done, oom = [], []
    for r in csv.DictReader(open(summary), delimiter="\t"):
        (oom if r["ca_ca_median"] == "out of memory" else done).append((int(r["length"]), float(r["seconds"]), float(r["gpu_gb"])))
    done = np.array(sorted(done[1:]))  # the first protein warms the model up
    n, secs, gb = done.T
    first = {}
    with gzip.open(lengths_tsv, "rt") as fh:
        next(fh)
        for line in fh:
            pid, length, md5 = line.rstrip("\n").split("\t")
            first.setdefault(md5, int(length))
    lengths = np.array(sorted(first.values()))
    # one point per length, for the interpolation
    un, idx = np.unique(n, return_inverse=True)
    us = np.exp(np.bincount(idx, np.log(secs)) / np.bincount(idx))
    # quadratic through memory, and the length where it reaches the GPU's memory
    q = np.polyfit(n, gb, 2)
    grid = np.arange(1, 6001)
    max_aa = int(grid[np.polyval(q, grid) <= GPU_GB].max())
    max_seen = int(n.max())
    est = np.exp(np.interp(np.log(lengths), np.log(un), np.log(us), left=np.log(us[0])))
    slope = np.polyfit(np.log(un[-8:]), np.log(us[-8:]), 1)
    over = lengths > un[-1]  # past the longest protein timed: follow the line of the last 8
    est[over] = np.exp(np.polyval(slope, np.log(lengths[over])))

    print("sequences %d (of %d protein IDs)" % (len(lengths), sum(1 for _ in gzip.open(lengths_tsv, "rt")) - 1))
    print("test proteins folded %d, out of memory %d (%s aa)" % (len(n) + 1, len(oom), ", ".join(str(o[0]) for o in oom)))
    print("seconds per length exponent over the last 8 lengths %.2f" % slope[0])
    print("longest length that fits in %d GB by the quadratic: %d aa; longest folded %d aa at %.0f GB" % (GPU_GB, max_aa, max_seen, gb.max()))
    rows = [("<= %d aa, Biomni has these" % BIOMNI_MAX_AA, lengths <= BIOMNI_MAX_AA),
            ("%d to %d aa, measured range" % (BIOMNI_MAX_AA + 1, max_seen), (lengths > BIOMNI_MAX_AA) & (lengths <= max_seen)),
            ("%d to %d aa, beyond the longest test, should fit" % (max_seen + 1, max_aa), (lengths > max_seen) & (lengths <= max_aa)),
            ("> %d aa, will not fit" % max_aa, lengths > max_aa)]
    print("group\tsequences\tGPU_hours")
    for name, m in rows:
        print("%s\t%d\t%.0f" % (name, m.sum(), est[m].sum() / 3600))
    todo = (lengths > BIOMNI_MAX_AA) & (lengths <= max_aa)
    print("to fold\t%d\t%.0f" % (todo.sum(), est[todo].sum() / 3600))
    for g in (1, 2, 4, 8):
        print("on %d GPUs at once: %.1f days" % (g, est[todo].sum() / 3600 / 24 / g))
    print("sum of seconds in the test %.0f" % (secs.sum()))

    pf.use_style()
    fig, (a, b) = pf.figure(pf.TWO_COLUMN_MM, 62, ncols=2)
    a.plot(n, secs / 60, "o", color=pf.OKABE_ITO["blue"], label="folded")
    for x, _, _ in oom:
        a.plot(x, 0.05, "x", color=pf.OKABE_ITO["vermillion"], label="out of memory" if x == oom[0][0] else None)
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xticks([100, 200, 500, 1000, 2000, 5000], ["100", "200", "500", "1000", "2000", "5000"])
    a.set_yticks([0.1, 1, 10, 30], ["0.1", "1", "10", "30"])
    a.set_xlabel("Protein length (aa)")
    a.set_ylabel("ESMFold2 time per protein (min)")
    a.legend(frameon=False)
    b.plot(n, gb, "o", color=pf.OKABE_ITO["blue"], label="folded")
    b.plot(grid[60:max_aa + 300], np.polyval(q, grid[60:max_aa + 300]), color=pf.GREY, label="quadratic through the points")
    b.axhline(GPU_GB, color=pf.OKABE_ITO["vermillion"], label="H100 memory, 80 GB")
    b.axvline(max_aa, color=pf.OKABE_ITO["vermillion"], ls=":", label="%d aa: predicted limit" % max_aa)
    b.set_xlabel("Protein length (aa)")
    b.set_ylabel("GPU memory used (GB)")
    b.legend(frameon=False)
    pf.panel_label(a, "a")
    pf.panel_label(b, "b")
    pf.save(fig, out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
