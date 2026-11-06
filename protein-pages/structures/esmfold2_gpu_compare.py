"""Compare a fold of test proteins on another GPU type with the same proteins on the H100.

    python esmfold2_gpu_compare.py <test folder> <H100 fused folder> <H100 default folder>

Each folder holds esmfold2_summary.tsv and <ID>.pdb from esmfold2_fold.py. For each protein of the
test folder's summary this reports, against the H100 fused-kernel fold: the CA RMSD (angstrom) after
the best superposition, over all residues and over the residues with H100 pLDDT of at least 70, the
change in mean pLDDT and pTM, and the seconds on each GPU. The H100 default-code fold gets the same
comparison, which shows how far two H100 folds of one protein are apart.

Writes <test folder>/esmfold2_gpu_compare.tsv (no residues) and prints whether the rules below pass.
These rules were set before the test ran. A chunk is at most 3.1 predicted H100 hours and the job
limit is 5 hours including 5 minutes of loading, so the other GPU may take at most 1.5 times as long.
"""

import csv
import sys
from pathlib import Path

import numpy as np

CONFIDENT = 70  # pLDDT
MAX_TIME_RATIO = 1.5
MAX_MEDIAN_RMSD = 0.5  # angstrom, over the residues with pLDDT of at least CONFIDENT
MAX_WORST_RMSD = 2.0
MAX_PLDDT_CHANGE = 5  # pLDDT points, for any one protein
COLUMNS = (
    "id",
    "length",
    "compared_with",
    "rmsd_all",
    "rmsd_confident",
    "n_confident",
    "plddt_change",
    "ptm_change",
    "seconds",
    "seconds_h100",
    "status",
)


def read_summary(folder):
    with open(folder / "esmfold2_summary.tsv") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def read_ca(path):
    """CA coordinates and pLDDT (B-factor column) of a .pdb file."""
    lines = [l for l in open(path) if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    xyz = np.array([[float(l[i : i + 8]) for i in (30, 38, 46)] for l in lines])
    return xyz, np.array([float(l[60:66]) for l in lines])


def rmsd(a, b):
    """RMSD of two equal-length point sets after the best rotation and shift (Kabsch)."""
    a, b = a - a.mean(0), b - b.mean(0)
    u, s, vt = np.linalg.svd(a.T @ b)
    s[-1] *= np.sign(np.linalg.det(u @ vt))  # a reflection is not allowed
    return float(np.sqrt(max(0.0, ((a**2).sum() + (b**2).sum() - 2 * s.sum()) / len(a))))


def number(text):
    return float(text) if text not in ("", None) else float("nan")


def compare(pid, test_row, test_dir, ref_row, ref_dir, name, rows):
    row = dict.fromkeys(COLUMNS, "")
    row.update(
        id=pid, length=test_row["length"], compared_with=name, seconds=test_row["seconds"], seconds_h100=ref_row["seconds"]
    )
    pdb = test_dir / (pid + ".pdb")
    if not pdb.exists() or not (ref_dir / (pid + ".pdb")).exists():
        row["status"] = test_row["ca_ca_median"] if test_row["ca_ca_median"] == "out of memory" else "no structure"
        rows.append(row)
        return
    a, plddt_a = read_ca(pdb)
    b, plddt_b = read_ca(ref_dir / (pid + ".pdb"))
    if len(a) != len(b):
        row["status"] = "different length"
        rows.append(row)
        return
    keep = plddt_b >= CONFIDENT
    row.update(
        rmsd_all="%.3f" % rmsd(a, b),
        n_confident=int(keep.sum()),
        status="ok",
        plddt_change="%.2f" % (number(test_row["mean_plddt"]) - number(ref_row["mean_plddt"])),
        ptm_change="%.3f" % (number(test_row["ptm"]) - number(ref_row["ptm"])),
    )
    if keep.sum() >= 3:
        row["rmsd_confident"] = "%.3f" % rmsd(a[keep], b[keep])
    rows.append(row)


def check(label, passed):
    print("  %s: %s" % ("PASS" if passed else "FAIL", label))
    return passed


def main():
    test_dir, fused_dir, default_dir = (Path(p) for p in sys.argv[1:4])
    test, fused, default = (read_summary(d) for d in (test_dir, fused_dir, default_dir))
    rows = []
    for pid, r in test.items():
        compare(pid, r, test_dir, fused[pid], fused_dir, "H100 fused", rows)
        compare(pid, default[pid], default_dir, fused[pid], fused_dir, "H100 default (reference)", rows)
    out = test_dir / "esmfold2_gpu_compare.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print("wrote", out)
    results = {}
    for name in ("H100 default (reference)", "H100 fused"):
        ok = [r for r in rows if r["compared_with"] == name and r["status"] == "ok"]
        confident = [float(r["rmsd_confident"]) for r in ok if r["rmsd_confident"]]
        change = [abs(float(r["plddt_change"])) for r in ok]
        results[name] = ok, confident, change
        print("%s: %d of %d proteins compared" % (name, len(ok), len(test)))
        if confident:
            print(
                "  RMSD over pLDDT >= %d: median %.3f, largest %.3f angstrom (%d proteins)"
                % (CONFIDENT, np.median(confident), max(confident), len(confident))
            )
            print("  RMSD over all residues: median %.3f angstrom" % np.median([float(r["rmsd_all"]) for r in ok]))
            print("  largest mean pLDDT change: %.2f" % max(change))
    ok, confident, change = results["H100 fused"]
    secs = sum(float(r["seconds"]) for r in ok)
    secs_h100 = sum(float(r["seconds_h100"]) for r in ok)
    ratio = secs / secs_h100 if secs_h100 else float("nan")
    memory = max((number(r["gpu_gb"]) for r in test.values() if r["gpu_gb"]), default=float("nan"))
    print("test GPU: %.0f s, H100: %.0f s, %.2f times as long; most GPU memory %.1f GB" % (secs, secs_h100, ratio, memory))
    print("rules:")
    tests = [
        check("all %d proteins folded" % len(test), len(ok) == len(test)),
        check("time at most %.1f times the H100 (%.2f)" % (MAX_TIME_RATIO, ratio), ratio <= MAX_TIME_RATIO),
    ]
    if confident:
        tests += [
            check(
                "median RMSD at most %.1f angstrom (%.3f)" % (MAX_MEDIAN_RMSD, np.median(confident)),
                np.median(confident) <= MAX_MEDIAN_RMSD,
            ),
            check(
                "largest RMSD at most %.1f angstrom (%.3f)" % (MAX_WORST_RMSD, max(confident)), max(confident) <= MAX_WORST_RMSD
            ),
            check("no mean pLDDT change above %d (%.2f)" % (MAX_PLDDT_CHANGE, max(change)), max(change) <= MAX_PLDDT_CHANGE),
        ]
    else:
        tests.append(check("structures to compare", False))
    print("all rules pass" if all(tests) else "at least one rule failed")


if __name__ == "__main__":
    main()
