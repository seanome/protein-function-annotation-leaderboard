#!/usr/bin/env python3
"""Checks that the submission scorer is the leaderboard's scorer and grades the same way.

1. Every file in scorer/ has the sha256 that scorer/SOURCE.json records.
2. With --scored: rows a submission run wrote (scripts/score_submission.py --out) are
   compared with the leaderboard's own rows for the same tool, setting, relative, truth set,
   dedup choice, stratum and metric. Run on the MMseqs2 worked example, every shared row
   must agree, because the example repeats the leaderboard's own MMseqs2 search.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parent.parent
KEY = ["tool", "setting", "target_db", "truth_source", "dedup_or_raw", "stratum_axis",
       "stratum_bin", "metric"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scored", type=Path, help="folder written by score_submission.py --out")
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    ap.add_argument("--tol", type=float, default=1e-9)
    a = ap.parse_args()
    bad = 0
    src = json.loads((ROOT / "scorer/SOURCE.json").read_text())
    for name, meta in src["files"].items():
        got = hashlib.sha256((ROOT / "scorer" / name).read_bytes()).hexdigest()
        if got != meta["sha256"]:
            print(f"FAIL scorer/{name}: sha256 {got}, SOURCE.json says {meta['sha256']}")
            bad += 1
    print(f"scorer files: {len(src['files'])} checked, {bad} differ from SOURCE.json")
    if a.scored:
        new = pl.read_parquet(a.scored / "results.parquet")
        site = pl.read_parquet(a.data / "results.parquet").filter(
            pl.col("experiment") == "qfo_region_midi_plus")
        j = new.select(*KEY, "value", "n_instances").join(
            site.select(*KEY, "value", "n_instances"), on=KEY, how="inner", suffix="_site")
        diff = j.filter(((pl.col("value") - pl.col("value_site")).abs() > a.tol)
                        | (pl.col("n_instances") != pl.col("n_instances_site")))
        print(f"submission rows: {new.height}; also on the leaderboard: {j.height}; "
              f"differ by more than {a.tol}: {diff.height}")
        if j.height == 0:
            print("FAIL no submission row has a leaderboard row to compare with")
            bad += 1
        if diff.height:
            print(diff.head(10))
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
