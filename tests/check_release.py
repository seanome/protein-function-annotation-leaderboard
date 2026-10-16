#!/usr/bin/env python3
"""Checks for the release table, its validator and the release dry run.

1. The release table built from data/ passes schema/results_row.json, keeps every score row of
   data/results.parquet, and carries the same value and n as the source row.
2. Its dollar column is the measured core-hours times the rate printed beside it.
3. The validator refuses: a missing column, a repeated key, a malformed digest, an interval that
   does not contain the value, dollars that are not hours x rates.
4. --submitted refuses a carried-over row (no digest, interval or run date) and accepts a
   complete one.
5. `release.py --dry-run` builds the files, sends nothing, and lists each file it would upload.
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import release_rows  # noqa: E402
import validate_rows  # noqa: E402

SCHEMA = json.loads((ROOT / "schema/results_row.json").read_text())
failures: list[str] = []


def check(cond: bool, what: str) -> None:
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        failures.append(what)


rows = release_rows.build(ROOT / "data")
src = pl.read_parquet(ROOT / "data/results.parquet")
errors, _ = validate_rows.validate(rows, SCHEMA, submitted=False)
check(not errors, f"release table passes the schema: {errors}")
check(rows.height == src.height, f"every score row kept: {rows.height:,} of {src.height:,}")
probe = src.filter(pl.col("experiment") == "qfo_region_midi_plus", pl.col("tool") == "hmmer3_phmmer",
                   pl.col("metric") == "fmax", pl.col("stratum_axis") == "all",
                   pl.col("truth_source") == "Pfam-A domain instances", pl.col("dedup_or_raw") == "raw",
                   pl.col("target_db") == "QfO zebrafish")
got = rows.filter(pl.col("experiment") == "qfo_region_midi_plus", pl.col("tool") == "hmmer3_phmmer",
                  pl.col("metric") == "fmax", pl.col("stratum_key") == "all",
                  pl.col("truth_source") == "Pfam-A domain instances", pl.col("dedup_or_raw") == "raw",
                  pl.col("target_db") == "QfO zebrafish")
check(probe.height == 1 and got.height == 1 and got["value"][0] == probe["value"][0]
      and got["n"][0] == probe["n_instances"][0],
      f"phmmer zebrafish Pfam Fmax and n carried over: {got.select('value', 'n').rows()}")
check(got["cpu_h_per_1000"][0] is not None and got["hardware"][0] is not None,
      f"phmmer on QfO has a measured cost: {got.select('cpu_h_per_1000', 'peak_rss_gb').rows()}")

# 2. dollars = hours x rates, recomputed here with this file's own arithmetic.
c = rows.filter(pl.col("usd_per_1000_from_rate").is_not_null()).head(1000)
diff = (c["usd_per_1000_from_rate"] - (c["cpu_h_per_1000"] * c["usd_rate_cpu_h"]
                                       + c["gpu_h_per_1000"] * c["usd_rate_gpu_h"])).abs().max()
check(c.height > 0 and diff < 1e-12, f"dollars = hours x rates on {c.height} rows (max difference {diff})")

# 3. Broken tables.
small = rows.head(50)
cases = {
    "missing column": small.drop("metric"),
    "repeated key": pl.concat([small.head(2), small.head(1)]),
    "malformed digest": small.with_columns(pl.lit("hmmer:3.4").alias("container_digest")),
    "interval misses value": small.with_columns(pl.lit(2.0).alias("ci95_low"), pl.lit(3.0).alias("ci95_high")),
    "dollars not hours x rates": rows.filter(pl.col("usd_per_1000_from_rate").is_not_null()).head(5)
        .with_columns((pl.col("usd_per_1000_from_rate") * 2).alias("usd_per_1000_from_rate")),
    "empty key column": small.with_columns(pl.lit(None, pl.Utf8).alias("tool")),
    "cost with no provenance": got.with_columns(pl.lit(None, pl.Utf8).alias("cost_provenance")),
    "unknown provenance": got.with_columns(pl.lit("estimated").alias("cost_provenance")),
}
for what, df in cases.items():
    e, _ = validate_rows.validate(df, SCHEMA, submitted=False)
    check(bool(e), f"validator refuses: {what} ({e[:1]})")

# 4. --submitted: a carried-over row fails, a complete row passes.
e, _ = validate_rows.validate(small, SCHEMA, submitted=True)
check(bool(e), "--submitted refuses a row with no digest, interval or run date")
full = got.with_columns(
    pl.lit("quay.io/biocontainers/hmmer@sha256:" + "a" * 64).alias("container_digest"),
    (pl.col("value") - 0.01).alias("ci95_low"), (pl.col("value") + 0.01).alias("ci95_high"),
    pl.lit("2026-10-16").alias("run_date"), pl.lit("3.4").alias("tool_version"))
e, _ = validate_rows.validate(full, SCHEMA, submitted=True)
check(not e, f"--submitted accepts a complete row: {e}")
e, _ = validate_rows.validate(full.with_columns(pl.lit("reported").alias("cost_provenance")), SCHEMA, submitted=True)
check(bool(e), "--submitted refuses a reported cost: a host-run row is measured")

# 5. Dry run: builds, sends nothing (no token needed), lists the uploads.
with tempfile.TemporaryDirectory() as tmp:
    env = {k: v for k, v in __import__("os").environ.items() if k not in ("ZENODO_TOKEN", "HF_TOKEN")}
    p = subprocess.run([sys.executable, str(ROOT / "scripts/release.py"), "--version", "v0-test", "--dry-run",
                        "--date", "2026-10-16", "--out", tmp], capture_output=True, text=True, env=env)
    out = p.stdout
    built = sorted(x.name for x in (Path(tmp) / "v0-test").iterdir())
    check(p.returncode == 0 and "nothing is sent" in out and "PUT  <bucket>/results.parquet" in out,
          f"dry run lists the uploads: rc {p.returncode} {p.stderr[-300:]}")
    check("results.parquet" in built and "README.md" in built and "results_row.schema.json" in built,
          f"dry run built the files: {built}")
    man = sorted(x.split()[0] for x in out.splitlines() if x.startswith("  ") and "sha256" in x)
    check(len(man) + 1 == len(built), f"every built file is listed: {len(man) + 1} listed, {len(built)} built")

# 6. The DIAMOND module as the worked example. tests/fixtures/annotator_stub/ holds the rows
#    and cost table of `make stub` for annotator_tools.nf (seanome/2024-kmerseek-analysis,
#    branch olgabot/submission-guide at 2be8c6b5, which merges G1's olgabot/annotator-tools
#    at fc5f1d8f): real DIAMOND and MMseqs2 output for 60 human queries against 120 yeast
#    proteins, with stub timings (all zero) and version "stub". The rows must build and pass
#    the schema, the image digest from their container_digest column must reach the release
#    row, and --submitted must still refuse them for the columns the stub cannot fill.
stub = release_rows.build(ROOT / "tests/fixtures/annotator_stub")
e, empty = validate_rows.validate(stub, SCHEMA, submitted=False)
check(not e and stub.height == 168, f"DIAMOND stub rows pass the schema: {stub.height} rows, {e}")
dia = stub.filter(pl.col("tool") == "diamond_besthit")
check(dia["container_digest"].unique().to_list()
      == ["quay.io/biocontainers/diamond@sha256:e02013280712f1295e59339dc6ea162adf8a6bbdd730e73a9c00700c03c3e556"],
      f"DIAMOND digest carried into container_digest: {dia['container_digest'].unique().to_list()}")
check(dia["cost_provenance"].unique().to_list() == ["measured"],
      f"DIAMOND search step gives a measured cost: {dia['cost_provenance'].unique().to_list()}")
e, _ = validate_rows.validate(stub, SCHEMA, submitted=True)
missing = sorted(x.split(":")[0] for x in e)
# peak_rss_gb is empty because Nextflow reports no memory for a task on macOS without a
# container; a Sherlock run fills it.
check(missing == ["ci95_high", "ci95_low", "peak_rss_gb", "run_date"],
      f"--submitted refuses the stub rows only for interval, peak memory and run date: {missing}")

if failures:
    print(f"{len(failures)} FAILURES")
    sys.exit(1)
print("all release checks passed")
