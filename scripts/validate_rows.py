#!/usr/bin/env python3
"""Check a results.parquet against schema/results_row.json.

    python scripts/validate_rows.py build/release/results.parquet
    python scripts/validate_rows.py --submitted rows_from_one_submission.parquet

Fails (exit 1) when a column is missing, extra or of the wrong type; when a column that may
not be empty has an empty cell; when the key columns repeat; when a value breaks its pattern,
enum or minimum; when ci95_low <= value <= ci95_high does not hold; or when the dollar column
is not the measured hours times the rates printed beside it.

With --submitted, every column marked submitted_required must be filled too: that is the bar
for a row the host ran from a submitted module. Without it, the empty cells in those columns
are counted and printed, since rows carried over from earlier runs have no digest, interval or
run date.
"""

import argparse
import json
import sys
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
TYPES = {"str": pl.Utf8, "f64": pl.Float64, "i64": pl.Int64}


def validate(df: pl.DataFrame, schema: dict, submitted: bool) -> tuple[list[str], dict]:
    errors: list[str] = []
    cols = schema["columns"]
    want = [c["name"] for c in cols]
    missing = [c for c in want if c not in df.columns]
    extra = [c for c in df.columns if c not in want]
    if missing:
        errors.append(f"missing columns: {missing}")
    if extra:
        errors.append(f"columns not in the schema: {extra}")
    if missing:
        return errors, {}
    empty_counts = {}
    for c in cols:
        name, s = c["name"], df[c["name"]]
        if s.dtype != TYPES[c["type"]] and not (s.dtype == pl.Null and c["nullable"]):
            errors.append(f"{name}: type {s.dtype}, schema says {c['type']}")
            continue
        nulls = s.null_count()
        if c.get("submitted_required"):
            empty_counts[name] = nulls
        if nulls and (not c["nullable"] or (submitted and c.get("submitted_required"))):
            errors.append(f"{name}: {nulls:,} empty cells; this column must be filled")
        v = s.drop_nulls()
        if "pattern" in c and v.len():
            bad = v.filter(~v.str.contains(c["pattern"]))
            if bad.len():
                errors.append(f"{name}: {bad.len():,} values break {c['pattern']}, e.g. {bad[0]!r}")
        if "enum" in c and v.len():
            bad = v.filter(~v.is_in(c["enum"]))
            if bad.len():
                errors.append(f"{name}: {bad.len():,} values not in {c['enum']}, e.g. {bad[0]!r}")
        if "min" in c and v.len() and v.min() < c["min"]:
            errors.append(f"{name}: minimum {v.min()} is below {c['min']}")
        if c["type"] == "f64" and v.len() and (v.is_nan() | v.is_infinite()).any():
            errors.append(f"{name}: NaN or infinite values")
    dup = df.group_by(schema["key"]).len().filter(pl.col("len") > 1)
    if dup.height:
        errors.append(f"{dup.height:,} key combinations repeat; key = {schema['key']}; "
                      f"first: {dup.row(0, named=True)}")
    ci = df.filter(pl.col("ci95_low").is_not_null() | pl.col("ci95_high").is_not_null())
    bad_ci = ci.filter(~((pl.col("ci95_low") <= pl.col("value")) & (pl.col("value") <= pl.col("ci95_high"))).fill_null(False))
    if bad_ci.height:
        errors.append(f"{bad_ci.height:,} rows where ci95_low <= value <= ci95_high does not hold")
    cost_cols = ["cpu_h_per_1000", "peak_rss_gb", "gpu_h_per_1000"]
    unlabelled = df.filter(pl.any_horizontal(pl.col(cost_cols).is_not_null()) & pl.col("cost_provenance").is_null())
    if unlabelled.height:
        errors.append(f"{unlabelled.height:,} rows have cost columns but no cost_provenance (measured or reported)")
    if submitted:
        rep = df.filter(pl.col("cost_provenance") != "measured")
        if rep.height:
            errors.append(f"{rep.height:,} rows have cost_provenance other than measured; a host-run row is measured")
    usd = df.filter(pl.col("usd_per_1000_from_rate").is_not_null())
    if usd.height:
        need = ["cpu_h_per_1000", "gpu_h_per_1000", "usd_rate_cpu_h", "usd_rate_gpu_h", "usd_rate_source"]
        no_rate = usd.filter(pl.any_horizontal(pl.col(need).is_null()))
        if no_rate.height:
            errors.append(f"{no_rate.height:,} rows have dollars without the hours, both rates and the rate source")
        recomputed = (pl.col("cpu_h_per_1000") * pl.col("usd_rate_cpu_h") + pl.col("gpu_h_per_1000") * pl.col("usd_rate_gpu_h"))
        off = usd.filter(((pl.col("usd_per_1000_from_rate") - recomputed).abs()
                          > 1e-9 * pl.max_horizontal(pl.lit(1.0), pl.col("usd_per_1000_from_rate").abs())).fill_null(True))
        if off.height:
            errors.append(f"{off.height:,} rows where dollars != hours x rates")
    return errors, empty_counts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("parquet", type=Path)
    ap.add_argument("--schema", type=Path, default=ROOT / "schema/results_row.json")
    ap.add_argument("--submitted", action="store_true", help="also require every submitted_required column")
    a = ap.parse_args()
    schema = json.loads(a.schema.read_text())
    df = pl.read_parquet(a.parquet)
    errors, empty = validate(df, schema, a.submitted)
    print(f"{a.parquet}: {df.height:,} rows, {df.width} columns, schema version {schema['version']}")
    if empty and not a.submitted:
        print("empty cells in columns a submitted row must fill (rows from earlier runs):")
        for k, v in empty.items():
            print(f"  {k:16s} {v:>10,} of {df.height:,}")
    for e in errors:
        print(f"ERROR {e}")
    print("ok" if not errors else f"{len(errors)} errors")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
