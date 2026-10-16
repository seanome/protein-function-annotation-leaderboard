#!/usr/bin/env python3
"""Build a release's results.parquet, one row per score in the shape of schema/results_row.json.

    python scripts/release_rows.py --data data --out build/release/results.parquet

Reads data/results.parquet (the scores) and data/compute.parquet (the Nextflow-trace costs).
Every score row is kept. Cost columns are filled for the QfO region runs where compute.parquet
has a search step for the same experiment, tool and setting; elsewhere they stay empty, never
estimated. The dollar column is the measured core-hours times a price-list rate, so it carries
the rate and the rate's source beside it.

Columns this source table does not have (container digest, 95% interval, run date) stay empty
for rows carried over from runs made before the submission path; scripts/validate_rows.py
counts them. A row the host runs from a submitted module fills them.
"""

import argparse
import json
import sys
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from build_site import BUILD_TOOL, HARDWARE, PRICE_FILE, PRICE_PAGE, R_CPU, R_GPU  # noqa: E402

RATE_SOURCE = (f"AWS on-demand, US East (N. Virginia), Linux, read 2026-10-15: c7i.large $0.08925/h for 2 vCPU, "
               f"g5.xlarge $1.006/h for 1 A10G GPU ({PRICE_PAGE}; {PRICE_FILE})")
COST_EXPERIMENTS = ["qfo_region_mini", "qfo_region_midi", "qfo_region_midi_plus"]
# The best-hit annotator runs (annotator_tools.nf): only the *BestHitSearch step is the search;
# placing domains and scoring are not part of what the tool costs.
ANNOTATOR_EXPERIMENTS = ["annotator_tools_midi_plus"]


def schema_columns() -> list[dict]:
    return json.loads((ROOT / "schema/results_row.json").read_text())["columns"]


def search_cost(comp: pl.DataFrame) -> pl.DataFrame:
    """One row per experiment, tool and setting: the search step's cost per relative per 1,000
    queries. Database builds are excluded, as on the Cost page."""
    # The same split into search and database build as the Cost page (build_site.load), so
    # it is only applied to the QfO region runs, whose step names that split was written for.
    is_build = (pl.col("setting") == "database build") | pl.col("step").str.contains("Index|Db|DB")
    qfo = pl.col("experiment").is_in(COST_EXPERIMENTS) & ~is_build
    annot = pl.col("experiment").is_in(ANNOTATOR_EXPERIMENTS) & pl.col("step").str.ends_with("BestHitSearch")
    s = comp.filter(qfo | annot, pl.col("n_queries").is_not_null(),
                    pl.col("n_target_dbs") > 0, pl.col("cpu_h").is_not_null())
    per = s.with_columns(
        (pl.col("cpu_h") / pl.col("n_target_dbs") * 1000 / pl.col("n_queries")).alias("cpu_h_per_1000"),
        (pl.col("gpu_h").fill_null(0.0) / pl.col("n_target_dbs") * 1000 / pl.col("n_queries")).alias("gpu_h_per_1000"),
    )
    return (per.group_by("experiment", "tool", "setting")
            .agg(pl.col("cpu_h_per_1000").sum(), pl.col("gpu_h_per_1000").sum(),
                 pl.col("peak_ram_gb").max().alias("peak_rss_gb"))
            .sort("experiment", "tool", "setting"))


def build(data: Path) -> pl.DataFrame:
    res = pl.read_parquet(data / "results.parquet")
    cost = search_cost(pl.read_parquet(data / "compute.parquet"))
    rows = (res.join(cost, on=["experiment", "tool", "setting"], how="left")
            .with_columns(
                pl.col("truth_version").alias("dataset_version"),
                pl.col("stratum_axis").alias("stratum_key"),
                pl.col("stratum_bin").alias("stratum_value"),
                pl.col("n_instances").alias("n"),
                # The annotator rows carry the image in their own container_digest column
                # (bin/annotator_tool_rows.py); PR 148's tables have none, so it stays empty.
                (pl.col("container_digest") if "container_digest" in res.columns
                 else pl.lit(None, pl.Utf8)).alias("container_digest"),
                pl.lit(None, pl.Float64).alias("ci95_low"),
                pl.lit(None, pl.Float64).alias("ci95_high"),
                pl.lit(None, pl.Utf8).alias("run_date"),
            )
            .with_columns(
                (pl.col("cpu_h_per_1000") * R_CPU + pl.col("gpu_h_per_1000") * R_GPU).alias("usd_per_1000_from_rate"),
                pl.when(pl.col("cpu_h_per_1000").is_not_null()).then(pl.lit(R_CPU)).alias("usd_rate_cpu_h"),
                pl.when(pl.col("cpu_h_per_1000").is_not_null()).then(pl.lit(R_GPU)).alias("usd_rate_gpu_h"),
                pl.when(pl.col("cpu_h_per_1000").is_not_null()).then(pl.lit(RATE_SOURCE)).alias("usd_rate_source"),
                pl.when(pl.col("cpu_h_per_1000").is_not_null()).then(pl.lit(HARDWARE)).alias("hardware"),
                pl.when(pl.col("cpu_h_per_1000").is_not_null()).then(pl.lit("measured")).alias("cost_provenance"),
            ))
    cols = [c["name"] for c in schema_columns()]
    key = json.loads((ROOT / "schema/results_row.json").read_text())["key"]
    return rows.select(cols).sort(key, nulls_last=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    df = build(a.data)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(a.out, compression="zstd", statistics=False)
    n_cost = df.filter(pl.col("cpu_h_per_1000").is_not_null()).height
    print(f"wrote {a.out}: {df.height:,} rows, {n_cost:,} with a measured cost")
    return 0


if __name__ == "__main__":
    sys.exit(main())
