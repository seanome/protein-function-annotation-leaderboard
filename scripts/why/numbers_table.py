"""results/why/numbers.parquet: one row per number, with its unit and the file it came from."""

import json
import pathlib

import polars as pl

SOURCES = pathlib.Path("results/why/sources.json")


class Numbers:
    def __init__(self):
        self.rows = []
        self.sources = json.loads(SOURCES.read_text())

    def add(self, id, step, description, value, unit, source, release, note=None):
        if source not in self.sources:
            raise KeyError(f"{id}: source {source!r} is not in {SOURCES}")
        s = self.sources[source]
        self.rows.append(
            {
                "id": id,
                "step": step,
                "description": description,
                "value": float(value),
                "unit": unit,
                "source": source,
                "source_file": s["file"],
                "source_url": s["url"],
                "source_sha256": s["sha256"],
                "release": release,
                "fetched_utc": s["fetched_utc"],
                "note": note,
            }
        )

    def frame(self):
        return pl.DataFrame(self.rows)

    def write(self, path, replace_steps):
        """Replace this script's steps in the table and keep rows from the other steps."""
        new = self.frame()
        path = pathlib.Path(path)
        if path.exists():
            old = pl.read_parquet(path).filter(~pl.col("step").is_in(list(replace_steps)))
            new = pl.concat([old, new], how="diagonal_relaxed")
        dupes = new.filter(pl.col("id").is_duplicated())
        if dupes.height:
            raise ValueError(f"duplicate ids: {dupes['id'].to_list()}")
        new.sort(["step", "id"]).write_parquet(path)
