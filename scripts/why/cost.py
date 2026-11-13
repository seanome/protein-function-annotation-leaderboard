"""Step 7: what labelling the two proteomes cost, from Slurm's own record of job 47961379.

On Sherlock, from the repo checkout (the date it was read goes in the first line):
    { echo "# read $(date -u +%FT%TZ) on $(hostname)"; \
      sacct -j 47961379 --format=JobID,JobName,Elapsed,AllocCPUS,TotalCPU,MaxRSS,State,NodeList --parsable2; } \
      > results/why/compute/sacct_47961379.txt
Then:
    python3 -I scripts/why/cost.py
Writes the step 7 rows of results/why/numbers.parquet and the sacct_47961379 entry of
results/why/sources.json.

Two different numbers per task, both kept:
  reserved core-hours = AllocCPUS x Elapsed / 3600: the cores held, busy or not.
  CPU-hours used      = TotalCPU / 3600: the time the cores spent computing.
Task 0 is MMseqs2 (human and oyster in one search), task 1 BUSCO human, task 2 BUSCO
oyster (scripts/why/sherlock_labels.sbatch). No prices here; the site's cost page has them.
"""

import datetime
import hashlib
import json
import pathlib
import re
import sys

import polars as pl

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from numbers_table import SOURCES, Numbers  # noqa: E402

JOB = "47961379"
SACCT = pathlib.Path(f"results/why/compute/sacct_{JOB}.txt")
STEP = 7
TASKS = {"0": ("mmseqs2", "MMseqs2 18.8cc5c against Swiss-Prot 2026_03", None),
         "1": ("busco", "BUSCO 6.1.0, metazoa_odb12", "human"),
         "2": ("busco", "BUSCO 6.1.0, metazoa_odb12", "oyster")}
QUERIES = {"human": 20652, "oyster": 30272}  # proteins in each one-protein-per-gene set (step 4)
# Task 0 searched both proteomes and the four proteins of step 6 in one MMseqs2 run, so its
# cost cannot be split by species; per 1,000 proteins it divides by all of them.
MMSEQS_QUERIES = sum(QUERIES.values()) + 4


def seconds(text):
    """Slurm times: [D-]HH:MM:SS, MM:SS.mmm or HH:MM:SS.mmm."""
    days = 0
    if "-" in text:
        d, text = text.split("-", 1)
        days = int(d)
    parts = [float(x) for x in text.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    h, m, s = parts
    return days * 86400 + h * 3600 + m * 60 + s


def gigabytes(text):
    """MaxRSS as sacct prints it (12345678K, 3.2G, or bytes with no letter). Slurm's K, M, G,
    T are powers of 1024 bytes; the result is in GB of 10^9 bytes."""
    if not text:
        return None
    m = re.fullmatch(r"([\d.]+)([KMGT]?)", text)
    if not m:
        raise ValueError(f"MaxRSS {text!r}")
    power = " KMGT".index(m.group(2) or " ")
    return float(m.group(1)) * 1024**power / 1e9


def read_sacct(path):
    lines = [x for x in path.read_text().splitlines() if x and not x.startswith("#")]
    t = pl.DataFrame([dict(zip(lines[0].split("|"), x.split("|"))) for x in lines[1:]])
    t = t.with_columns(task=pl.col("JobID").str.extract(rf"^{JOB}_(\d+)"), step=pl.col("JobID").str.extract(r"\.(\w+)$"))
    out = []
    for task, g in t.filter(pl.col("task").is_not_null()).group_by("task", maintain_order=True):
        main = g.filter(pl.col("step").is_null()).row(0, named=True)
        rss = [gigabytes(x) for x in g["MaxRSS"].to_list() if x]
        out.append({
            "task": task[0], "state": main["State"], "node": main.get("NodeList") or "not recorded",
            "alloc_cpus": int(main["AllocCPUS"]), "elapsed_s": seconds(main["Elapsed"]),
            "total_cpu_s": seconds(main["TotalCPU"]), "max_rss_gb": max(rss) if rss else None,
        })
    return pl.DataFrame(out)


def record_source():
    sources = json.loads(SOURCES.read_text())
    first = SACCT.read_text().splitlines()[0]
    sources[f"sacct_{JOB}"] = {
        "file": SACCT.name, "url": f"results/why/compute/{SACCT.name} (sacct on Sherlock, {first.lstrip('# ')})",
        "sha256": hashlib.sha256(SACCT.read_bytes()).hexdigest(), "release": "computed",
        "fetched_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    SOURCES.write_text(json.dumps(sources, indent=1) + "\n")


def main():
    record_source()
    t = read_sacct(SACCT)
    print(t.to_pandas().to_string())
    N = Numbers()
    src, rel = f"sacct_{JOB}", f"Sherlock job {JOB}"
    for r in t.iter_rows(named=True):
        tool, what, sp = TASKS[r["task"]]
        if r["state"] != "COMPLETED":
            raise ValueError(f"task {r['task']} is {r['state']}")
        pid = f"cost_{tool}" + (f"_{sp}" if sp else "")
        who = f"{what}, " + (f"{sp} proteome" if sp else "human and oyster proteomes and the four proteins in one search") + f" (task {r['task']})"
        reserved = r["alloc_cpus"] * r["elapsed_s"] / 3600
        used = r["total_cpu_s"] / 3600
        n_q = QUERIES[sp] if sp else MMSEQS_QUERIES
        N.add(f"{pid}_wall_h", STEP, f"{who}: wall-clock time", round(r["elapsed_s"] / 3600, 3), "hours", src, rel)
        N.add(f"{pid}_cpus", STEP, f"{who}: CPUs reserved", r["alloc_cpus"], "CPUs", src, rel, note=f"node {r['node']}")
        N.add(f"{pid}_reserved_core_h", STEP, f"{who}: reserved core-hours (CPUs reserved x wall-clock hours)", round(reserved, 3), "core-hours", src, rel)
        N.add(f"{pid}_used_cpu_h", STEP, f"{who}: CPU-hours used (Slurm TotalCPU)", round(used, 3), "CPU-hours", src, rel)
        N.add(f"{pid}_reserved_core_h_per_1000", STEP, f"{who}: reserved core-hours per 1,000 query proteins", round(reserved * 1000 / n_q, 4), "core-hours per 1,000 proteins", src, rel, note=f"{n_q:,} query proteins")
        N.add(f"{pid}_used_cpu_h_per_1000", STEP, f"{who}: CPU-hours used per 1,000 query proteins", round(used * 1000 / n_q, 4), "CPU-hours per 1,000 proteins", src, rel, note=f"{n_q:,} query proteins")
        if r["max_rss_gb"] is not None:
            N.add(f"{pid}_peak_rss_gb", STEP, f"{who}: peak memory (Slurm MaxRSS)", round(r["max_rss_gb"], 2), "GB", src, rel)
    N.write("results/why/numbers.parquet", replace_steps={STEP})
    print(N.frame().select("id", "value", "unit").to_pandas().to_string())


if __name__ == "__main__":
    main()
