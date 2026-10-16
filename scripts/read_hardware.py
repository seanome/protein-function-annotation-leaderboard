#!/usr/bin/env python3
"""Which machine ran each task that data/compute.parquet counts, read on Sherlock.

The run's Nextflow traces have no `hostname` column, so the trace cannot say which node a
task ran on. The trace does name the SLURM partition (`queue`) and the SLURM job id
(`native_id`). This script reads those, asks SLURM's job record (`sacct`) which node each
job ran on, and asks SLURM's node record (`scontrol show node`) which CPU and GPU that node
has. Sherlock documents the feature names at NODE_FEATURES_PAGE.

It picks the same task rows as seanome/2024-kmerseek-analysis
scripts/build_leaderboard_long.py (summarise_tasks): COMPLETED rows of the two run folders,
for the species the results scored, mapped to (tool, setting, step, target), the newest completion kept per task. It then checks
that its task count per (tool, setting, step) equals n_tasks in data/compute.parquet.

Needs `ssh sherlock` to work (Olga's control socket: `ssh -fN sherlock`). Light reads only:
awk over the traces, sacct, scontrol, all on the login node.

    python scripts/read_hardware.py --data data
writes data/hardware.parquet and data/hardware_provenance.json.
"""

import argparse
import datetime as dt
import io
import json
import re
import subprocess
import sys
from pathlib import Path

import polars as pl

EXPERIMENT = "qfo_region_midi_plus"
RUN_DIR = "/scratch/users/olgabot/2024-kmerseek-analysis/nextflow-runs/qfo-pfam-region-benchmark"
RUN_FOLDERS = ["run-midi-plus", "run-midi-plus-all-qfo"]
NODE_FEATURES_PAGE = "https://www.sherlock.stanford.edu/docs/advanced-topics/node-features/"
TECH_SPECS_PAGE = "https://www.sherlock.stanford.edu/docs/tech/"
# The run's traces span 2026-09-01 to 2026-09-13; sacct looks for jobs in this window.
SACCT_FROM, SACCT_TO = "2026-08-31", "2026-09-14"
SLURM_USER = "olgabot"

# Same maps as build_leaderboard_long.py (commit 4a4b1171), style "metrics".
QFO_SEARCH = {
    "phmmerSearch": ("hmmer3_phmmer", "default"),
    "jackhmmerSearch": ("hmmer3_jackhmmer", "n3"),
    "hhblitsSearch": ("hhblits", "single_seq"),
    "foldseekSearch": ("foldseek", "3di_aa"),
    "prostt5Search": ("prostt5", "3di_from_seq"),
    "reseekSearch": ("reseek", "default"),
    "hmmscanAnnotate": ("hmmscan", "direct"),
    "folddiscoQuery": ("folddisco", "default"),
}
QFO_DB = {"foldseekDb": "foldseek", "prostt5Db": "prostt5", "hhblitsBuildDB": "hhblits",
          "mmseqsDb": "mmseqs2", "reseekConvert": "reseek"}
KMERSEEK_TAG = re.compile(r"([a-z0-9]+)_(.+?)_k(\d+)(?:\.s(\d+))?_lc(false|true)")


def task_key(step: str, tag: str):
    tag = tag or ""
    if step in ("kmerseekSearch", "kmerseekIndex", "kmerseekIndexAndSearch"):
        m = KMERSEEK_TAG.fullmatch(tag)
        if not m:
            return None
        species, alphabet, k, scaled, lc = m.groups()
        s = f".s{scaled}" if scaled else ""
        return "kmerseek", f"{alphabet}_k{k}{s}_lc{lc.capitalize()}", species
    if step == "mmseqs2Search":
        m = re.fullmatch(r"human_vs_(\w+) \[(mmseqs2_seqseq|mmseqs2_iterative)\]", tag)
        return (m[2], "s7", m[1]) if m else None
    if step == "hmmscanAnnotate":
        return "hmmscan", "direct", "all"
    if step in QFO_SEARCH:
        m = re.fullmatch(r"human_vs_(\w+)", tag)
        return (*QFO_SEARCH[step], m[1]) if m else None
    if step in QFO_DB:
        return QFO_DB[step], "database build", tag
    return None


def ssh(cmd: str) -> str:
    p = subprocess.run(["ssh", "-o", "BatchMode=yes", "sherlock", cmd], capture_output=True, text=True)
    if p.returncode != 0:
        err = "\n".join(x for x in p.stderr.splitlines() if not x.startswith("Automated tools:"))
        raise SystemExit(f"ssh sherlock failed ({p.returncode}) on {cmd[:120]}...: {err.strip()[-600:]}\n"
                         "Open the socket with `ssh -fN sherlock` and run again.")
    return p.stdout


def read_trace_rows() -> pl.DataFrame:
    # Columns by name from each file's header, so a file with a different order still reads.
    awk = (r"""awk -F'\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;i++)c[$i]=i; next}"""
           r"""{print FILENAME,$c["native_id"],$c["process"],$c["tag"],$c["status"],$c["complete"],$c["queue"]}'""")
    files = " ".join(f"{f}/*trace*.txt" for f in RUN_FOLDERS)
    text = ssh(f"cd {RUN_DIR} && {awk} {files}")
    df = pl.read_csv(io.StringIO(text), separator="\t", has_header=False, infer_schema=False, quote_char=None,
                     new_columns=["trace_file", "native_id", "process", "tag", "status", "complete", "queue"])
    return df.with_columns(pl.col("process").str.split(":").list.last().alias("step"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data"))
    a = ap.parse_args()

    # The traces also hold searches of species the results never scored (the all-QfO folder
    # searched 77 targets). Keep the scored species, the pooled hmmscan step ("all") and the
    # query-side database builds ("human-<hash>"), as the table builder does.
    scored = set(pl.read_parquet(a.data / "results.parquet").filter(pl.col("experiment") == EXPERIMENT)
                 ["target_db"].str.strip_prefix("QfO ").unique().to_list()) | {"all"}
    tr = read_trace_rows().filter(pl.col("status") == "COMPLETED")
    rows = []
    for r in tr.iter_rows(named=True):
        k = task_key(r["step"], r["tag"])
        if k and (k[2] in scored or k[2].startswith("human")):
            rows.append({"tool": k[0], "setting": k[1], "step": r["step"], "target": k[2],
                         "complete": r["complete"], "native_id": r["native_id"], "partition": r["queue"],
                         "trace_file": f"{RUN_DIR}/{r['trace_file']}"})
    tasks = (pl.DataFrame(rows).sort("complete", nulls_last=True)
             .unique(["tool", "setting", "step", "target"], keep="last", maintain_order=True))
    # A tool's cost is kept only where the experiment scored it (Reseek and Folddisco searched
    # but were not scored); a build named after a tool family ("mmseqs2") belongs to its modes.
    tools = set(pl.read_parquet(a.data / "results.parquet").filter(pl.col("experiment") == EXPERIMENT)
                ["tool"].unique().to_list())
    keep = {t for t in tasks["tool"].unique().to_list() if any(x == t or x.startswith(t + "_") for x in tools)}
    tasks = tasks.filter(pl.col("tool").is_in(sorted(keep)))

    # Check: the same task rows as compute.parquet counts.
    comp = pl.read_parquet(a.data / "compute.parquet").filter(pl.col("experiment") == EXPERIMENT)
    mine = tasks.group_by("tool", "setting", "step").agg(pl.len().alias("n_mine"))
    cmp = comp.select("tool", "setting", "step", "n_tasks").join(mine, on=["tool", "setting", "step"], how="full",
                                                                coalesce=True)
    bad = cmp.filter(pl.col("n_tasks").fill_null(-1) != pl.col("n_mine").fill_null(-1))
    if bad.height:
        print(bad, file=sys.stderr)
        raise SystemExit(f"{bad.height} (tool, setting, step) groups differ from compute.parquet's n_tasks")

    # sacct: node and allocated resources of each job, in chunks.
    ids = sorted(set(tasks["native_id"].to_list()))
    acct = []
    # One call for every job of the run's user in the run's dates (about 1 s), filtered here.
    # Asking for the job ids with -j took about 3 minutes per 400 ids (measured 2026-10-16).
    want = set(ids)
    out = ssh(f"sacct -u {SLURM_USER} -X -n -P -S {SACCT_FROM} -E {SACCT_TO} "
              "--format=JobID,NodeList,Partition,AllocTRES")
    for line in out.splitlines():
        jid, node, part, tres = line.split("|")
        if jid in want:
            acct.append({"native_id": jid, "node": node or None, "sacct_partition": part, "alloc_tres": tres})
    acct = pl.DataFrame(acct).unique("native_id")

    # scontrol: each node's features (CPU_SKU, CPU_MNF, GPU_SKU ...).
    nodes = sorted({n for n in acct["node"].drop_nulls().to_list() if n != "None assigned"})
    # One node at a time: a node retired since the run is no longer in SLURM, and asking for
    # it fails the whole call. Such a node keeps a null CPU and GPU ("not recorded").
    out = ssh("for n in " + " ".join(nodes) + "; do scontrol -o show node $n 2>/dev/null; done; true")
    feats = []
    for line in out.splitlines():
        kv = dict(re.findall(r"(\w+)=(\S*)", line))
        if "NodeName" not in kv:
            print(f"skipped a line that names no node: {line[:120]}", file=sys.stderr)
            continue
        f = dict(x.split(":", 1) for x in kv.get("AvailableFeatures", "").split(",") if ":" in x)
        feats.append({"node": kv["NodeName"], "cpu_mnf": f.get("CPU_MNF"), "cpu_sku": f.get("CPU_SKU"),
                      "cpu_gen": f.get("CPU_GEN"), "gpu_brd": f.get("GPU_BRD"), "gpu_sku": f.get("GPU_SKU"),
                      "node_class": f.get("CLASS")})
    feats = pl.DataFrame(feats)

    hw = (tasks.join(acct, on="native_id", how="left").join(feats, on="node", how="left")
          .with_columns(
              # Sherlock's own spelling: manufacturer and SKU ("AMD 8224P", "INTEL 5118").
              pl.when(pl.col("cpu_sku").is_not_null())
              .then(pl.concat_str(pl.col("cpu_mnf"), pl.col("cpu_sku"), separator=" "))
              .alias("cpu_model"),
              pl.col("alloc_tres").str.extract(r"gres/gpu=(\d+)", 1).cast(pl.Int64).fill_null(0).alias("gpus_allocated"),
              pl.lit(None, dtype=pl.String).alias("trace_hostname"),
          )
          .with_columns(
              pl.when(pl.col("gpus_allocated") > 0).then(pl.col("gpu_sku")).alias("gpu_model"))
          .sort("tool", "setting", "step", "target"))
    hw.write_parquet(a.data / "hardware.parquet")
    prov = {
        "read_on": dt.date.today().isoformat(),
        "experiment": EXPERIMENT,
        "traces": f"{RUN_DIR}/{{{','.join(RUN_FOLDERS)}}}/*trace*.txt",
        "trace_has_hostname": False,
        "node_from": f"sacct -u {SLURM_USER} -X -n -P -S {SACCT_FROM} -E {SACCT_TO} "
                     "--format=JobID,NodeList,Partition,AllocTRES, matched on the trace's native_id",
        "hardware_from": "scontrol -o show node <node>: AvailableFeatures CPU_MNF, CPU_SKU, GPU_SKU",
        "feature_names_documented_at": NODE_FEATURES_PAGE,
        "partition_specs": TECH_SPECS_PAGE,
        "n_tasks": hw.height,
        "n_tasks_with_node": int(hw["node"].drop_nulls().len()),
        "n_nodes": len(nodes),
    }
    (a.data / "hardware_provenance.json").write_text(json.dumps(prov, indent=1) + "\n")
    print(json.dumps(prov, indent=1))
    print(hw.group_by("tool", "cpu_model", "gpu_model").len().sort("tool", "cpu_model"))


if __name__ == "__main__":
    main()
