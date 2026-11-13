"""Predict each proteome's CPU-hours, wall time and memory from Slurm's record of job 47961379.

Job 47961379 (scripts/why/sherlock_labels.sbatch, 2026-10-15) is the only measured history:
  task 0, MMseqs2: 50,928 queries (human 20,652 + oyster 30,272 + 4 proteins), 16 CPUs;
  task 1, BUSCO metazoa_odb12 on human (20,652 proteins);
  task 2, BUSCO metazoa_odb12 on oyster (30,272 proteins).
MMseqs2 is scaled in proportion to the number of proteins. BUSCO is the straight line through
its two measured points. BUSCO on fungi_odb12, viridiplantae_odb12 and eukaryota_odb12 was
never run here, so it is predicted with the metazoa_odb12 line and marked unmeasured.

Download time (fetch.sbatch) uses the slowest UniProt paged download of the leaderboard's
fetch.py, on the laptop on 2026-10-15: oyster_interpro_xrefs, 57,874 rows between
14:12:39 and 14:18:25 UTC (results/why/sources.json), so 167 rows a second. Not measured on
Sherlock.

    python3 -I predict_resources.py --proteomes proteomes.tsv --sacct sacct_47961379.txt \
        --out resources.tsv            # table; prints the totals
    python3 -I predict_resources.py --proteomes proteomes.tsv --sacct sacct_47961379.txt --time-ask label
    python3 -I predict_resources.py --proteomes proteomes.tsv --sacct sacct_47961379.txt --time-ask fetch

Python 3.6 standard library, so submit_all.sh can call it on Sherlock.
"""

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect import sacct_by_task  # noqa: E402

JOB = 47961379
MMSEQS_QUERIES = 20652 + 30272 + 4  # results/why/SUMMARY.md: human and oyster sets, four.fasta
BUSCO_PROTEINS = {1: 20652, 2: 30272}  # task -> proteins (human, oyster)
FETCH_ROWS_PER_S = 57874 / (5 * 60 + 46)
# Proteins in the FASTA, counted by the leaderboard (results/why/SUMMARY.md). For human the
# README of 2026_03 lists 147,503 canonical entries, but the FASTA (replaced on 2026-09-15)
# holds 20,652; for every other proteome the README count is used.
FASTA_COUNTED = {"UP000005640": 20652, "UP000005408": 30272}
SAFETY = 3  # time ask = 3 x the largest prediction
MIN_ASK_H = 1.0


def read_rows(path):
    with open(path, newline="") as fh:
        delim = "\t" if path.endswith(".tsv") else ","
        return list(csv.DictReader(fh, delimiter=delim))


def model(sacct_path):
    t = sacct_by_task(sacct_path, JOB)
    mm = t[0]
    m = {
        "mmseqs_cpu_h_per_protein": mm["cpu_h"] / MMSEQS_QUERIES,
        "mmseqs_wall_h_per_protein": mm["wall_h"] / MMSEQS_QUERIES,
        "mmseqs_cpus": mm["cpus"],
        "peak_gib": max(x["maxrss_gib"] for x in t.values()),
    }
    (n1, n2) = BUSCO_PROTEINS[1], BUSCO_PROTEINS[2]
    for k in ("cpu_h", "wall_h"):
        slope = (t[2][k] - t[1][k]) / (n2 - n1)
        m["busco_" + k] = (t[1][k] - slope * n1, slope)  # intercept, slope
    return m, t


def predict(row, m):
    n = FASTA_COUNTED.get(row["upid"]) or int(row["readme_canonical_proteins"])
    lineages = list(dict.fromkeys([row["busco_lineage"], row["busco_lineage_all"]]))
    bi, bs = m["busco_cpu_h"]
    wi, ws = m["busco_wall_h"]
    mm_cpu = m["mmseqs_cpu_h_per_protein"] * n
    mm_wall = m["mmseqs_wall_h_per_protein"] * n
    b_cpu = sum(bi + bs * n for _ in lineages)
    b_wall = sum(wi + ws * n for _ in lineages)
    entries = int(row.get("uniprot_entries") or 0)
    return {
        "group": row["group"],
        "species": row["species"],
        "upid": row["upid"],
        "proteins": n,
        "uniprot_entries": entries,
        "busco_lineages": " + ".join(lineages),
        "busco_lineages_unmeasured": " + ".join(x for x in lineages if x != "metazoa_odb12"),
        "label_cpu_h": mm_cpu + b_cpu,
        "label_wall_h": mm_wall + b_wall,
        "fetch_wall_h": entries / FETCH_ROWS_PER_S / 3600,
    }


def ask(pred_h):
    return max(MIN_ASK_H, math.ceil(SAFETY * pred_h))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteomes", required=True, nargs="+")
    ap.add_argument("--sacct", required=True)
    ap.add_argument("--out")
    ap.add_argument("--time-ask", choices=("label", "fetch"))
    ap.add_argument("--concurrent", type=int, default=20, help="label tasks running at once (submit_all.sh's %%20)")
    args = ap.parse_args()
    m, measured = model(args.sacct)
    rows = [r for p in args.proteomes for r in read_rows(p) if r.get("upid")]
    preds = [predict(r, m) for r in rows]

    if args.time_ask:
        key = args.time_ask + "_wall_h"
        h = ask(max(p[key] for p in preds))
        if args.time_ask == "fetch":
            h = max(h, 2)  # download speed on Sherlock is unmeasured (laptop rate used)
        print("{:02d}:00:00".format(int(h)))
        return

    if args.out:
        with open(args.out, "w") as fh:
            w = csv.DictWriter(fh, list(preds[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            for p in preds:
                w.writerow({k: round(v, 4) if isinstance(v, float) else v for k, v in p.items()})
    tot_cpu = sum(p["label_cpu_h"] for p in preds)
    tot_wall = sum(p["label_wall_h"] for p in preds)
    print("measured, job {}: MMseqs2 {:.4f} CPU-h and {:.5f} wall-h per 1,000 proteins at {} CPUs; peak memory {:.2f} GiB".format(
        JOB, 1000 * m["mmseqs_cpu_h_per_protein"], 1000 * m["mmseqs_wall_h_per_protein"], m["mmseqs_cpus"], m["peak_gib"]))
    print("BUSCO metazoa_odb12 CPU-h = {:.4f} + {:.3e} x proteins; wall-h = {:.4f} + {:.3e} x proteins".format(*(m["busco_cpu_h"] + m["busco_wall_h"])))
    print("{} proteomes, {:,} proteins: {:.1f} CPU-h predicted for the label tasks ({:.1f} reserved core-hours at 16 CPUs)".format(
        len(preds), sum(p["proteins"] for p in preds), tot_cpu, 16 * tot_wall))
    print("label tasks, wall clock if {} run at once: {:.2f} h; largest task {:.2f} h, ask {} h".format(
        args.concurrent, tot_wall / args.concurrent, max(p["label_wall_h"] for p in preds), ask(max(p["label_wall_h"] for p in preds))))
    print("fetch tasks: largest {:.2f} h, total {:.2f} h one at a time, ask {} h".format(
        max(p["fetch_wall_h"] for p in preds), sum(p["fetch_wall_h"] for p in preds), ask(max(p["fetch_wall_h"] for p in preds))))
    print("BUSCO on any lineage but metazoa_odb12 is unmeasured; predicted with the metazoa_odb12 line")


if __name__ == "__main__":
    main()
