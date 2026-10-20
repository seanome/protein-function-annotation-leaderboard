"""Solvers unreadable, step 4: one table row per proteome, with CPU-hours from Slurm.

Python 3.6 standard library. Reads every $RUN/out/<upid>/summary.json and the sacct text
written by collect.sbatch, and writes $RUN/unreadable_table.tsv. A proteome whose task
failed keeps its row with status "failed" and empty numbers, so nothing drops out unseen.

    python3 -I collect.py --proteomes proteomes.tsv --run $RUN --sacct $RUN/sacct.txt \
        --fetch-job 123 --label-job 124
"""

import argparse
import csv
import json
import os


def hours(t):
    """Slurm time ([D-]HH:MM:SS, MM:SS.mmm or HH:MM:SS) to hours."""
    days = 0
    if "-" in t:
        d, t = t.split("-", 1)
        days = int(d)
    parts = [float(x) for x in t.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    h, m, s = parts
    return days * 24 + h + m / 60 + s / 3600


def gib(rss):
    """MaxRSS ("4448664K") to GiB; Slurm's K, M and G are powers of 1024."""
    if not rss:
        return None
    unit = {"K": 1 / 1024 ** 2, "M": 1 / 1024, "G": 1.0}[rss[-1]]
    return float(rss[:-1]) * unit


def sacct_by_task(path, job):
    """{task index: {state, wall_h, cpu_h, maxrss_gib, cpus}} for one array job."""
    out = {}
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="|"))
    for r in rows:
        jid = r["JobID"]
        if not jid.startswith(str(job) + "_"):
            continue
        task, _, step = jid[len(str(job)) + 1:].partition(".")
        if not task.isdigit():
            continue
        t = out.setdefault(int(task), {})
        if step == "":
            t.update(state=r["State"], wall_h=hours(r["Elapsed"]), cpu_h=hours(r["TotalCPU"]), cpus=int(r["AllocCPUS"]))
        elif step == "batch":
            t["maxrss_gib"] = gib(r["MaxRSS"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteomes", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--sacct", required=True)
    ap.add_argument("--fetch-job", required=True)
    ap.add_argument("--label-job", required=True)
    args = ap.parse_args()
    rows = list(csv.DictReader(open(args.proteomes), delimiter="\t"))
    fetch = sacct_by_task(args.sacct, args.fetch_job)
    label = sacct_by_task(args.sacct, args.label_job)
    cols = [
        "group", "species", "upid", "proteome_taxid", "divergence_my", "status", "n_proteins",
        "n_interpro_match", "n_sprot_hit", "n_unlabelled", "unlabelled_pct", "n_uncharacterized_name",
        "uncharacterized_name_pct", "busco_lineage", "busco_completeness_pct", "busco_lineage_all",
        "busco_all_completeness_pct", "label_state", "label_wall_h", "label_cpu_h", "label_cpus",
        "label_maxrss_gib", "fetch_state", "fetch_wall_h", "fetch_cpu_h",
    ]
    out_rows = []
    for i, r in enumerate(rows):
        o = {k: r.get(k, "") for k in ("group", "species", "upid", "proteome_taxid", "divergence_my", "busco_lineage", "busco_lineage_all")}
        path = os.path.join(args.run, "out", r["upid"], "summary.json")
        if os.path.exists(path):
            s = json.load(open(path))
            o["status"] = "ok"
            for k in ("n_proteins", "n_interpro_match", "n_sprot_hit", "n_unlabelled", "unlabelled_pct", "n_uncharacterized_name", "uncharacterized_name_pct"):
                o[k] = s[k]
            by = {b["lineage"]: b["completeness_pct"] for b in s["busco"]}
            o["busco_completeness_pct"] = by.get(r["busco_lineage"], "")
            o["busco_all_completeness_pct"] = by.get(r["busco_lineage_all"], "")
        else:
            o["status"] = "failed"
        lt, ft = label.get(i, {}), fetch.get(i, {})
        o.update(
            label_state=lt.get("state", ""), label_wall_h=lt.get("wall_h", ""), label_cpu_h=lt.get("cpu_h", ""),
            label_cpus=lt.get("cpus", ""), label_maxrss_gib=lt.get("maxrss_gib", ""),
            fetch_state=ft.get("state", ""), fetch_wall_h=ft.get("wall_h", ""), fetch_cpu_h=ft.get("cpu_h", ""),
        )
        out_rows.append(o)
    dest = os.path.join(args.run, "unreadable_table.tsv")
    with open(dest, "w") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    failed = [o["upid"] for o in out_rows if o["status"] != "ok"]
    print(dest, len(out_rows), "proteomes,", len(failed), "failed:", " ".join(failed))


if __name__ == "__main__":
    main()
