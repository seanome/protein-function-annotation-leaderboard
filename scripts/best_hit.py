"""The best-hit transfer group on the leaderboard page (Table 2), from data/best_hit/.

Source: seanome/2024-kmerseek-analysis PR 160 (head job 48547660 on Sherlock). Each human
query is searched against Swiss-Prot 2026_03 with every mammal removed; the Pfam domains of
its single best hit are copied onto the query, residue by residue, through the alignment.
That is a different call rule from Table 1's, where every hit's whole aligned region is the
call, so these rows are shown apart and never ranked against Table 1. The same best hits
scored with Table 1's rule (call_rules.parquet, call_rule "whole_region") sit beside them.

The target is Swiss-Prot without mammals, not one relative's proteome: no time since the
split from human, and Table 1's filters do not apply.

`values(data)` returns every number the table prints, keyed by setting; tests/check_site.py
recomputes them from the same parquet files with its own code.
"""

import html
import json
import shutil
from pathlib import Path

import polars as pl

PR160 = "https://github.com/seanome/2024-kmerseek-analysis/pull/160"
# Fixed order, as run: the table is not sorted by score.
SETTINGS = [
    ("diamond_besthit", "sensitive", "DIAMOND", "--sensitive"),
    ("diamond_besthit", "ultra-sensitive", "DIAMOND", "--ultra-sensitive"),
    ("mmseqs2_besthit", "s7.5", "MMseqs2", "-s 7.5, one round"),
]
SEARCH_STEP = {"diamond_besthit": "diamondBestHitSearch", "mmseqs2_besthit": "mmseqsBestHitSearch"}
FIGURES = ["fmax_best_hit_vs_qfo_relatives_midi_plus.png",
           "median_iou_best_hit_vs_qfo_relatives_midi_plus.png"]


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def values(data: Path) -> dict | None:
    d = data / "best_hit"
    if not (d / "rows.parquet").exists():
        return None
    rows = pl.read_parquet(d / "rows.parquet").filter(
        (pl.col("stratum_axis") == "all") & (pl.col("dedup_or_raw") == "raw")
        & (pl.col("truth_source") == "Pfam-A domain instances"))
    rules = pl.read_parquet(d / "call_rules.parquet")
    comp = pl.read_parquet(d / "compute.parquet")
    out = {}
    for tool, setting, _, _ in SETTINGS:
        r = rows.filter((pl.col("tool") == tool) & (pl.col("setting") == setting))
        m = {x["metric"]: x["value"] for x in r.iter_rows(named=True)}
        n = r["n_instances"].unique().to_list()
        if len(n) != 1:
            raise SystemExit(f"best_hit {tool} {setting}: truth instances {n}")
        rule = {x["call_rule"]: x for x in rules.filter(
            (pl.col("tool") == tool) & (pl.col("setting") == setting)).iter_rows(named=True)}
        if abs(rule["placed"]["fmax"] - m["fmax"]) > 1e-12:
            raise SystemExit(f"best_hit {tool} {setting}: call_rules placed Fmax "
                             f"{rule['placed']['fmax']} != rows {m['fmax']}")
        card = json.loads((d / "cards" / f"{tool}.{setting}.json").read_text())
        s = comp.filter((pl.col("tool") == tool) & (pl.col("setting") == setting)
                        & (pl.col("step") == SEARCH_STEP[tool])).row(0, named=True)
        out[f"{tool}|{setting}"] = {
            "fmax_placed": m["fmax"],
            "fmax_whole": rule["whole_region"]["fmax"],
            "precision": m["fmax_precision"],
            "recall": m["fmax_recall"],
            "iou_placed": m["median_iou_tp"],
            "n_tp_placed": rule["placed"]["n_tp_calls"],
            "iou_whole": rule["whole_region"]["median_iou_tp"],
            "n_tp_whole": rule["whole_region"]["n_tp_calls"],
            "n_truth": n[0],
            "queries_with_call": card["queries"]["n_queries_with_call"],
            "n_queries": card["queries"]["n_queries"],
            "cpu_h_per_1000": s["cpu_h"] * 1000 / s["n_queries"],
            "peak_ram_gb": s["peak_ram_gb"],
            "version": r["tool_version"].unique().to_list()[0],
        }
    shared = comp.filter(pl.col("tool") == "annotator_shared")
    out["_shared_cpu_h"] = float(shared["cpu_h"].sum())
    out["_shared_steps"] = sorted(shared["step"].to_list())
    # Database builds, per tool: left out of the rows' search core-hours too.
    db = comp.filter(pl.col("setting") == "database build")
    out["_db_cpu_h"] = {t: float(db.filter(pl.col("tool") == t)["cpu_h"].sum())
                        for t in ("diamond_besthit", "mmseqs2_besthit")}
    out["_n_queries"] = next(v["n_queries"] for k, v in out.items() if not k.startswith("_"))
    return out


def copy_figures(data: Path, out: Path) -> None:
    dest = out / "assets" / "best_hit"
    dest.mkdir(parents=True, exist_ok=True)
    for f in FIGURES:
        shutil.copyfile(data / "best_hit" / f, dest / f)


def section(v: dict | None, src: dict) -> str:
    if v is None:
        return ""
    commit = src.get("besthit_commit") or ""
    body = []
    for tool, setting, name, words in SETTINGS:
        x = v[f"{tool}|{setting}"]

        def td(k, text):
            return f'<td data-k="{k}" data-v="{esc(x[k])}">{text}</td>'

        body.append(
            f'<tr data-tool="{tool}" data-setting="{esc(setting)}">'
            f'<th scope="row">{name} <span class="muted">{esc(words)}</span></th>'
            + td("fmax_placed", f"<b>{x['fmax_placed']:.3f}</b>")
            + td("fmax_whole", f"{x['fmax_whole']:.3f}")
            + td("precision", f"{x['precision']:.3f}")
            + td("recall", f"{x['recall']:.3f}")
            + td("iou_placed", f"{x['iou_placed']:.3f} <span class=\"muted\">({x['n_tp_placed']:,})</span>")
            + td("iou_whole", f"{x['iou_whole']:.3f} <span class=\"muted\">({x['n_tp_whole']:,})</span>")
            + td("n_truth", f"{x['n_truth']:,}")
            + td("queries_with_call", f"{x['queries_with_call']:,} of {x['n_queries']:,}")
            + td("cpu_h_per_1000", f"{x['cpu_h_per_1000']:.3f}")
            + td("peak_ram_gb", f"{x['peak_ram_gb']:.1f}")
            + td("version", esc(x["version"]))
            + "</tr>")
    figs = "".join(
        f'<figure class="lbfig"><img src="assets/best_hit/{f}" alt="{esc(alt)}" '
        f'style="max-width:100%;height:auto" loading="lazy"><figcaption>{cap}</figcaption></figure>'
        for f, alt, cap in (
            (FIGURES[0],
             "Fmax per relative: phmmer, MMseqs2 and kmerseek against each QfO proteome, the three "
             "best-hit settings under both call rules, and the hmmscan ceiling",
             "<b>Figure 18.</b> Fmax, the same 2,435 Pfam domains in every row. Circles: every "
             "hit's whole aligned region (filled: duplicate calls kept; open: merged), the rule of "
             "Table 1. Triangles: the best hit's whole aligned region. Squares: the best hit's Pfam "
             "domains copied through the alignment. Only the squares jump; kmerseek's light marks "
             "are its best of 15 settings for each relative, picked after seeing the result."),
            (FIGURES[1],
             "Median IoU of correct calls for the same rows, colours and shapes",
             "<b>Figure 19.</b> Median IoU of correct calls: the overlap between a call and the "
             "true domain, divided by their union. Beside each row, the fewest and most correct "
             "calls its medians rest on. A tool that finds few domains can still score a high IoU, "
             "so read it next to Fmax. The squares sit at 1.0 because the hit proteins' domains "
             "were found with the same Pfam models that define the truth."),
        ))
    return f"""
<section class="besthit" id="best-hit" aria-labelledby="best-hit-h">
<h2 id="best-hit-h">A different call rule: the best Swiss-Prot hit's Pfam domains, copied through the alignment</h2>
<p class="prose">These three rows are not ranked against Table 1, because they make their calls
by a different rule. Each human protein is searched against Swiss-Prot 2026_03 with every mammal
removed. The Pfam domains of its single best hit are copied onto the query, residue by residue,
through the alignment. Table 1 uses the whole-region rule: for every hit, the whole aligned part of
the query is the call. The column "Fmax, whole aligned region" scores the same best hits with
Table 1's rule. The gap between a row's first two numbers is what the call rule adds.</p>
<p class="prose">The target is Swiss-Prot without mammals, not one relative's proteome. So these rows
have no time since the split from human, and the choices above Table 1 do not change them. They
score Pfam domains on every protein, with duplicate calls kept. The rows are in a fixed order,
not sorted by score. The copied domains' median IoU is 1.000 because the hit proteins' domains
were found with the same Pfam models that define the truth.
Source: <a href="{PR160}">PR 160</a>, copied at commit <code>{esc(commit[:8])}</code>.</p>
<figure class="lbfig">
<figcaption><b>Table 2.</b> Best Swiss-Prot hit (Mammalia removed), {v['_n_queries']:,} human proteins.</figcaption>
<div class="tablebox"><table class="lb" id="besthit-table"><thead><tr>
<th scope="col">Search</th>
<th scope="col">Fmax, domains copied through the alignment</th>
<th scope="col">Fmax, best hit's whole aligned region (Table 1's rule)</th>
<th scope="col">Precision at Fmax</th>
<th scope="col">Recall at Fmax</th>
<th scope="col">Median IoU, copied domains (correct calls)</th>
<th scope="col">Median IoU, best hit's whole aligned region (correct calls)</th>
<th scope="col">Pfam domains in the truth</th>
<th scope="col">Proteins with a call</th>
<th scope="col">Core-hours per 1,000 proteins, search alone</th>
<th scope="col">Peak memory of the search, GB</th>
<th scope="col">Version</th>
</tr></thead><tbody>
{''.join(body)}
</tbody></table></div>
</figure>
<p class="prose oneline">Three steps are shared by the three rows and left out of their core-hours:
building the reference, listing the hit proteins, and running hmmscan on the hit proteins. Together
they took <span id="besthit-shared" data-v="{v['_shared_cpu_h']}">{v['_shared_cpu_h']:.3f}</span>
core-hours for all {v['_n_queries']:,} proteins, not per 1,000. Building the search databases is left
out of the rows too:
<span id="besthit-db-diamond" data-v="{v['_db_cpu_h']['diamond_besthit']}">{v['_db_cpu_h']['diamond_besthit']:.4f}</span>
core-hours for DIAMOND and
<span id="besthit-db-mmseqs2" data-v="{v['_db_cpu_h']['mmseqs2_besthit']}">{v['_db_cpu_h']['mmseqs2_besthit']:.4f}</span>
for MMseqs2. Core-hours are cores asked times wall time, from the run's task records.</p>
{figs}
</section>
"""
