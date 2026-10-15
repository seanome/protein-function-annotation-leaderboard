#!/usr/bin/env python3
"""Build the static site in docs/ from the tables in data/.

Every number on every page comes from data/results.parquet or data/compute.parquet. A
filter combination with no row has no entry in the JSON the pages read, and the pages show
"—" with "not run" for it. Nothing is averaged, filled or estimated here except
cost_usd_est, whose formula and prices are printed beside it.
"""

import argparse
import html
import json
import math
import re
import sys
from pathlib import Path

import markdown
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tools import KINDS, SPECIES, TOOLS  # noqa: E402

SITE = "Protein annotation leaderboard"
REPO = "seanome/protein-function-annotation-leaderboard"

# Default prices for cost_usd_est. Both from AWS's public on-demand price list for
# US East (N. Virginia), Linux, read 2026-10-15 from the file the price page loads
# (publication date 2026-10-12): c7i.large $0.08925 per hour for 2 vCPU, and g5.xlarge
# $1.006 per hour for 1 NVIDIA A10G GPU.
PRICE_PAGE = "https://aws.amazon.com/ec2/pricing/on-demand/"
PRICE_FILE = (
    "https://b0.p.awsstatic.com/pricing/2.0/meteredUnitMaps/ec2/USD/current/"
    "ec2-ondemand-without-sec-sel/US%20East%20(N.%20Virginia)/Linux/index.json"
)
R_CPU = 0.08925 / 2
R_GPU = 1.006

TRUTH = {"swissprot": "Swiss-Prot features", "pfam": "Pfam", "pfamn": "Pfam-N"}
DEFAULTS = {
    "truth": "pfam", "calls": "raw", "species": "zebrafish", "axis": "all",
    "stratum": "all", "metric": "fmax", "sort": "fmax", "dir": "desc",
    "r_cpu": R_CPU, "r_gpu": R_GPU,
}
DEFAULT_WORDS = "Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched"
# A stratum enters a card's comparisons and questions only with at least this many true
# instances; below it one or two calls swing Fmax from 0 to 1.
MIN_N = 30

IDENTITY = [("0-20%", "0–20%"), ("20-30%", "20–30%"), ("30-40%", "30–40%"),
            ("40-60%", "40–60%"), ("60-100%", "60–100%"),
            ("no_homolog", "no same-family domain")]
LENGTH_SHORT = [("1", "1 aa"), ("2-15", "2–15 aa"), ("16-30", "16–30 aa"), ("31-60", "31–60 aa")]
LENGTH_DOMAIN = [("61-120", "61–120 aa"), ("121-250", "121–250 aa"), ("251+", "251 aa or more")]
PLDDT = [("0-50", "0–50"), ("50-70", "50–70"), ("70-90", "70–90"), ("90-100", "90–100")]
FEATURES = {
    "ACT_SITE": "active site residue", "BINDING": "binding site residue",
    "COILED": "coiled coil", "DNA_BIND": "DNA-binding region", "DOMAIN": "domain",
    "INTRAMEM": "region inside the membrane that does not cross it", "MOTIF": "short motif",
    "REGION": "region of interest", "REPEAT": "repeat", "SITE": "other single-residue site",
    "TRANSMEM": "membrane-spanning segment", "ZN_FING": "zinc finger",
}
SPECIES_GROUPS = [
    ("mammals", ["ptroglodytes", "rnorvegicus", "clupus", "btaurus", "mouse", "mdomestica"]),
    ("other vertebrates", ["chicken", "xtropicalis", "zebrafish"]),
    ("invertebrates", ["ciona", "fly", "worm"]),
    ("fungus, plant, bacterium", ["yeast", "arabidopsis", "ecoli"]),
    ("no relative", ["all"]),
]
# What the setting names mean. kmerseek alphabets are listed in the kmerseek README
# (seanome/kmerseek), with their residue groups.
ALPHABETS = {
    "protein20": "the 20 amino acids, not reduced",
    "gbmr4": "4 letters: ADKERNTSQ, YFLIVMCWH, G, P (Solis and Rackovsky 2000)",
    "wwmj5": "5 letters: CMFILVWY, ATH, GP, DE, SNQRK (Wang and Wang 1999)",
    "hp_thomas_dill2": "2 letters, hydrophobic ACFILMVWY and polar DEGHKNPQRST (Thomas and Dill 1996)",
}
SETTING_WORDS = {
    ("hmmer3_phmmer", "default"): "phmmer's default options, E-value cutoff 10",
    ("hmmer3_jackhmmer", "n3"): "three rounds of search (-N 3)",
    ("mmseqs2_seqseq", "s7"): "sensitivity 7 (-s 7), one round",
    ("mmseqs2_iterative", "s7"): "sensitivity 7 (-s 7), three rounds",
    ("hhblits", "single_seq"): "profiles built from single proteins, no family database",
    ("foldseek", "3di_aa"): "3Di shape letters and amino acids together",
    ("prostt5", "3di_from_seq"): "3Di letters predicted from sequence by ProstT5",
    ("hmmscan", "direct"): "scans Pfam's models; no relative searched",
}


def setting_words(tool: str, setting: str) -> str:
    if (tool, setting) in SETTING_WORDS:
        return SETTING_WORDS[(tool, setting)]
    m = re.fullmatch(r"(.+)_k(\d+)_lc(False|True)", setting)
    if tool == "kmerseek" and m:
        lc = "low-complexity k-mers kept" if m[3] == "False" else "low-complexity k-mers removed"
        return f"alphabet {m[1]} ({ALPHABETS.get(m[1], 'see the kmerseek README')}), k = {m[2]}, {lc}"
    return setting


# Marker shape for each kmerseek alphabet on Figure 1; other tools take theirs from TOOLS.
KMERSEEK_SHAPES = {"protein20": "circle", "gbmr4": "square", "wwmj5": "triangle",
                   "hp_thomas_dill2": "diamond"}

METRICS = [
    {"id": "fmax", "label": "Fmax", "field": "fmax"},
    {"id": "auprc", "label": "AUPRC", "field": "auprc"},
    {"id": "iou", "label": "Median IoU", "field": "iou"},
]
STRATUM_WORDS = {
    "identity": dict(IDENTITY),
    "feature_length_bin": dict(LENGTH_SHORT + LENGTH_DOMAIN + [("whole_protein", "whole protein")]),
    "plddt": dict(PLDDT),
    "feature_type": {k: f"{k} ({v})" for k, v in FEATURES.items()},
}
AXIS_WORDS = {
    "identity": "identity to the closest same-family domain in the relative",
    "feature_length_bin": "length of the annotated region",
    "plddt": "AlphaFold confidence (mean pLDDT) of the human protein",
    "feature_type": "Swiss-Prot feature type",
}


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def chiprows(species_mya: dict) -> list[dict]:
    def c(v, label, title=None):
        d = {"value": v, "label": label}
        if title:
            d["title"] = title
        return d

    return [
        {"id": "truth", "kind": "plain", "label": "Truth set",
         "help": "the known annotations calls are graded against. Pfam: protein domain "
                 "families. Pfam-N: regions a deep-learning model added where Pfam's own "
                 "models found nothing. Swiss-Prot features: 12 kinds of curated feature, "
                 "such as active sites and membrane-spanning segments",
         "groups": [{"chips": [c("swissprot", "Swiss-Prot features", "curated features from Swiss-Prot entries"),
                               c("pfam", "Pfam", "Pfam-A domain families"),
                               c("pfamn", "Pfam-N", "regions added by a deep-learning model where Pfam-A's models found nothing (Pfam 35.0)")]}]},
        {"id": "calls", "kind": "plain", "label": "Dedup or raw",
         "help": "raw keeps every call; dedup merges calls of one family that overlap on the "
                 "same protein, so one domain hit by many relatives' proteins counts once",
         "groups": [{"chips": [c("raw", "raw"), c("dedup", "dedup")]}]},
        {"id": "species", "kind": "plain", "label": "Nearest reviewed relative: time since split",
         "help": "the species whose annotated proteome is searched and whose annotations are "
                 "copied onto the human proteins; the number is how many million years (Mya) "
                 "ago its line split from ours",
         "groups": [{"label": g, "chips": [
             c(s, SPECIES[s] if s == "all" else f"{SPECIES[s]} {int(species_mya[s])}")
             for s in members]} for g, members in SPECIES_GROUPS]},
        {"id": "identity", "kind": "stratum", "axis": "identity",
         "label": "Nearest reviewed relative: % identity",
         "help": "% identity of each human domain to the closest domain of the same family "
                 "in that relative",
         "groups": [{"chips": [c("all", "any")] + [c(v, w) for v, w in IDENTITY]}]},
        {"id": "length", "kind": "stratum", "axis": "feature_length_bin",
         "label": "Region type and length",
         "help": "length of the annotated region the call must find",
         "groups": [{"chips": [c("all", "any")]},
                    {"label": "short feature, 60 aa or less",
                     "chips": [c(v, w) for v, w in LENGTH_SHORT]},
                    {"label": "folded domain, 61 aa or more",
                     "chips": [c(v, w) for v, w in LENGTH_DOMAIN]},
                    {"chips": [c("whole_protein", "whole protein")]}]},
        {"id": "plddt", "kind": "stratum", "axis": "plddt", "label": "Disorder",
         "help": "mean AlphaFold confidence (pLDDT, 0–100) of the human protein; "
                 "low values mean disordered or unknown shape",
         "groups": [{"chips": [c("all", "any")] + [c(v, w) for v, w in PLDDT]}]},
        {"id": "ftype", "kind": "stratum", "axis": "feature_type", "label": "Feature type",
         "help": "Swiss-Prot feature key; only the Swiss-Prot truth set has these",
         "groups": [{"chips": [c("all", "any")] + [c(k, k, v) for k, v in FEATURES.items()]}]},
    ]


def mark_empty(rows: list[dict], res: pl.DataFrame) -> list[dict]:
    """Flag stratum chips that no row of the results table has, for any tool or filter."""
    present = {(a, v) for a, v in res.select("stratum_axis", "stratum").unique().iter_rows()}
    for row in rows:
        if row["kind"] != "stratum":
            continue
        for g in row["groups"]:
            for c in g["chips"]:
                if c["value"] != "all" and (row["axis"], c["value"]) not in present:
                    c["empty"] = True
    return rows


def layout(title: str, page: str, body: str, src: dict, depth: int = 0,
           head: str = "", description: str = "") -> str:
    up = "../" * depth
    nav = [("index.html", "Leaderboard", "leaderboard"), ("cost.html", "Score against cost", "cost"),
           ("pileup.html", "Region pile-up", "pileup"),
           ("tools/index.html", "Tool cards", "tools"), ("coverage.html", "Coverage", "coverage"),
           ("methods.html", "Methods", "methods")]
    links = "".join(
        f'<a href="{up}{h}"{" data-keep-filters" if p in ("leaderboard", "cost") else ""}'
        f'{" aria-current=\"page\"" if p == page else ""}>{t}</a>' for h, t, p in nav)
    commit = src["commit"]
    src_line = (
        f'Numbers from <a href="https://github.com/{src["repo"]}/tree/{commit}/{src["path"]}">'
        f'{src["repo"]}</a> at commit <a class="mono" href="https://github.com/{src["repo"]}/commit/{commit}">'
        f'{commit[:8]}</a> (<span class="mono">{src["path"]}</span>), copied {src["copied_at"][:10]}.'
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<meta name="description" content="{esc(description)}">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Source+Sans+3:wght@400;600;700&display=swap">
<link rel="stylesheet" href="{up}assets/site.css">
<link rel="icon" href="{up}assets/favicon.svg" type="image/svg+xml">
{head}
</head>
<body data-page="{page}">
<header class="site"><div class="wrap">
<a class="brand" href="{up}index.html">{SITE}</a>
<nav class="site" aria-label="Pages">{links}</nav>
</div></header>
<main class="wrap" data-page="{page}">
<p class="source" id="source-commit" data-commit="{commit}">{src_line}</p>
{body}
</main>
<footer class="site"><div class="wrap">
<p>{src_line}</p>
<p>Site source: <a href="https://github.com/{REPO}">{REPO}</a>. Code font: Fantasque Sans Mono 1.8.0 (SIL Open Font License).</p>
</div></footer>
</body>
</html>
"""


def legend_kinds(tools_used: set[str]) -> str:
    items = []
    for kind, k in KINDS.items():
        names = [TOOLS[t]["name"] for t in TOOLS if TOOLS[t]["kind"] == kind and t in tools_used]
        if names:
            items.append(f'<span class="item"><span class="kind k-{kind}" aria-hidden="true"></span>'
                         f'{esc(k["label"])}: {esc(", ".join(names))}</span>')
    return '<div class="legend" aria-label="Colour key">' + "".join(items) + "</div>"


def svg_marker(shape: str, color: str, filled: bool = True) -> str:
    fill = color if filled else "var(--surface)"
    if shape == "square":
        m = f'<rect x="3" y="3" width="10" height="10" fill="{fill}" stroke="{color}" stroke-width="1.6"/>'
    elif shape == "diamond":
        m = f'<path d="M8,1.5 L14.5,8 L8,14.5 L1.5,8 Z" fill="{fill}" stroke="{color}" stroke-width="1.6"/>'
    elif shape == "triangle":
        m = f'<path d="M8,2 L14,13 L2,13 Z" fill="{fill}" stroke="{color}" stroke-width="1.6"/>'
    else:
        m = f'<circle cx="8" cy="8" r="5.5" fill="{fill}" stroke="{color}" stroke-width="1.6"/>'
    return f'<svg width="16" height="16" viewBox="0 0 16 16" aria-hidden="true">{m}</svg>'


# --------------------------------------------------------------------------- data

EXPERIMENT = "qfo_region_midi_plus"
DATA_DIR = Path("data")
TRUTH_FROM_LONG = {"Pfam-A domain instances": "pfam", "Pfam-N domain instances": "pfamn",
                   "Swiss-Prot features": "swissprot"}
HARDWARE = "Stanford Sherlock, SLURM partition hns"


def load(data: Path):
    """Read the long tables (one row per experiment, tool, setting, target, truth, stratum and
    metric) and shape the leaderboard's experiment into the columns the pages use."""
    long = pl.read_parquet(data / "results.parquet")
    exp = long.filter(pl.col("experiment") == EXPERIMENT)
    res = (
        exp.filter(pl.col("metric").is_in(["fmax", "fmax_precision", "fmax_recall", "auprc", "median_iou_tp"]))
        .with_columns(
            pl.col("truth_source").replace_strict(TRUTH_FROM_LONG).alias("truth_set"),
            (pl.col("dedup_or_raw") == "dedup").alias("dedup"),
            pl.col("dedup_or_raw").alias("calls"),
            pl.when(pl.col("target_db") == "Pfam-A HMM library").then(pl.lit("all"))
            .otherwise(pl.col("target_db").str.replace("^QfO ", "")).alias("species"),
            pl.col("target_mya_to_human").alias("species_mya"),
            pl.col("stratum_bin").alias("stratum"),
            pl.col("n_instances").alias("n_truth_instances"),
        )
        .pivot(on="metric", index=["truth_set", "dedup", "calls", "tool", "setting", "species", "species_mya",
                                   "stratum_axis", "stratum", "n_truth_instances"], values="value")
        .with_columns(pl.lit(None, dtype=pl.Int64).alias("n_tp_calls"))
    )
    # Compute: one row per tool, setting and step, summed over the relatives searched.
    # cpu_h there is reserved core-hours (cores x wall time). The leaderboard shows the search
    # step per relative: cpu_h / n_target_dbs, per 1,000 queries.
    c = pl.read_parquet(data / "compute.parquet").filter(pl.col("experiment") == EXPERIMENT)
    search = c.filter(pl.col("setting") != "database build").filter(~pl.col("step").str.contains("Index|Db|DB"))
    build = c.filter((pl.col("setting") == "database build") | pl.col("step").str.contains("Index|Db|DB"))
    rows = []
    for r in search.iter_rows(named=True):
        b = build.filter((pl.col("tool") == r["tool"]) & pl.col("setting").is_in([r["setting"], "database build"]))
        rows.append({
            "tool": r["tool"], "setting": r["setting"], "step": r["step"],
            "n_relatives": r["n_target_dbs"], "n_queries": r["n_queries"], "cores": r["cores"],
            "wall_s": r["wall_s"], "cpu_h_total": r["cpu_h"],
            "cpu_h_per_1000_queries": r["cpu_h"] / r["n_target_dbs"] * 1000 / r["n_queries"],
            "peak_rss_gb": r["peak_ram_gb"],
            "gpu": "none" if r["gpu_type"] in (None, "none") else r["gpu_type"],
            "gpu_h_per_1000_queries": 0.0 if r["gpu_h"] is None else r["gpu_h"] / r["n_target_dbs"] * 1000 / r["n_queries"],
            "db_cpu_h": b["cpu_h"].sum() if b.height else None,
            "db_n": int(b["n_target_dbs"].sum()) if b.height else None,
            "db_gpu_h": b["gpu_h"].drop_nulls().sum() if b.height and b["gpu_h"].drop_nulls().len() else None,
            "db_gpu_type": next((g for g in b["gpu_type"].to_list() if g and g != "none"), None) if b.height else None,
            "hardware": HARDWARE, "cpu_h_reported": None, "n_completed_runs": 1,
        })
    comp_ts = pl.DataFrame(rows, infer_schema_length=None)
    # One row per relative this tool and setting searched, all carrying the same mean, so the
    # pages can key cost by relative; a relative it did not search gets no row ("not run").
    searched = res.select("tool", "setting", "species").unique()
    comp = searched.join(comp_ts, on=["tool", "setting"], how="inner")
    src = json.loads((data / "SOURCE.json").read_text())
    prov = json.loads((data / "provenance.json").read_text())
    prov["n_queries"] = int(comp_ts["n_queries"].max())
    cov = (data / "coverage.md").read_text()
    return res, comp, src, prov, cov


def row_order(res: pl.DataFrame) -> list[dict]:
    kinds = list(KINDS)
    pairs = res.select("tool", "setting").unique().to_dicts()
    pairs.sort(key=lambda r: (kinds.index(TOOLS[r["tool"]]["kind"]), list(TOOLS).index(r["tool"]),
                              r["setting"]))
    out = []
    for i, p in enumerate(pairs):
        alpha = re.fullmatch(r"(.+)_k\d+_lc(?:False|True)", p["setting"])
        shape = KMERSEEK_SHAPES.get(alpha[1], "circle") if p["tool"] == "kmerseek" and alpha else TOOLS[p["tool"]]["shape"]
        out.append(dict(p, order=i, words=setting_words(p["tool"], p["setting"]), shape=shape))
    return out


def results_json(res: pl.DataFrame) -> dict:
    out: dict[str, dict] = {}
    for r in res.iter_rows(named=True):
        key = "|".join([r["truth_set"], r["calls"], r["species"], r["stratum_axis"], r["stratum"]])
        k2 = f'{r["tool"]}|{r["setting"]}'
        if k2 in out.setdefault(key, {}):
            raise SystemExit(f"two rows for {key} {k2}")
        out[key][k2] = [r["fmax"], r["auprc"], r["median_iou_tp"], r["n_truth_instances"], r["n_tp_calls"]]
    return out


def compute_json(comp: pl.DataFrame) -> dict:
    out = {}
    for r in comp.iter_rows(named=True):
        rep = r["cpu_h_reported"]
        out[f'{r["tool"]}|{r["setting"]}|{r["species"]}'] = {
            "cpu_h_per_1000": r["cpu_h_per_1000_queries"],
            "cpu_h_reported_per_1000": None if rep is None else rep * 1000 / r["n_queries"],
            "peak_rss_gb": r["peak_rss_gb"],
            "gpu": r["gpu"],
            "gpu_h_per_1000": r["gpu_h_per_1000_queries"],
            "hardware": r["hardware"],
            "n_relatives": r["n_relatives"],
        }
    return out


# --------------------------------------------------------------------------- pages

def page_leaderboard(meta_json: str, src: dict, n_queries: int, tools_used: set[str], dd: dict) -> str:
    body = f"""
<h1>{SITE}</h1>
<p class="lede">Which tool finds the known domains and features of a protein, and what does it
cost? Each row is one tool at one setting. It searched {n_queries:,} human proteins against the
proteome of one relative and copied the relative's annotations onto them. The scores say how many
of the human proteins' known annotations it put in the right place. The queries are human proteins,
964 of them from chromosome 6; how the ranking holds for queries from other species is not
measured.</p>
<details class="prose"><summary>How to choose between Fmax and AUPRC, and raw against dedup</summary>
<p>Fmax is the score at the single best cutoff on the tool's own score, protein by protein: use it
when you will keep every call above one cutoff. AUPRC sums precision over every cutoff, counting
each call: use it when you will rank calls and read down the list. The two can rank tools in
opposite orders, because AUPRC pools every call of every protein while Fmax averages protein by
protein. A domain reported many times counts against both. Dedup merges those copies; for phmmer against zebrafish (Pfam) it moves Fmax from
{dd['raw_fmax']:.3f} to {dd['dedup_fmax']:.3f} and AUPRC from {dd['raw_auprc']:.3f} to
{dd['dedup_auprc']:.3f}. Definitions: <a href="methods.html">Methods</a>.</p></details>
<form class="filters" id="filters-form" onsubmit="return false" aria-label="Filters">
<div id="filters"></div>
<p class="axisnote">The benchmark split the proteins along one of the last four rows at a time,
so picking a value in one of them sets the other three back to "any".</p>
</form>
<div class="prices" aria-label="Prices used for the cost column">
<label>r_cpu, USD per CPU-hour <input id="r_cpu" type="number" min="0" step="0.001" inputmode="decimal"></label>
<label>r_gpu, USD per GPU-hour <input id="r_gpu" type="number" min="0" step="0.01" inputmode="decimal"></label>
<span class="muted">Default prices are AWS on-demand rates for US East (N. Virginia), Linux, read
2026-10-15 from the <a href="{PRICE_PAGE}">EC2 on-demand price page</a> (<a href="{PRICE_FILE}">its
price file</a>). CPU: a c7i.large costs $0.08925 per hour for 2 vCPU, so $0.044625 per CPU-hour.
GPU: a g5.xlarge with one NVIDIA A10G costs $1.006 per hour.</span>
</div>
<figure>
<figcaption><b>Table 1.</b> Tools ranked under the filters above: <span id="filterline"></span>
Fmax is the best balance of precision and recall at any one score cutoff; AUPRC is precision
summed over all cutoffs. Both run from 0 to 1, and higher is better. Median IoU is how well a
correct call's ends match the true region: the overlap divided by the combined span, 1 for a
perfect match. Click a column name to sort; rows with no value always sort last.
<span id="count"></span> On a narrow screen, scroll the table sideways; the tool name stays in
place.</figcaption>
<div class="legend" aria-label="Key to the marks in Table 1">
<span class="item"><span class="prov measured">measured</span> measured on Stanford's Sherlock
cluster, from the run's own task records</span>
<span class="item"><span class="prov reported">reported</span> reported by the tool's submitter, on the
hardware they name (none yet)</span>
<span class="item">— not run: no row in the results table for these filters</span>
<span class="item">— not measured: the search ran, but no task record holds its cost</span>
<span class="item"><i>no correct calls</i>: a median IoU needs at least one correct call on a region
longer than one residue</span>
<span class="item"><span class="tag">needs AlphaFold model</span> <span class="tag">needs GPU</span>
what the tool needs besides the sequence</span>
<span class="item"><span class="tag warn">graded on its own models</span> hmmscan under the Pfam truth
set: Pfam's domains are defined by the models hmmscan runs, so read its row as an upper limit</span>
</div>
{legend_kinds(tools_used)}
<div class="tablebox"><table class="lb" id="lb"><thead><tr></tr></thead><tbody></tbody></table></div>
</figure>
<p class="prose muted">cost_usd_est = core-hours × r_cpu + GPU-hours × r_gpu, per 1,000 queries,
pricing a reserved core-hour as one CPU-hour. Core-hours are cores reserved × wall time, for the
search alone, as the mean over the relatives each tool searched; the run recorded them per tool,
not per relative, so the same number appears for every relative. Building databases is on the
tool cards. No search used a GPU; ProstT5's step that predicts 3Di shape letters used one while
building its databases (see its card). Setting names, decoded: hover a setting, or see the tool's card. Definitions and worked
examples: <a href="methods.html">Methods</a>.</p>
<script type="application/json" id="lb-data">{meta_json}</script>
<script src="assets/app.js"></script>
"""
    return layout(SITE, "leaderboard", body, src,
                  description="Protein annotation tools ranked by Fmax, AUPRC and boundary accuracy, with measured cost.")


def page_cost(meta_json: str, src: dict, n_queries: int, tools_used: set[str]) -> str:
    legend = []
    for kind, k in KINDS.items():
        for t, info in TOOLS.items():
            if info["kind"] != kind or t not in tools_used:
                continue
            if t == "kmerseek":
                for alpha, shape in KMERSEEK_SHAPES.items():
                    legend.append(f'<span class="item">{svg_marker(shape, k["color"])}kmerseek, {esc(alpha)}</span>')
            else:
                legend.append(f'<span class="item">{svg_marker(info["shape"], k["color"])}{esc(info["name"])}</span>')
    body = f"""
<h1>Score against cost</h1>
<p class="lede">The same tools and filters as the leaderboard, drawn as score against compute time.
The line joins the tools that no other tool beats on both: none is both cheaper and better.
Colour is the kind of method; within a kind, marker shape is the tool, and for kmerseek the
alphabet. Peak memory is not on this figure; it is in Table 1.</p>
<form class="filters" onsubmit="return false" aria-label="Filters"><div id="filters"></div>
<p class="axisnote">The benchmark split the proteins along one of the last four rows at a time,
so picking a value in one of them sets the other three back to "any".</p></form>
<div class="metricpick" role="group" aria-label="Score on the vertical axis">
<span>Score on the vertical axis:</span>
<button type="button" class="chip" data-metric="fmax">Fmax</button>
<button type="button" class="chip" data-metric="auprc">AUPRC</button>
<button type="button" class="chip" data-metric="iou">Median IoU</button>
</div>
<figure>
<figcaption><b>Figure 1.</b> Score against core-hours per 1,000 queries for searching one relative
(the mean over the relatives each tool searched), one marker per tool and setting: <span id="filterline"></span> Right is more expensive, up is better. Each step along the
horizontal axis is ten times more CPU time than the last. Hover or focus a marker to see its row;
click it to open the tool's card.</figcaption>
<div class="legend" aria-label="Key to Figure 1">
{"".join(legend)}
<span class="item"><svg width="28" height="10" aria-hidden="true"><line x1="1" y1="5" x2="27" y2="5" stroke="var(--ink)" stroke-width="1.6"/></svg>
no other tool is both cheaper and better</span>
<span class="item">{svg_marker("circle", "var(--ink)")} filled: cost measured on Sherlock</span>
<span class="item">{svg_marker("circle", "var(--ink)", False)} open: cost reported by a submitter (none yet)</span>
</div>
<div class="plotbox" id="plot"></div>
</figure>
<p id="notplotted" class="prose muted"></p>
<script type="application/json" id="lb-data">{meta_json}</script>
<script src="assets/app.js"></script>
"""
    return layout("Score against cost", "cost", body, src,
                  description="Fmax, AUPRC or median IoU against measured CPU-hours per 1,000 queries.")


def fnum(x, d=3):
    return f"{x:.{d}f}"


def best_by(res: pl.DataFrame, **filt) -> pl.DataFrame:
    q = res
    for k, v in filt.items():
        q = q.filter(pl.col(k) == v)
    return q


def card_numbers(tool: str, res: pl.DataFrame, comp: pl.DataFrame) -> dict:
    """Strong and weak strata, cost and the three questions, all read from the tables."""
    species = "all" if tool == "hmmscan" else "zebrafish"
    base = dict(truth_set="pfam", calls="raw", species=species)
    mine_all = best_by(res, tool=tool, stratum_axis="all", **base)
    out: dict = {"species": species, "rows": []}
    if mine_all.height == 0:
        return out
    best = mine_all.sort(["fmax", "setting"], descending=[True, False]).row(0, named=True)
    out["best_setting"] = best["setting"]
    out["all_settings"] = res.filter(pl.col("tool") == tool).select("setting").unique().to_dicts()
    out["best_all"] = best
    others = res.filter(pl.col("tool") != tool)

    comparisons = []
    for axis in ("identity", "plddt", "feature_length_bin"):
        mine = best_by(res, tool=tool, setting=best["setting"], stratum_axis=axis, **base)
        for r in mine.filter(pl.col("n_truth_instances") >= MIN_N).iter_rows(named=True):
            o = best_by(others, stratum_axis=axis, stratum=r["stratum"], **base)
            if o.height == 0:
                continue
            ob = o.sort(["fmax", "tool", "setting"], descending=[True, False, False]).row(0, named=True)
            comparisons.append({
                "axis": axis, "stratum": r["stratum"], "fmax": r["fmax"],
                "n": r["n_truth_instances"], "other_tool": ob["tool"],
                "other_setting": ob["setting"], "other_fmax": ob["fmax"],
                "margin": r["fmax"] - ob["fmax"],
            })
    comparisons.sort(key=lambda c: (-c["margin"], c["axis"], c["stratum"]))
    out["strong"] = comparisons[:2]
    out["weak"] = list(reversed(comparisons[-2:])) if len(comparisons) > 2 else []

    # Swiss-Prot feature types, same relative, this tool's best setting there.
    sp = best_by(res, tool=tool, truth_set="swissprot", calls="raw", species=species,
                 stratum_axis="feature_type")
    if sp.height:
        s_all = best_by(res, tool=tool, truth_set="swissprot", calls="raw", species=species,
                        stratum_axis="all").sort(["fmax", "setting"], descending=[True, False]).row(0, named=True)
        sp = sp.filter((pl.col("setting") == s_all["setting"]) & (pl.col("n_truth_instances") >= MIN_N)) \
            .sort(["fmax", "stratum"], descending=[True, False])
    if sp.height:
        out["sp_setting"] = s_all["setting"]
        out["sp_best"] = sp.row(0, named=True)
        out["sp_worst"] = sp.row(sp.height - 1, named=True)

    c = comp.filter((pl.col("tool") == tool) & (pl.col("setting") == best["setting"]))
    here = c.filter(pl.col("species") == species)
    out["cost_here"] = here.row(0, named=True) if here.height else None


    ident = {r["stratum"]: r for r in best_by(res, tool=tool, setting=best["setting"],
                                              stratum_axis="identity", **base).iter_rows(named=True)}
    out["identity"] = ident
    length = {r["stratum"]: r for r in best_by(res, tool=tool, setting=best["setting"],
                                               stratum_axis="feature_length_bin", **base).iter_rows(named=True)}
    out["length"] = length
    return out


def stratum_words(axis, stratum):
    v = STRATUM_WORDS[axis].get(stratum, stratum)
    if axis == "identity":
        return "no same-family domain in the relative" if stratum == "no_homolog" else f"identity {v} to the closest same-family domain"
    if axis == "feature_length_bin":
        return f"regions {v} long"
    if axis == "plddt":
        return f"AlphaFold confidence (pLDDT) {v}"
    return v


def card_markdown(tool: str, nums: dict, src: dict) -> str:
    t = TOOLS[tool]
    sp_name = SPECIES[nums["species"]]
    where = ("Pfam truth set, raw calls, no relative (hmmscan searches the Pfam models)"
             if nums["species"] == "all" else DEFAULT_WORDS)
    L = [f"# {t['name']}", "",
         f"Tool id `{tool}`. Numbers from {src['repo']} commit `{src['commit'][:8]}`.", "",
         "## What it compares", "", t["compares"], "",
         "## What it needs", ""]
    for k, v in t["needs"].items():
        L.append(f"- {k}: {'yes' if v else 'no'}")
    L += ["", "## Reference searched", "", t["reference"], "",
          "## How to read its score", "", t["score"], ""]

    L += ["## Strong and weak strata", ""]
    if "best_setting" not in nums:
        L += ["Not run under the Pfam truth set for this relative.", ""]
    else:
        b = nums["best_all"]
        L += [f"Under the {where}, its best setting is `{nums['best_setting']}`: Fmax "
              f"{fnum(b['fmax'])} over {b['n_truth_instances']:,} true instances. Below, each "
              "stratum compares this setting with the best setting of any other tool in that "
              "same stratum, so the other tool's setting can change from line to line. "
              f"Only strata with at least {MIN_N} true instances are compared.", ""]
        if nums["strong"]:
            L.append("Where it is closest to, or furthest ahead of, the best other tool:")
            L.append("")
            for c in nums["strong"]:
                L.append(f"- {stratum_words(c['axis'], c['stratum'])}: Fmax {fnum(c['fmax'])} "
                         f"(n = {c['n']:,}); best other: {TOOLS[c['other_tool']]['name']} "
                         f"`{c['other_setting']}`, {fnum(c['other_fmax'])}.")
            L.append("")
        if nums["weak"]:
            L.append("Where it is furthest behind the best other tool:")
            L.append("")
            for c in nums["weak"]:
                L.append(f"- {stratum_words(c['axis'], c['stratum'])}: Fmax {fnum(c['fmax'])} "
                         f"(n = {c['n']:,}); best other: {TOOLS[c['other_tool']]['name']} "
                         f"`{c['other_setting']}`, {fnum(c['other_fmax'])}.")
            L.append("")
        if not nums["strong"] and not nums["weak"]:
            L += ["No other tool was scored in the same strata, so no comparison: not run.", ""]
        if "sp_best" in nums:
            sb, sw = nums["sp_best"], nums["sp_worst"]
            L += [f"Swiss-Prot features, raw calls, {sp_name}, setting `{nums['sp_setting']}`: "
                  f"best feature type {sb['stratum']} ({FEATURES.get(sb['stratum'], '')}), Fmax "
                  f"{fnum(sb['fmax'])} (n = {sb['n_truth_instances']:,}); worst {sw['stratum']} "
                  f"({FEATURES.get(sw['stratum'], '')}), Fmax {fnum(sw['fmax'])} "
                  f"(n = {sw['n_truth_instances']:,}).", ""]
        else:
            L += ["Swiss-Prot feature types: not run.", ""]

    L += ['## Known failure modes {#failure-modes}', ""]
    L += [f"- {f}" for f in t["failure_modes"]]
    L.append("")

    L += ["## Cost per 1,000 queries", ""]
    ch = nums.get("cost_here")
    if ch:
        gpu = ch["gpu_h_per_1000_queries"]
        cost = ch["cpu_h_per_1000_queries"] * R_CPU + gpu * R_GPU
        target = ("scanning Pfam's models" if nums["species"] == "all"
                  else f"searching one relative (mean over the {ch['n_relatives']} relatives it searched)")
        L += [f"Setting `{nums['best_setting']}`, {target}: {ch['cpu_h_per_1000_queries']:.3g} "
              f"core-hours reserved per 1,000 queries ({ch['cores']} cores for the task's wall time), "
              f"measured on Sherlock. Peak memory {ch['peak_rss_gb']:.1f} GB. That is about ${cost:.3g} "
              f"at ${R_CPU:.6g} per CPU-hour."]
        if ch["db_cpu_h"] is not None:
            gtxt = (f", plus {ch['db_gpu_h']:.3g} GPU-hours ({ch['db_gpu_type']})" if ch["db_gpu_h"] else "")
            L += [f"Building its databases took {ch['db_cpu_h']:.3g} core-hours over {ch['db_n']} "
                  f"database{'s' if ch['db_n'] != 1 else ''}{gtxt}. A database is built once and reused, "
                  "so this is not in the per-query cost."]
        else:
            L += ["Building its databases: not measured."]
        L += [""]
    elif "best_setting" in nums:
        L += [f"Setting `{nums['best_setting']}`: not measured (no task record holds its search).", ""]
    else:
        L += ["Not run.", ""]
    L += ["The per-query cost counts the search only. The run recorded cost per tool and setting, "
          "summed over all relatives, so the cost for one relative is the mean.", ""]

    L += ["## Licence", "",
          f"- Licence: {t['licence']} ([source]({t['licence_url']}))",
          f"- Commercial use allowed: {t['commercial']}", ""]

    settings = sorted({r["setting"] for r in nums.get("all_settings", [])})
    if settings:
        L += ["## Settings in the table", ""]
        L += [f"- `{st}`: {setting_words(tool, st)}" for st in settings]
        L += [""]
    L += ["## Version, container, command", "",
          f"- Version: {t['version']}", f"- Container: `{t['container']}`", "",
          "```bash", t["command"], "```", "", t["command_note"], ""]

    L += ["## Will it work for you?", ""]
    L += questions(tool, nums)
    return "\n".join(L) + "\n"


def questions(tool: str, nums: dict) -> list[str]:
    """Three yes/no questions. Each answer is the Fmax range the table shows for it."""
    t = TOOLS[tool]
    q = []

    def answer(rows: dict, yes: list[str], no: list[str], words: dict) -> str | None:
        ok = {k: r for k, r in rows.items() if r["n_truth_instances"] >= MIN_N}
        y = [k for k in yes if k in ok]
        n = [k for k in no if k in ok]
        if not y or not n:
            return None
        show = lambda ks: ", ".join(  # noqa: E731
            f"{words[k]} {fnum(ok[k]['fmax'])} (n = {ok[k]['n_truth_instances']:,})" for k in ks)
        return f"If yes, Fmax here is {show(y)}. If no, {show(n)}."

    base_words = ("Pfam truth set, raw calls, no relative searched" if nums["species"] == "all"
                  else DEFAULT_WORDS)
    where = f"{base_words}, setting `{nums.get('best_setting')}`"
    a1 = answer(nums.get("identity", {}), ["40-60%", "60-100%"], ["0-20%", "20-30%", "30-40%"],
                STRATUM_WORDS["identity"])
    if a1:
        q.append("1. **Is the closest same-family domain in the relative you will search at least "
                 f"40% identical, for most of your proteins?** {a1} ({where}; strata with at least "
                 f"{MIN_N} true instances.)")
    elif nums.get("identity"):
        q.append("1. **Is the closest same-family domain in the relative you will search at least "
                 f"40% identical?** Too few true instances (under {MIN_N}) in the identity strata to say.")
    else:
        q.append("1. **Is every family you expect to find already in Pfam?** If no, this tool "
                 "cannot find it: it only reports families that have a Pfam model.")
    if t["requirement_question"]:
        a, b = t["requirement_question"]
        q.append(f"2. **{a}** {b}")
    else:
        a2 = answer(nums.get("length", {}), ["61-120", "121-250", "251+"],
                    ["2-15", "16-30", "31-60"], STRATUM_WORDS["feature_length_bin"])
        q.append("2. **Are the features you care about longer than 60 residues?** "
                 + (f"{a2} ({where}.)" if a2 else "Fmax by region length for this tool: not run."))
    ch = nums.get("cost_here")
    if ch:
        gpu = ch["gpu_h_per_1000_queries"]
        extra = f" and {gpu:.3g} GPU-hours" if gpu else ""
        if nums["species"] == "all":
            q.append(f"3. **Can you spend about {ch['cpu_h_per_1000_queries']:.3g} core-hours per "
                     "1,000 proteins?** That is the measured cost of scanning them against Pfam's "
                     "models.")
        else:
            q.append(f"3. **Can you spend about {ch['cpu_h_per_1000_queries']:.3g} core-hours{extra} "
                     "per 1,000 proteins for each relative you search?** That is the measured mean "
                     "per relative; see \"Cost per 1,000 queries\" above.")
    else:
        q.append("3. **Can you afford the search?** Its cost was not measured in this run.")
    return q


def page_card(tool: str, md: str, src: dict) -> str:
    body_html = markdown.markdown(md, extensions=["fenced_code", "attr_list", "tables"])
    body = f'<div class="card">{body_html}</div>\n<p><a href="index.html">All tool cards</a> · ' \
           f'<a href="{tool}.md">this card as Markdown</a> · <a href="{tool}.json">as JSON</a></p>'
    return layout(f"{TOOLS[tool]['name']} card", "tools", body, src, depth=1,
                  description=f"What {TOOLS[tool]['name']} compares, what it needs, where it is strong and weak, and what it costs.")


def page_tools_index(src: dict, tools_used: set[str]) -> str:
    items = []
    for kind in KINDS:
        for t, info in TOOLS.items():
            if info["kind"] == kind and t in tools_used:
                items.append(f'<li><span class="kind k-{kind}" aria-hidden="true"></span>'
                             f'<a href="{t}.html"><b>{esc(info["name"])}</b></a><br>'
                             f'<span class="muted">{esc(KINDS[kind]["label"])}</span></li>')
    body = f"""<h1>Tool cards</h1>
<p class="lede">One card per tool: what it compares, what it needs, how to read its score, where it
is strong and weak in the table, what it costs, and the exact command that was run.</p>
{legend_kinds(tools_used)}
<ul class="tools">{"".join(items)}</ul>"""
    return layout("Tool cards", "tools", body, src, depth=1,
                  description="One card per protein annotation tool.")


def parse_tables(md: str) -> list[tuple[list[str], list[list[str]]]]:
    """Every markdown table in coverage.md as (header, rows)."""
    out, cur = [], None
    for ln in md.splitlines():
        if ln.startswith("|"):
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            if cur is None:
                cur = (cells, [])
            elif not set("".join(cells)) <= set("-: "):
                cur[1].append(cells)
        elif cur is not None:
            out.append(cur)
            cur = None
    if cur is not None:
        out.append(cur)
    return out


def parse_coverage(md: str):
    """First table: tool by stratum axis. Returns (axes, rows) where rows are [tool, *counts]."""
    head, rows = parse_tables(md)[0]
    return head[1:], rows


AXIS_TITLES = {"target species": "whole set, per relative", "identity": "% identity to nearest relative",
               "plddt": "AlphaFold pLDDT", "feature_type": "Swiss-Prot feature type",
               "feature_length_bin": "region length"}


def page_coverage(cov_md: str, src: dict, species_mya: dict) -> str:
    tables = parse_tables(cov_md)
    kinds = list(KINDS)

    def tool_cell(t):
        if t in TOOLS:
            return (f'<span class="kind k-{TOOLS[t]["kind"]}" aria-hidden="true"></span>'
                    f'<a href="tools/{t}.html">{esc(TOOLS[t]["name"])}</a>')
        return f'{esc(t)} <span class="muted">(no card yet)</span>'

    def order(t):
        return (kinds.index(TOOLS[t]["kind"]), list(TOOLS).index(t)) if t in TOOLS else (len(kinds), 0)

    figs, todo, n_on, n_off = [], {}, 0, 0
    for k, (head, rows) in enumerate(tables[:2]):
        cols = head[1:]
        labels = [AXIS_TITLES.get(c, c) for c in cols]
        th = "".join(f'<th scope="col">{esc(x)}</th>' for x in labels)
        body = []
        for r in sorted(rows, key=lambda r: order(r[0])):
            tds = []
            for c, v in zip(cols, r[1:]):
                on = v != ""
                n_on += on
                n_off += not on
                # The shuffled-alphabet SCOPe run is a null control for kmerseek, not a gap.
                if not on and not c.endswith("shuffled_k2_seed1"):
                    todo.setdefault(r[0], []).append(AXIS_TITLES.get(c, c) if k == 0 else f"experiment {c}")
                tds.append(f'<td data-measured="{"1" if on else "0"}"><span class="sq{" on" if on else ""}" role="img" '
                           f'aria-label="{"measured" if on else "not run"}"></span>'
                           f'{"<br><span class=mono>" + esc(v) + "</span>" if on else ""}</td>')
            body.append(f'<tr><td class="txt">{tool_cell(r[0])}</td>{"".join(tds)}</tr>')
        title = ("Figure 2a. Tool by stratum axis" if k == 0 else "Figure 2b. Tool by experiment")
        summed = ("experiments, query sets, relatives and truth sets" if k == 0
                  else "the query sets, relatives, truth sets and strata of that experiment")
        figs.append(f"""<figure>
<figcaption><b>{title}.</b> A filled square means the tables hold measured rows for that tool there.
The number under it is the count of true instances scored, summed over {summed}. An empty square
means not run.{" The last column, scope40_v040_shuffled_k2_seed1, is a null control for kmerseek (alphabet letters shuffled), not a run other tools need." if k == 1 else ""}</figcaption>
<div class="tablebox"><table class="cov"><thead><tr><th class="txt" scope="col">Tool</th>{th}</tr></thead>
<tbody>{"".join(body)}</tbody></table></div>
</figure>""")
    todo_html = "".join(f"<li><b>{esc(TOOLS[t]['name'] if t in TOOLS else t)}</b>: {esc('; '.join(v))}.</li>"
                        for t, v in sorted(todo.items(), key=lambda kv: order(kv[0])))
    body = f"""<h1>Coverage</h1>
<p class="lede">What has been measured, across every experiment in the tables, and what has not.
Every empty square is a run still to do.</p>
<div class="legend" aria-label="Key to Figure 2">
<span class="item"><span class="sq on"></span> measured, with the number of true instances scored</span>
<span class="item"><span class="sq"></span> not run</span>
</div>
{legend_kinds(set(TOOLS))}
{"".join(figs)}
<p class="prose muted">{n_on} squares are filled and {n_off} are empty. Where each experiment's output
came from, and which outputs were left out and why: <a href="data/inventory.md">inventory.md</a>
(copied with the tables).</p>
<h2>What to run next</h2>
<ul class="prose">{todo_html}</ul>"""
    return layout("Coverage", "coverage", body, src,
                  description="Which tool, stratum and experiment combinations have been measured and which have not been run.")


def page_methods(res: pl.DataFrame, comp: pl.DataFrame, src: dict, prov: dict) -> str:
    ex = best_by(res, tool="hmmer3_phmmer", setting="default", truth_set="pfam", calls="raw",
                 species="zebrafish", stratum_axis="all", stratum="all").row(0, named=True)
    P, R = ex["fmax_precision"], ex["fmax_recall"]
    f_check = 2 * P * R / (P + R)
    assert abs(f_check - ex["fmax"]) < 1e-9, (f_check, ex["fmax"])
    cx = comp.filter((pl.col("tool") == "hmmer3_phmmer") & (pl.col("species") == "zebrafish")).row(0, named=True)
    cpu_h = cx["cores"] * cx["wall_s"] / 3600
    assert abs(cpu_h - cx["cpu_h_total"]) < 1e-9, (cpu_h, cx["cpu_h_total"])
    per1k = cpu_h / cx["n_relatives"] * 1000 / cx["n_queries"]
    cost = per1k * R_CPU
    ded = best_by(res, tool="hmmer3_phmmer", setting="default", truth_set="pfam", calls="dedup",
                  species="zebrafish", stratum_axis="all", stratum="all").row(0, named=True)
    n_rows = res.height
    pr5 = comp.filter(pl.col("tool") == "prostt5").row(0, named=True)
    n_long = pl.read_parquet(Path(DATA_DIR) / "results.parquet").height
    exps = sorted(pl.read_parquet(Path(DATA_DIR) / "results.parquet")["experiment"].unique().to_list())
    pfn_all = best_by(res, tool="hmmer3_phmmer", truth_set="pfamn", calls="raw", species="zebrafish",
                      stratum_axis="all").row(0, named=True)["n_truth_instances"]
    pfn_nh = best_by(res, tool="hmmer3_phmmer", truth_set="pfamn", calls="raw", species="zebrafish",
                     stratum_axis="identity", stratum="no_homolog").row(0, named=True)["n_truth_instances"]
    pfn = {"all": pfn_all, "nohom": pfn_nh}

    body = rf"""<h1>Methods</h1>
<div class="prose">
<p>The tables hold {n_long:,} rows from {len(exps)} experiments ({", ".join(exps)}); the
coverage page counts all of them. Everything else on this site comes from one of those
experiments, <span class="mono">{EXPERIMENT}</span>, a single benchmark run. It searched {prov['n_queries']:,} human
proteins against the proteomes of 15 relatives from the Quest for Orthologs reference proteomes,
release 2020_04. Of the human proteins, 964 are encoded on chromosome 6; the other 34 are immune
proteins (KIR and LILR receptors, CD1 and MR1, and β2-microglobulin). The run finished on 2026-09-13 on Stanford's Sherlock
cluster. Its tables are <span class="mono">results.parquet</span> ({n_rows:,} rows) and
<span class="mono">compute.parquet</span>, built by
<span class="mono">make leaderboard-tables</span> in {src['repo']} at commit
<a class="mono" href="https://github.com/{src['repo']}/commit/{src['commit']}">{src['commit'][:8]}</a>.
These 998 proteins are not a random sample of the human proteome.</p>

<h2>No tool wins everywhere</h2>
<p>No tool wins everywhere. So this board ranks tools within a stratum, a slice of the proteins
by relative, identity, region length, disorder or feature type, and gives no overall score. Two
earlier community benchmarks found the same. The Quest for Orthologs benchmark ran 15 orthology
methods on 20 benchmarks. Among the balanced methods, "no method consistently obtained a better
balance than the other methods across all data sets". Its authors chose not to compute a combined
score, because the right trade-off depends on the application (Altenhoff et al. 2016). The Massive Text Embedding Benchmark (MTEB), for text embedding models,
reported that "no particular text embedding method dominates across all tasks" (Muennighoff et al.
2023).</p>
<p class="refs">Altenhoff, Adrian M., Brigitte Boeckmann, Salvador Capella-Gutierrez, et al.
"Standardized Benchmarking in the Quest for Orthologs." <i>Nature Methods</i> 13, no. 5 (2016):
425–30. <a href="https://doi.org/10.1038/nmeth.3830">https://doi.org/10.1038/nmeth.3830</a>.<br>
Muennighoff, Niklas, Nouamane Tazi, Loic Magne, and Nils Reimers. "MTEB: Massive Text Embedding
Benchmark." In <i>Proceedings of the 17th Conference of the European Chapter of the Association for
Computational Linguistics</i>, 2014–37. Dubrovnik, Croatia, 2023.
<a href="https://doi.org/10.18653/v1/2023.eacl-main.148">https://doi.org/10.18653/v1/2023.eacl-main.148</a>.</p>

<h2>From a search hit to a call</h2>
<p>Each tool reports aligned regions: a stretch of the human protein matched to a stretch of a
protein in the relative. A region becomes a <b>call</b> when it covers at least half of an annotated
domain or feature on the relative's protein. The call names that annotation and sits on the human
protein where the region is. A call is <b>correct</b> when the human protein carries that
annotation and the call overlaps it with an IoU of at least 0.5. <b>IoU</b> (intersection over
union) is the number of residues the call and the true region share, divided by the number of
residues covered by either:</p>
<p>\[\text{{IoU}} = \frac{{|\text{{call}} \cap \text{{true region}}|}}{{|\text{{call}} \cup \text{{true region}}|}}\]</p>
<p>For example, a call on residues 30–129 and a true domain on residues 10–109 share 80 residues
(30–109). Together they cover 120 (10–129). So IoU = 80 / 120 = 0.67, and the call counts as
correct. A call of the right family that misses the region counts as wrong.</p>
<p>A true annotation is <b>reachable</b> when its family occurs somewhere in the relative's
proteome. A family the relative lacks cannot be copied from it by any tool.</p>

<h2>The scores</h2>
<h3>Fmax</h3>
<p><b>Fmax</b> is the best balance of precision and recall a tool reaches at any score cutoff,
counted protein by protein as in the CAFA function-prediction challenge. At a cutoff \(t\), the
precision of protein \(i\) is its correct calls divided by all its calls scoring at least \(t\). Its
recall is the true annotations found divided by its true annotations. Precision is averaged over
the \(m(t)\) proteins with at least one call at or above \(t\). Recall is averaged over all \(n\)
proteins with an annotation, and counts every true annotation, whether or not the relative could
supply it:</p>
<p>\[\text{{Pr}}(t) = \frac{{1}}{{m(t)}}\sum_{{i}} \frac{{\text{{correct}}_i(t)}}{{\text{{calls}}_i(t)}}, \qquad
\text{{Rc}}(t) = \frac{{1}}{{n}}\sum_{{i}} \frac{{\text{{found}}_i(t)}}{{\text{{true}}_i}}, \qquad
F_{{\max}} = \max_t \frac{{2\,\text{{Pr}}(t)\,\text{{Rc}}(t)}}{{\text{{Pr}}(t)+\text{{Rc}}(t)}}\]</p>
<p>The cutoffs tried are 101 evenly spaced quantiles of the tool's call scores. Worked example:
phmmer searching zebrafish, Pfam truth set, raw calls, all proteins. At the best cutoff (bit score
on the tool's bit score) precision is {P:.4f} and recall is {R:.4f}. So
F = 2 × {P:.4f} × {R:.4f} / ({P:.4f} + {R:.4f}) = {f_check:.4f}, the Fmax in Table 1.</p>

<h3>AUPRC</h3>
<p><b>AUPRC</b>, the area under the precision-recall curve, sums precision over every cutoff,
weighted by the recall gained at that cutoff:</p>
<p>\[\text{{AUPRC}} = \sum_k (R_k - R_{{k-1}})\,P_k\]</p>
<p>Unlike Fmax, AUPRC pools all proteins. \(P_k\) is correct calls ÷ all judged calls at cutoff
\(k\). \(R_k\) is reachable annotations found ÷ all reachable annotations. A tool that makes many calls per protein pays for each wrong one. A call in a
part of a protein that no annotation covers is left out of \(P_k\), since nobody has looked there;
Fmax counts it. Worked
example, same row: AUPRC = {ex['auprc']:.4f}. The curve itself is in the source run
(<span class="mono">all_domain_curves.parquet</span>), not in these tables, so the sum cannot be
recomputed here.</p>

<h3>Median IoU</h3>
<p><b>Median IoU</b> is the middle IoU of all correct calls, with single-residue annotations left
out (an IoU against one residue says nothing about where the ends are). Worked example, same row:
{ex['median_iou_tp']:.3f}, so half of phmmer's correct calls on zebrafish overlap their true domain
by at least that much. A correct call on a region longer than one residue has an IoU of at least
0.5. So the run stores 0 only when there is no such call. Either there are no correct calls at
all, or the correct calls are on single residues only (the 1 aa length bin and the ACT_SITE,
BINDING and SITE features).
The pages show those as <i>no correct calls</i>, not as 0.</p>

<h3>True instances</h3>
<p><b>True instances</b> (n) is the number of annotations on the human proteins in the stratum,
reachable or not. Same row: {ex['n_truth_instances']:,}.</p>

<h3>Core-hours per 1,000 queries</h3>
<p><b>Core-hours</b> are taken from the run's Nextflow task records: the cores a task reserved
multiplied by its wall time, so a task that sat partly idle still counts all its cores. The tables
record them per tool and setting, summed over every relative searched, with building the databases
as separate rows. The leaderboard divides the search's core-hours by the number of relatives
searched and by the {cx['n_queries']} queries, times 1,000: the mean cost of searching one relative.
Worked example: phmmer reserved {cx['cores']} cores for {cx['wall_s']:,.0f} s over its
{cx['n_relatives']} searches. That is {cx['cores']} × {cx['wall_s']:,.0f} / 3,600 = {cpu_h:.4f}
core-hours, and {cpu_h:.4f} / {cx['n_relatives']} × 1,000 / {cx['n_queries']} = {per1k:.4f}
core-hours per 1,000 queries per relative. Peak RAM is the largest resident memory of any of
those tasks. The tool cards give the database builds separately.</p>
<p>Numbers marked <span class="prov measured">measured</span> come from these records, on Sherlock
(partition hns). No search used a GPU. ProstT5 has to turn sequences into 3Di shape letters
before it can search; that step ran on one GPU as part of building its databases, queries
included, and is on the ProstT5 card ({pr5['db_gpu_h']:.3g} GPU-hours over {pr5['db_n']}
databases), not in the per-query cost. Numbers marked
<span class="prov reported">reported</span> would come from a tool's submitter, on hardware they
name; none have been submitted yet.</p>

<h3>Estimated cost</h3>
<p>\[\text{{cost\_usd\_est}} = \text{{CPU-hours}} \times r_\text{{cpu}} + \text{{GPU-hours}} \times r_\text{{gpu}}\]</p>
<p>per 1,000 queries. The default \(r_\text{{cpu}}\) = ${R_CPU:.6g} per CPU-hour is half of AWS's
on-demand price for a c7i.large (2 vCPU, $0.08925 per hour). The default \(r_\text{{gpu}}\) =
${R_GPU} per GPU-hour is a g5.xlarge (one NVIDIA A10G). Both are for US East (N. Virginia), Linux,
from the <a href="{PRICE_PAGE}">EC2 on-demand price page</a> (<a href="{PRICE_FILE}">price
file</a>, published 2026-10-12, read 2026-10-15). A core-hour is priced as one CPU-hour. Memory is not
priced. Worked example: {per1k:.4f} × {R_CPU:.6g} + 0 ×
{R_GPU} = ${cost:.5f} for phmmer against zebrafish. Change either price on the leaderboard page.</p>

<h2>Truth sets</h2>
<h3>Pfam</h3>
<p>Every Pfam-A domain with coordinates on the human proteins, from the Pfam regions file for
UniProt reference proteomes. The same file supplies the domains on each relative's proteins, which
are what a call copies. Pfam families are defined by profile models, so tools that build profiles
(jackhmmer, HH-suite) are graded against the kind of model they build, and hmmscan against its own
models. Release of the regions file: not recorded in the run; the hmmscan models are Pfam-A 38.2.
hmmscan was scored on raw calls only, so it shows "not run" under dedup.</p>
<h3>Pfam-N</h3>
<p>Pfam-N is an extension of Pfam made by a deep-learning model: it labels regions the Pfam
profile models missed. It is only published for Pfam 35.0. This is the one truth set the profile
models did not produce. Its regions were filtered to the human proteins queried here. Few of them
can be copied from any relative: against zebrafish, {pfn['nohom']:,} of the {pfn['all']:,} Pfam-N
regions have no domain of the same family in the zebrafish proteome. So every tool scores low on
Pfam-N, and the set says little about which tool is better.</p>
<h3>Swiss-Prot features</h3>
<p>Curated features from Swiss-Prot entries, of 12 types (the feature keys on the leaderboard,
such as TRANSMEM for a membrane-spanning segment and ACT_SITE for an active site residue). They
come from the literature, not from any sequence or structure model. Single-residue features
(ACT_SITE, BINDING, SITE) are scored by whether the call covers the residue, since an IoU against
one residue cannot reach 0.5. A Swiss-Prot feature counts as reachable when the relative has a
protein that carries the same feature type and shares a Pfam family with the human protein.</p>

<h2>Dedup against raw</h2>
<p><b>Raw</b> keeps every call a tool makes. One human region that matches 300 mouse kinases
becomes 300 kinase calls on nearly the same residues, and all but one count as wrong.
<b>Dedup</b> merges calls of the same family on the same human protein that overlap each other with
an IoU of 0.5 or more, keeping the highest-scoring one. Two copies of a domain side by side
barely overlap, so both are kept. Comparing the two separates finding a domain from reporting it many times.
For phmmer on zebrafish (Pfam, all proteins): Fmax {ex['fmax']:.3f} raw and {ded['fmax']:.3f} dedup.</p>

<h2>Strata</h2>
<ul>
<li><b>Nearest reviewed relative, % identity</b>: for each human domain, the percent identity to the
closest domain of the same family in the relative. Bins: 0–20, 20–30, 30–40, 40–60 and 60–100%. A
domain whose family is absent from the relative is "no same-family domain".</li>
<li><b>Nearest reviewed relative, Mya</b>: the relative searched, labelled with how many million
years ago its line split from ours.</li>
<li><b>Region type and length</b>: the length of the annotated region, in residues. Bins up to 60
aa are short features; 61 aa and up is the size of a folded domain. Whole-protein annotations
were not scored in this run.</li>
<li><b>Disorder</b>: the mean AlphaFold confidence (pLDDT, 0–100) of the human protein. Low pLDDT
marks regions without a fixed shape.</li>
<li><b>Feature type</b>: the Swiss-Prot feature key, for the Swiss-Prot truth set only.</li>
</ul>
<p>The run computed each stratum on its own, so two strata cannot be combined; a combination with
no row shows "—" and "not run".</p>

<h2>Checks on this site</h2>
<p>A headless browser renders every filter combination and confirms that each "—" on the
leaderboard is a missing row in results.parquet or compute.parquet. The source table checksums are
in <span class="mono">data/SOURCE.json</span> of the site repository.</p>
</div>
<link rel="stylesheet" href="assets/katex/katex.min.css">
<script defer src="assets/katex/katex.min.js"></script>
<script defer src="assets/katex/contrib/auto-render.min.js"
 onload="renderMathInElement(document.body,{{delimiters:[{{left:'\\[',right:'\\]',display:true}},{{left:'\\(',right:'\\)',display:false}}]}})"></script>
"""
    return layout("Methods", "methods", body, src,
                  description="How each score is defined and computed, how the truth sets were built, and what dedup means.")


def page_pileup(src: dict, res: pl.DataFrame, prov: dict, data: dict, check: dict) -> str:
    regions = data["regions"]
    import statistics
    med_len = statistics.median(r["len"] for r in regions)
    n_all = len(regions)
    n_past = sum(r["past_end"] for r in regions)
    n_outside = sum(r["c0"] is None or r["c0"] > 199 for r in regions)
    n_in = n_all - sum(1 for r in regions if r["s"] > r["plen"])
    n_fam = len({r["fam"] for r in regions})
    past = "".join(f'<li>{r["acc"]} {r["fam"]}: {r["s"]}–{r["e"]} on a {r["plen"]} aa sequence</li>'
                   for r in regions if r["past_end"])
    km = best_by(res, tool="kmerseek", truth_set="pfam", calls="raw", species="zebrafish",
                 stratum_axis="all").sort(["fmax", "setting"], descending=[True, False]).row(0, named=True)
    dashes = {"hmmer3_phmmer": "", "mmseqs2_seqseq": "6 3", "hmmer3_jackhmmer": "", "mmseqs2_iterative": "6 3",
              "hhblits": "2 2", "foldseek": "", "prostt5": "6 3", "kmerseek": "", "hmmscan": ""}
    hist = [{"tool": t, "setting": st, "dash": dashes[t]} for t, st in
            [("hmmer3_phmmer", "default"), ("mmseqs2_seqseq", "s7"), ("hmmer3_jackhmmer", "n3"),
             ("mmseqs2_iterative", "s7"), ("hhblits", "single_seq"), ("foldseek", "3di_aa"),
             ("prostt5", "3di_from_seq"), ("kmerseek", km["setting"]), ("hmmscan", "direct")]]
    # Figure 4's numbers for its caption, from the same per-region offsets the page draws.
    early = {}
    for h in hist:
        pnl = next(q for q in data["panels"] if q["tool"] == h["tool"] and q["setting"] == h["setting"])
        ds = [row["d"] for row in pnl["rows"] if row["d"] is not None]
        h["n"] = len(ds)
        early[h["tool"]] = sum(d < -100 for d in ds) / len(ds)
    aligners = ["hmmer3_phmmer", "mmseqs2_seqseq", "hmmer3_jackhmmer", "mmseqs2_iterative", "hhblits", "foldseek", "prostt5"]
    lo_t = min(aligners, key=lambda t: early[t]); hi_t = max(aligners, key=lambda t: early[t])
    scores = {t: ("region_enrichment (shared k-mers over the number expected by chance)" if t == "kmerseek" else "bit score")
              for t in TOOLS}
    words = {f'{q["tool"]}|{q["setting"]}': setting_words(q["tool"], q["setting"]) for q in data["panels"]}
    meta = {
        "kinds": {k: v["color"] for k, v in KINDS.items()},
        "names": {t: TOOLS[t]["name"] for t in TOOLS}, "words": words, "scores": scores,
        "median_len": med_len, "focus_default": "hmmer3_phmmer", "hist": hist,
    }
    hist_legend = "".join(
        f'<span class="item"><svg width="30" height="10" aria-hidden="true"><line x1="1" y1="5" x2="29" y2="5" '
        f'stroke="{KINDS[TOOLS[h["tool"]]["kind"]]["color"]}" stroke-width="2" stroke-dasharray="{h["dash"]}"/></svg>'
        f'{esc(TOOLS[h["tool"]]["name"])}{" " + esc(h["setting"]) if h["tool"] == "kmerseek" else ""} '
        f'(n = {h["n"]:,})</span>' for h in hist)
    body = f"""
<h1>Region pile-up</h1>
<p class="source">The pile-up's own tables (regions and calls) come from {src['repo']} at commit
<a class="mono" href="https://github.com/{src['repo']}/commit/{src.get('pileup_commit') or src['commit']}">{(src.get('pileup_commit') or src['commit'])[:8]}</a>
(<span class="mono">results/leaderboard/pileup</span>).</p>
<p class="lede">Where each tool's calls fall around a known domain. A <b>call</b> is a stretch of a
human protein that a tool matched to a zebrafish protein carrying a Pfam domain. The call names that
domain and spans the whole matched stretch, which can be much wider than the domain itself.
Every row below is one <b>region</b>: one of the {n_all:,} Pfam domains on the 998 human query
proteins. Every panel is one tool and setting searching the zebrafish proteome (zebrafish and
human last shared an ancestor about 430 million years ago), and every panel has the same rows.</p>
<p class="prose">A residue counts as <b>covered</b> when a call naming that region's own Pfam family
spans it, at any score the tool reported (E-value cutoff 10). Covered does not mean correctly
placed: a call that runs far past the domain still covers it. Each row shows 100 residues before
the region, then the region itself stretched or squeezed into 100 equal steps, then 100 residues
after it, as deepTools does for genes. Each step shows the residue at its centre; the median region
is {med_len:.0f} aa, so one step is about {med_len / 100:.1f} aa.</p>
<div class="metricpick" role="group" aria-label="Row order">
<span>Order rows by:</span>
<button type="button" class="chip" data-sort="coverage">coverage in the chosen panel, most first (missed regions last)</button>
<button type="button" class="chip" data-sort="identity">% identity to the closest zebrafish domain of the same family, highest first</button>
<button type="button" class="chip" data-sort="length">region length, longest first</button>
<button type="button" class="chip" data-sort="plddt">mean pLDDT over the region, highest first</button>
</div>
<p class="axisnote">The order is set by one panel and used in all of them, so the rows line up across
panels. In the other panels the same order can look like scattered stripes. Identity is the best
pairwise identity between the human domain and any zebrafish domain of its family. Regions with no
such domain, or no matching AlphaFold model, sort last.</p>
<p id="pu-loading" class="muted">Loading 23 panels…</p>
<figure>
<figcaption><b>Figure 3.</b> Calls piled up around {n_all:,} Pfam regions, one panel per tool and
setting. At the top of each panel is the share of regions covered at each position (the column
average). A tool that covers the region and nothing around it draws a high step between the start
and end lines. Below it, every region as a row, squeezed to fit; darker means more of the rows in
that band are covered. Click a panel to draw it full size below, one pixel row per region.</figcaption>
{legend_kinds({q["tool"] for q in data["panels"]})}
<div class="legend" aria-label="Key to Figure 3">
<span class="item"><span class="sw" style="background:var(--surface);border:1px solid var(--hair)"></span> not covered</span>
<span class="item"><span class="sw" style="background:var(--null)"></span> past the protein's end (no residue there)</span>
<span class="item"><svg width="28" height="10" aria-hidden="true"><line x1="1" y1="5" x2="27" y2="5" stroke="var(--k-letters)" stroke-width="2"/></svg>
coloured curve: share of all regions covered</span>
<span class="item"><svg width="28" height="10" aria-hidden="true"><line x1="1" y1="5" x2="27" y2="5" stroke="var(--ink)" stroke-width="1.6"/></svg>
black curve: the same with one region per Pfam family, since regions from one protein or family are not independent</span>
</div>
<div class="prose muted pu-readkey">
<p>The numbers under each panel: <b>missed</b> is the share of regions with no covered residue
inside them (lower is better). <b>Inside</b> is the mean coverage inside the region (higher is
better). <b>Flanks</b> is the mean coverage in the 100 residues either side. For placing domain
boundaries lower is better, though a second copy of the same domain nearby counts as covered too.
<b>0.8→0.2</b> is how many residues the curve takes to fall from 0.8 to 0.2 at each end of the region:
fewer residues means a sharper edge. "Not reached" means the curve never gets to 0.8 inside the
region; "does not fall that far" means it stays above 0.2 for all 100 residues outside.</p>
<p>Counts: {n_in:,} of the {n_all:,} regions have at least one residue inside their protein's
sequence; the readouts use those. With one region per Pfam family ({n_fam:,} families), the region
kept is the first by protein accession, then start. {n_past} Pfam regions run partly or wholly past the
end of the protein sequence this run searched, so their Pfam coordinates belong to a different
version of the sequence.</p>
<details><summary>The {n_past} regions that run past their protein's end</summary><ul class="mono">{past}</ul></details>
</div>
<div id="pu-grid" class="pu-grid"></div>
</figure>
<figure>
<figcaption><b>Figure 3, chosen panel:</b> <span id="pu-ftitle"></span>. <span id="pu-caption-numbers"></span>
Hover a row to see its protein, domain, coordinates, calls and residues.</figcaption>
<div class="legend" aria-label="Key to the chosen panel">
<span class="item"><span class="sw" style="background:var(--strip-a)"></span><span class="sw" style="background:var(--strip-b)"></span>
left strip: rows from one protein share a shade; the shade changes where the protein changes</span>
<span class="item">Colour cells by: <button type="button" class="chip" data-color="coverage">covered or not</button>
<button type="button" class="chip" data-color="score">call score</button></span>
<span class="item" id="pu-scorebar" hidden><span class="name"></span><span class="lo mono"></span><span class="bar"></span><span class="hi mono"></span>
highest score of the calls over the residue, 5th to 95th percentile of this tool's scores</span>
</div>
<div class="plotbox pu-focus">
<svg id="pu-fcurve" class="pu-fcurve" role="img" aria-label="share of regions covered at each position"></svg>
<div class="pu-scroll"><canvas id="pu-fheat" class="pu-fheat" aria-label="one row per region"></canvas>
<div id="pu-tip" class="tip pu-tip" hidden></div></div>
</div>
<p id="pu-fro" class="prose muted"></p>
</figure>
<h2>Where calls start</h2>
<figure>
<figcaption><b>Figure 4.</b> For each region with a call, take the call that overlaps it most (ties go to
the higher score, then the earlier start). Subtract the region's first residue from the call's first
residue. 0 means the call starts on the domain's first residue; negative values start early. One line
per tool, as a share of that tool's regions with a call (n in the key), in 5-residue steps; the end
bins hold everything beyond 100 residues. For kmerseek the line is its setting with the highest Fmax
(the best balance of precision and recall, see <a href="methods.html">Methods</a>) against
zebrafish, <span class="mono">{esc(km["setting"])}</span>. Among the aligners, between
{100 * early[lo_t]:.0f}% ({esc(TOOLS[lo_t]["name"])}) and {100 * early[hi_t]:.0f}%
({esc(TOOLS[hi_t]["name"])}) of calls start more than 100 residues early. That is because a call
spans the whole matched stretch. hmmscan's calls start on the domain because Pfam's domains are defined by the
models hmmscan runs. Where calls end is not drawn.</figcaption>
<div class="legend" aria-label="Key to Figure 4">{hist_legend}</div>
<div class="plotbox"><svg id="pu-hist" role="img" aria-label="histogram of call start minus true start"></svg></div>
</figure>
<script type="application/json" id="pu-meta">{json.dumps(meta)}</script>
<script src="assets/pileup.js"></script>
"""
    return layout("Region pile-up", "pileup", body, src,
                  description="Where each tool's calls fall around known Pfam domains, residue by residue.")
from html.parser import HTMLParser


class _Text(HTMLParser):
    """Headings, paragraphs, captions, list items and table cells of <main>, in order."""
    KEEP = {"h1", "h2", "h3", "p", "figcaption", "li", "dt", "dd", "summary", "pre"}

    def __init__(self):
        super().__init__()
        self.out, self.stack, self.in_main, self.skip = [], [], False, 0

    def handle_starttag(self, tag, attrs):
        if tag == "main":
            self.in_main = True
        if tag in ("script", "style", "svg"):
            self.skip += 1
        if self.in_main and tag in self.KEEP:
            self.stack.append([tag, ""])

    def handle_endtag(self, tag):
        if tag in ("script", "style", "svg"):
            self.skip -= 1
        if tag == "main":
            self.in_main = False
        if self.stack and self.stack[-1][0] == tag:
            t, text = self.stack.pop()
            text = " ".join(text.split())
            if text:
                self.out.append({"tag": t, "text": text})

    def handle_data(self, data):
        if self.stack and not self.skip:
            self.stack[-1][1] += data


def page_json(out: Path, rel: str, src: dict, data) -> None:
    """Write <page>.json beside <page>: the page's text in reading order plus its data, so an
    agent gets what a person gets."""
    html_text = (out / rel).read_text()
    parser = _Text()
    parser.feed(html_text)
    title = re.search(r"<title>(.*?)</title>", html_text).group(1)
    doc = {"page": rel, "title": html.unescape(title), "source": src, "text": parser.out, "data": data}
    (out / (rel + ".json")).write_text(json.dumps(doc, indent=1, ensure_ascii=False, allow_nan=False) + "\n")
    name = rel.split("/")[-1]
    link = f'<p>This page as JSON, for scripts and agents: <a class="mono" href="{name}.json">{name}.json</a></p>\n'
    (out / rel).write_text(html_text.replace("</div></footer>", link + "</div></footer>", 1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data"))
    ap.add_argument("--out", type=Path, default=Path("docs"))
    a = ap.parse_args()
    global DATA_DIR
    DATA_DIR = a.data
    res, comp, src, prov, cov = load(a.data)

    unknown = set(res["tool"].unique()) - set(TOOLS)
    if unknown:
        raise SystemExit(f"tools with no card facts: {unknown}")
    species_mya = {r["species"]: r["species_mya"] for r in
                   res.select("species", "species_mya").unique().iter_rows(named=True)
                   if r["species_mya"] is not None and not math.isnan(r["species_mya"])}
    tools_used = set(res["tool"].unique())
    rows = row_order(res)
    searched = sorted({f'{r["tool"]}|{r["setting"]}|{r["species"]}'
                       for r in res.select("tool", "setting", "species").unique().iter_rows(named=True)})
    meta = {
        "defaults": DEFAULTS,
        "chiprows": mark_empty(chiprows(species_mya), res),
        "rows": rows,
        "tools": {t: {"name": TOOLS[t]["name"], "kind": TOOLS[t]["kind"], "shape": TOOLS[t]["shape"],
                      "tags": [w for k, w in (("AlphaFold model", "needs AlphaFold model"), ("GPU", "needs GPU"))
                               if TOOLS[t]["needs"][k]]}
                  for t in tools_used},
        "kinds": KINDS,
        "species": {s: {"name": SPECIES[s], "mya": int(species_mya[s]) if s in species_mya else None}
                    for s in SPECIES},
        "truthNames": TRUTH,
        "metrics": METRICS,
        "searched": searched,
        "source": src,
    }
    payload = json.dumps({"meta": meta, "results": results_json(res), "compute": compute_json(comp)},
                         separators=(",", ":"), allow_nan=False)
    payload = payload.replace("</", "<\\/")
    n_q = prov["n_queries"]

    out = a.out
    (out / "tools").mkdir(parents=True, exist_ok=True)
    (out / "data").mkdir(parents=True, exist_ok=True)
    (out / ".nojekyll").write_text("")
    phm = res.filter((pl.col("tool") == "hmmer3_phmmer") & (pl.col("truth_set") == "pfam")
                     & (pl.col("species") == "zebrafish") & (pl.col("stratum_axis") == "all"))
    dd = {f"{k}_{m}": phm.filter(pl.col("calls") == k)[m][0] for k in ("raw", "dedup") for m in ("fmax", "auprc")}
    (out / "index.html").write_text(page_leaderboard(payload, src, n_q, tools_used, dd))
    (out / "cost.html").write_text(page_cost(payload, src, n_q, tools_used))
    (out / "coverage.html").write_text(page_coverage(cov, src, species_mya))
    (out / "methods.html").write_text(page_methods(res, comp, src, prov))
    (out / "tools" / "index.html").write_text(page_tools_index(src, tools_used))
    for t in sorted(tools_used):
        nums = card_numbers(t, res, comp)
        md = card_markdown(t, nums, src)
        (out / "tools" / f"{t}.md").write_text(md)
        card = {
            "tool": t, "name": TOOLS[t]["name"], "kind": KINDS[TOOLS[t]["kind"]]["label"],
            "compares": TOOLS[t]["compares"], "needs": TOOLS[t]["needs"],
            "reference": TOOLS[t]["reference"], "score": TOOLS[t]["score"],
            "failure_modes": TOOLS[t]["failure_modes"], "version": TOOLS[t]["version"],
            "licence": TOOLS[t]["licence"], "commercial_use": TOOLS[t]["commercial"],
            "licence_url": TOOLS[t]["licence_url"],
            "container": TOOLS[t]["container"], "command": TOOLS[t]["command"],
            "numbers_from": {"repo": src["repo"], "commit": src["commit"],
                             "filters": "pfam, raw, " + nums["species"]},
            "best_setting": nums.get("best_setting"),
            "strong": nums.get("strong", []), "weak": nums.get("weak", []),
            "cost_per_1000_queries": None if not nums.get("cost_here") else {
                "core_h_reserved_per_relative": nums["cost_here"]["cpu_h_per_1000_queries"],
                "n_relatives": nums["cost_here"]["n_relatives"],
                "database_build_core_h": nums["cost_here"]["db_cpu_h"],
                "database_build_gpu_h": nums["cost_here"]["db_gpu_h"],
                "peak_rss_gb": nums["cost_here"]["peak_rss_gb"],
                "hardware": nums["cost_here"]["hardware"]},
            "questions": [re.sub(r"\*\*", "", q[3:]) for q in questions(t, nums)],
            "markdown": f"{t}.md",
        }
        (out / "tools" / f"{t}.json").write_text(json.dumps(card, indent=1, allow_nan=False) + "\n")
        (out / "tools" / f"{t}.html").write_text(page_card(t, md, src))
    (out / "data" / "SOURCE.json").write_text(json.dumps(src, indent=1) + "\n")
    if (a.data / "pileup/regions.parquet").exists():
        import pileup
        calls = pl.read_parquet(a.data / "pileup/calls.parquet").select("tool", "setting").unique()
        kinds = list(KINDS)
        panels = sorted(({"tool": t, "setting": st, "kind": TOOLS[t]["kind"]} for t, st in calls.iter_rows()),
                        key=lambda q: (kinds.index(q["kind"]), list(TOOLS).index(q["tool"]), q["setting"]))
        data, check = pileup.build(a.data, panels)
        (out / "data" / "pileup.json").write_text(json.dumps(data, separators=(",", ":"), allow_nan=False))
        (out / "data" / "pileup_check.json").write_text(json.dumps(check, indent=1) + "\n")
        prov["pileup"] = json.loads((a.data / "pileup/provenance.json").read_text())
        (out / "pileup.html").write_text(page_pileup(src, res, prov, data, check))

    # Every page also as JSON at its path plus .json.
    full = json.loads(payload.replace("<\\/", "</"))
    page_json(out, "index.html", src, {"note": "The table is drawn from these rows: results keyed "
              "truth|calls|species|stratum_axis|stratum, then tool|setting -> [fmax, auprc, median_iou_tp, "
              "n_truth_instances, n_tp_calls]; compute keyed tool|setting|species.", **full,
              "prices": {"r_cpu": R_CPU, "r_gpu": R_GPU, "source": PRICE_PAGE, "price_file": PRICE_FILE}})
    page_json(out, "cost.html", src, {"note": "Same rows as index.html.json; the plot is score against "
              "cpu_h_per_1000 for every tool and setting that has both.", "results": full["results"],
              "compute": full["compute"], "rows": full["meta"]["rows"]})
    page_json(out, "coverage.html", src, {"tables": [
        {"columns": h[1:], "rows": [{"tool": r[0], "instances": {c: (int(v.replace(",", "")) if v else None)
                                                             for c, v in zip(h[1:], r[1:])}} for r in rows]}
        for h, rows in parse_tables(cov)[:2]]})
    if (a.data / "inventory.md").exists():
        (out / "data" / "inventory.md").write_text((a.data / "inventory.md").read_text())
    page_json(out, "methods.html", src, None)
    page_json(out, "tools/index.html", src, {"cards": [f"{t}.json" for t in sorted(tools_used)]})
    if (out / "pileup.html").exists():
        page_json(out, "pileup.html", src, {
            "note": "Readouts computed by scripts/pileup.py readouts(), the same rules as the page's "
                    "readouts(); the matrix itself is data/pileup.json.",
            "matrix": "data/pileup.json", "per_residue_recall_check": check, **pileup.readouts(data)})
    for t in sorted(tools_used):
        page_json(out, f"tools/{t}.html", src, json.loads((out / "tools" / f"{t}.json").read_text()))
    print(f"built {out}: {len(rows)} tool-setting rows, {len(tools_used)} cards, payload {len(payload) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
