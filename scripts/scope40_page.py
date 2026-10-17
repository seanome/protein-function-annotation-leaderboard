"""The SCOPe40 ranking page (scope40.html) and its data file (data/scope40.json).

Reads data/scope40_long.parquet and data/scope40_summary.json (scripts/scope40_score.py)
and data/scope40_testbed.json (scripts/scope40_testbed.py). Every number in the page text is
read from those files at build time; the figures draw from data/scope40.json in the browser.
"""

import html
import json
from pathlib import Path

import polars as pl

import scope40_rules as R
from tools import KINDS, TOOLS

FIG_CURVES, FIG_CLASSES, TABLE = 19, 20, 3
LEVELS = ["family", "superfamily", "fold"]
CLASS_NAMES = {"a": "all alpha", "b": "all beta", "c": "alpha/beta", "d": "alpha+beta",
               "e": "multi-domain proteins", "f": "membrane and cell surface", "g": "small proteins"}
# The analysis-repository PR that holds the pipeline for the tools not run yet.
PIPELINE_PR = "https://github.com/seanome/2024-kmerseek-analysis/pull/166"

SCORED_NOTE = "from the per-query file other authors published; our own run is in the pipeline"
NOT_RUN = "no per-query hit list on SCOPe40 yet; pipeline written, not launched"
SCORED_NOTES = {"kmerseek": "from the pair tables the PR 107 combination runs kept on Sherlock, ranked by "
                            "region E-value; scored with Foldseek's rules"}
# Rows of Table 3 in the order of Table 1, then the two no-search rows.
ROWS = [
    ("kmerseek", NOT_RUN + "; the PR 107 runs kept one area per setting, not the hits"),
    ("hmmer3_phmmer", NOT_RUN),
    ("hmmer3_jackhmmer", NOT_RUN),
    ("mmseqs2_seqseq", SCORED_NOTE),
    ("mmseqs2_iterative", NOT_RUN),
    ("hhblits", NOT_RUN),
    ("foldseek", SCORED_NOTE),
    ("prostt5", NOT_RUN),
    ("hmmscan", "it scores a sequence against Pfam family models, not one domain against another"),
]
NOSEARCH = [
    ("kyte_doolittle", "Kyte-Doolittle composition scan (notebook 231)",
     "not run: not yet run on SCOPe40; it would rank targets by its own score"),
    ("shuffled", "Dipeptide-preserving shuffle of every query (chance level)",
     "not run: needs every tool run again on the shuffled queries, same settings"),
]


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def pct(x: float) -> str:
    return f"{100 * x:.1f}%"


def best_setting(summ: dict, tool: str) -> str | None:
    """The tool's setting with the highest family AUC (the site's best-setting convention)."""
    rows = [r for r in summ["overall"] if r["tool"] == tool and r["level"] == "family"]
    return max(rows, key=lambda r: (r["auc"], r["setting"]))["setting"] if rows else None


def row_info(tool: str, setting: str, summ: dict) -> dict:
    """Name, colour and marks for one (tool, setting) row on the page."""
    if tool.endswith("_shuffled"):
        base = tool[: -len("_shuffled")]
        return {"name": f"{TOOLS[base]['name']} {setting}, queries shuffled (chance level)", "kind": "nosearch",
                "shape": "circle", "color": "var(--k-nosearch)", "dash": "5 3",
                "primary": setting == best_setting(summ, base), "chance": True}
    best = best_setting(summ, tool)
    name = TOOLS[tool]["name"] + (f" {setting}" if tool == "kmerseek" else "")
    return {"name": name, "kind": TOOLS[tool]["kind"], "shape": TOOLS[tool]["shape"],
            "color": f"var(--k-{TOOLS[tool]['kind']})", "dash": "", "primary": setting == best, "chance": False}


def data_json(data_dir: Path) -> dict:
    long = pl.read_parquet(data_dir / "scope40_long.parquet")
    summ = json.loads((data_dir / "scope40_summary.json").read_text())
    tools = long.select("tool", "setting").unique(maintain_order=True).rows()
    q = (long.filter((pl.col("level") == "family") & (pl.col("tool") == tools[0][0])
                     & (pl.col("setting") == tools[0][1]))
         .sort("query").select("query", "scope_sccs", "scope_class"))
    out = {"queries": q["query"].to_list(), "sccs": q["scope_sccs"].to_list(), "cls": q["scope_class"].to_list(),
           "tools": [], "sens": {}, "hits": {}, "summary": summ, "levels": LEVELS,
           "classNames": {c: n for c, n in CLASS_NAMES.items() if c in set(q["scope_class"].to_list())}}
    for t, s in sorted(tools, key=lambda x: (x[0].endswith("_shuffled"), x[0], x[1])):
        key = f"{t}|{s}"
        out["tools"].append({"key": key, "tool": t, "setting": s, **row_info(t, s, summ)})
        g = long.filter((pl.col("tool") == t) & (pl.col("setting") == s))
        out["sens"][key] = {}
        for lvl in LEVELS:
            gl = g.filter(pl.col("level") == lvl).sort("query")
            assert gl["query"].to_list() == out["queries"], (t, s, lvl)
            out["sens"][key][lvl] = [round(v, 6) for v in gl["sensitivity"].to_list()]
        out["hits"][key] = [int(x) for x in g.filter(pl.col("level") == "family").sort("query")["has_hits"].to_list()]
    return out


def auc_cells(by: dict, tool: str, setting: str) -> tuple[list[str], float | None]:
    cells, cov = [], None
    for lvl in LEVELS:
        r = by.get((tool, setting, lvl))
        if r:
            cov = r["coverage"]
            cells.append(f'<td class="num" data-v="{r["auc"]:.6f}">{r["auc"]:.3f}'
                         f'<span class="key" style="display:block">{r["lo"]:.3f}–{r["hi"]:.3f}</span></td>')
        else:
            cells.append('<td class="num" data-v="">—</td>')
    return cells, cov


def table_rows(summ: dict) -> str:
    by = {(r["tool"], r["setting"], r["level"]): r for r in summ["overall"]}
    settings = {}
    for r in summ["overall"]:
        settings.setdefault(r["tool"], set()).add(r["setting"])
    out = []

    def row(tool_id, name, kind, setting, note, status, cls=""):
        cells, cov = auc_cells(by, tool_id, setting)
        measured = cov is not None
        return (f'<tr data-tool="{tool_id}" data-setting="{esc(setting)}"{f" class={chr(34)}{cls}{chr(34)}" if cls else ""}>'
                f'<td class="txt"><span class="kind k-{kind}" aria-hidden="true"></span> {esc(name)}</td>'
                f'<td class="txt">{esc(setting if measured else "—")}</td>' + "".join(cells)
                + (f'<td class="num" data-v="{cov:.6f}">{pct(cov)}</td>' if measured else '<td class="num" data-v="">—</td>')
                + ('<td class="num" data-v="">not measured</td><td class="num" data-v="">not measured</td>'
                   if measured else '<td class="num" data-v="">—</td><td class="num" data-v="">—</td>')
                + '<td class="num cost" data-v="">—</td>'
                + f'<td class="txt">{esc(status)}: {esc(note)}</td></tr>')

    for t, note in ROWS:
        info = TOOLS[t]
        if t in settings:
            best = best_setting(summ, t)
            out.append(row(t, info["name"], info["kind"], best, SCORED_NOTES.get(t, SCORED_NOTE), "scored"))
            others = sorted(settings[t] - {best}, key=lambda s: -by[(t, s, "family")]["auc"])
            if others:
                out.append(f'<tr class="more-toggle" data-tool="{t}"><td class="txt" colspan="10">'
                           f'<button type="button" class="s40-more" data-tool="{t}" aria-expanded="false">'
                           f'+{len(others)} more {esc(info["name"])} settings</button></td></tr>')
                for s in others:
                    out.append(row(t, info["name"], info["kind"], s, SCORED_NOTES.get(t, SCORED_NOTE), "scored",
                                   cls=f"more more-{t}"))
        else:
            out.append(row(t, info["name"], info["kind"], "—", note,
                           "not applicable" if t == "hmmscan" else "not run"))
    for t, name, note in NOSEARCH:
        if t == "shuffled":
            for tool in sorted(x for x in settings if x.endswith("_shuffled")):
                base = tool[: -len("_shuffled")]
                best = best_setting(summ, base)
                for s in sorted(settings[tool], key=lambda s: (s != best, s)):
                    out.append(row(tool, f"{TOOLS[base]['name']}, queries shuffled", "nosearch", s,
                                   "each query shuffled with its dipeptide counts kept (k = 2, seed 1), same "
                                   "setting; other tools: not run yet", "chance level",
                                   cls="pinned" + ("" if s == best else f" more more-{tool}")))
                others = len(settings[tool]) - 1
                if others:
                    out.append(f'<tr class="more-toggle pinned" data-tool="{tool}"><td class="txt" colspan="10">'
                               f'<button type="button" class="s40-more" data-tool="{tool}" aria-expanded="false">'
                               f'+{others} more shuffled kmerseek settings</button></td></tr>')
            continue
        out.append(f'<tr data-tool="{t}" class="pinned"><td class="txt"><span class="kind k-nosearch" aria-hidden="true"></span> '
                   f'{esc(name)}</td><td class="txt">—</td>' + '<td class="num" data-v="">—</td>' * 4
                   + '<td class="num" data-v="">—</td><td class="num" data-v="">—</td><td class="num cost" data-v="">—</td>'
                   + f'<td class="txt">{esc(note)}</td></tr>')
    return "".join(out)


def page(data_dir: Path, r_cpu: float, r_gpu: float) -> tuple[str, dict]:
    dj = data_json(data_dir)
    summ = dj["summary"]
    tb = json.loads((data_dir / "scope40_testbed.json").read_text())
    scored_tools = sorted({r["tool"] for r in summ["overall"]})
    best = {t: best_setting(summ, t) for t in scored_tools}
    by = {(r["tool"], r["level"]): r for r in summ["overall"] if r["setting"] == best[r["tool"]]}
    fs, mm = {lvl: by[("foldseek", lvl)] for lvl in LEVELS}, {lvl: by[("mmseqs2_seqseq", lvl)] for lvl in LEVELS}
    km = {lvl: by[("kmerseek", lvl)] for lvl in LEVELS} if "kmerseek" in best else None
    ch = {lvl: by[("kmerseek_shuffled", lvl)] for lvl in LEVELS} if "kmerseek_shuffled" in best else None
    n_km = len({r["setting"] for r in summ["overall"] if r["tool"] == "kmerseek"})
    n_q, n_dom = summ["n_queries"], summ["n_domains"]
    miss = summ["missing"]["foldseek"]
    n_miss = len(miss)
    cls_n = {c: r["n_queries"] for r in summ["by_class"] if r["tool"] == "foldseek" and r["level"] == "family"
             for c in [r["scope_class"]]}
    absent = [c for c in CLASS_NAMES if c not in cls_n]
    rows = R.worked_example()
    fam, sf, fold, nest = rows
    prim = [r for r in summ["overall"] if r["setting"] == best[r["tool"]]]
    max_delta = max(abs(r["auc_hits_only"] - r["auc"]) for r in prim)
    min_cov = min(r["coverage"] for r in prim if not r["tool"].endswith("_shuffled"))
    kind_legend = "".join(
        f'<span class="item"><span class="kind k-{k}" aria-hidden="true"></span>{esc(KINDS[k]["label"])}</span>'
        for k in ("kmers", "letters", "profiles", "shapes", "models", "nosearch"))
    not_scored = [TOOLS[t]["name"] for t, _ in ROWS if t not in best and t != "hmmscan"]
    km_text = "" if km is None else (
        f" kmerseek is scored from its own hit lists: the one-row-per-pair tables that its SCOPe 2.08 runs "
        f"kept on Sherlock (PR 107), {n_km} kmerseek settings, ranked by region E-value and scored here "
        f"with Foldseek's rules. The chance level is the same kmerseek settings run on every query "
        f"shuffled with its dipeptide counts kept.")
    km_sentence = "" if km is None else (
        f" kmerseek's best setting, {best['kmerseek']}, reaches {km['family']['auc']:.3f} at the family "
        f"level and {km['superfamily']['auc']:.3f} at the superfamily level; the same setting on shuffled "
        f"queries reaches {ch['family']['auc']:.3f} and {ch['superfamily']['auc']:.3f}, the chance level."
        if ch else "")

    body = f"""<h1>SCOPe40: does a tool rank a domain's relatives first?</h1>
<p class="lede">Given one known domain, does the tool rank the domain's known relatives above the
first unrelated hit?</p>
<div class="prose">
<p>Table 1 places annotations on human proteins by searching one relative's proteome, and scores
where the ends of each call land. This page does something else: it ranks pairs of single
structural domains from SCOPe, and never scores boundaries. SCOPe sorts domains into
class.fold.superfamily.family; an id such as <span class="mono">b.1.1.1</span> is class b (all
beta), fold 1, superfamily 1, family 1. Domains in one superfamily are taken to share an ancestor;
domains in one fold share a shape, with or without one. The set searched is SCOPe 2.08 at 40%
identity, {n_dom:,} domains. Every domain is searched against all of them. The queries are the
{n_q:,} domains that have another member in the set at all three levels, as in Foldseek's
benchmark (<a href="methods.html#scope40">Methods</a>). Foldseek's own set was SCOPe 2.01, 11,211
domains, so these numbers are not the paper's.</p>
<p>For one query, a <b>true positive</b> (TP) is a hit to another domain at the level being scored.
At the family level that is the same family. At the superfamily level it is the same superfamily but
a different family, and at the fold level the same fold but a different superfamily. A <b>false positive</b> (FP) is a hit to a domain
of a different fold. Any other hit is <b>not counted</b>. Sort the hits by the tool's own score,
best first, after dropping the query's hit to itself, and stop at the first FP. The
<b>sensitivity up to the first false positive</b> is the TPs above that cut divided by the TPs the
query has in the set. Foldseek's scoring script, which this page follows, adds 1 to that count at
the superfamily and fold levels (the query itself stays in it); Figure {R.FIG} shows both.</p>
</div>

<figure class="ex-fig" id="fig-scope40-labels">
<figcaption><span class="ex-made-up">made-up example</span><b>Figure {R.FIG}.</b> One query,
<span class="mono">{R.EX_QUERY}</span>, and its eight hits, labelled at each level. The made-up set
holds {R.EX_N['family']} other domains in the query's family, {R.EX_N['superfamily']} in its
superfamily but outside the family, and {R.EX_N['fold']} in its fold but outside the superfamily.
Hits 1, 5 and 7 are three different domains of the query's family; the query's hit to itself was
dropped before hit 1. The small grey line under the superfamily and fold rows divides by one more,
as Foldseek's script does. The bottom row is the rule this site first planned, which also counts
the family at the superfamily level; it is not used.</figcaption>
{R.legend_html()}
<div class="ex-box">{R.figure_svg(rows)}</div>
<p class="ex-read">Hit {fam['cut'] + 1}, <span class="mono">{R.EX_HITS[fam['cut']]}</span>, is in
a different fold, so the walk stops there at every level. Family: hits 1 and 5,
{fam['above']} of {fam['n_possible']} = {fam['sens']:.2f}. Superfamily: hits 2 and 4,
{sf['above']} of {sf['n_possible']} = {sf['sens']:.2f} ({sf['above']} of {sf['n_script']} =
{sf['above'] / sf['n_script']:.2f} in Foldseek's script). Fold: hit 3, {fold['above']} of
{fold['n_possible']} = {fold['sens']:.2f} ({fold['above']} of {fold['n_script']} =
{fold['above'] / fold['n_script']:.2f}). Hit 7 is below the cut and adds nothing.</p>
</figure>

<div class="prose">
<p>Over many queries, sort the queries by sensitivity, highest first, and plot sensitivity against
the fraction of queries. The <b>area under the curve</b> (AUC) is the mean sensitivity over all
queries: the AUC and the mean are the same number. <b>Coverage</b> is the share of queries with at
least one hit at the settings the tool was run with.</p>
<p><b>Every query counts.</b> A query with no hits scores 0 and stays in every mean, every curve and
every area. An earlier kmerseek claim, that it "outperforms Foldseek on SCOPe", computed the area
over only the queries kmerseek returned hits for, 3.3% of them. That number is retired.</p>
<p><b>Inputs.</b> Foldseek runs on the SCOPe domain structures themselves: experimental
coordinates from the PDB, not predicted. ProstT5 + Foldseek and every other tool run on the domain
sequences. Foldseek and MMseqs2 are read from per-query files that other
authors published on SCOPe 2.08 (folder
<span class="mono">data/tea_scope40_rocx_files</span> of seanome/2024-kmerseek-analysis); those
files do not record which structures or settings the runs used. Their own scoring script dropped
each self-hit. {n_miss} of the {n_q:,} queries are missing from both files; they score 0 and are
listed under Table {TABLE}.{km_text} Not scored yet: {esc(", ".join(not_scored))}; their pipeline keeps
per-query hit lists and is written but not launched (<a href="{PIPELINE_PR}">pipeline PR</a>).</p>
</div>

<figure class="fgfig" id="fig-scope40-curves">
<figcaption><b>Figure {FIG_CURVES}.</b> Sensitivity up to the first false positive, one line per
tool, queries sorted from highest to lowest; a, family; b, superfamily; c, fold. The area under
each line is its AUC. n = {n_q:,} queries, SCOPe 2.08 at 40% identity.</figcaption>
<div class="legend" aria-label="Colour key">{kind_legend}</div>
<div class="s40-controls" id="s40-curve-controls">
<fieldset><legend>Tools</legend><span id="s40-tool-boxes"></span></fieldset>
<fieldset><legend>Queries with no hits</legend>
<label><input type="radio" name="s40-nohit" value="zero" checked> counted as 0 (default)</label>
<label><input type="radio" name="s40-nohit" value="drop"> dropped</label></fieldset>
</div>
<div id="s40-curve-legend"></div>
<div class="fg-grid-panels" id="s40-curves"></div>
<p class="fg-sentence">Foldseek, comparing shapes, finds family members for most queries
(AUC {fs['family']['auc']:.3f}) and keeps much of its lead at the superfamily level
({fs['superfamily']['auc']:.3f} against MMseqs2's {mm['superfamily']['auc']:.3f}). At the fold
level both are near the floor ({fs['fold']['auc']:.3f} and {mm['fold']['auc']:.3f}). A fold-level
TP is a domain of the same shape, which need not share an ancestor (Cheng et al. 2008).{km_sentence}
Dropping the no-hit queries moves a drawn AUC by at most {max_delta:.3f} here, because every
search tool drawn returns hits for at least {pct(min_cov)} of queries.</p>
</figure>

<figure class="fgfig" id="fig-scope40-classes">
<figcaption><b>Figure {FIG_CLASSES}.</b> Mean sensitivity up to the first false positive by SCOPe
class, with the 95% bootstrap interval as a bar (1,000 resamples of the queries, seed
{summ['seed']}); a, family; b, superfamily; c, fold. Classes:
{"; ".join(f"{c}, {CLASS_NAMES[c]} (n = {cls_n[c]:,})" for c in CLASS_NAMES if c in cls_n)}.
{"".join(f"Class {c} ({CLASS_NAMES[c]}) has no query with another member at all three levels, so it has no column. " for c in absent)}Click
a class to list its queries under the figure.</figcaption>
<div class="legend" aria-label="Marks" id="s40-class-legend"></div>
<div id="s40-classes"></div>
<p class="fg-sentence">Classes f (membrane and cell surface) and g (small proteins): composition
and cysteine spacing can match without homology there, so read those two columns with the
composition-scan row beside them. That row is not run yet.</p>
<div class="fg-list" id="s40-class-list"></div>
</figure>

<h2 id="s40-table">Table {TABLE}. Tool by level</h2>
<p class="prose">AUC with its 95% bootstrap interval under it, coverage, and cost. Click a heading
to sort; rows with no value sort last; the two no-search rows stay at the bottom. Core-hours and
GPU-hours come from a run's task records; the published files have none, so cost reads
<i>not measured</i>. Prices, as on Table 1: <label>r<sub>cpu</sub> $ per CPU-hour
<input type="number" id="s40-rcpu" step="0.001" min="0" value="{r_cpu:g}"></label>
<label>r<sub>gpu</sub> $ per GPU-hour <input type="number" id="s40-rgpu" step="0.01" min="0"
value="{r_gpu:g}"></label>.</p>
<div class="tablebox"><table class="lb" id="s40-tab"><thead><tr>
<th class="txt" scope="col"><button data-col="0">Tool<span class="arrow"></span></button></th>
<th class="txt" scope="col">Setting</th>
<th scope="col"><button data-col="2">Family AUC<span class="arrow"></span></button></th>
<th scope="col"><button data-col="3">Superfamily AUC<span class="arrow"></span></button></th>
<th scope="col"><button data-col="4">Fold AUC<span class="arrow"></span></button></th>
<th scope="col"><button data-col="5">Coverage<span class="arrow"></span></button></th>
<th scope="col"><button data-col="6">Core-hours<span class="arrow"></span></button></th>
<th scope="col"><button data-col="7">GPU-hours<span class="arrow"></span></button></th>
<th scope="col"><button data-col="8">Cost per 1,000 queries<span class="arrow"></span></button></th>
<th class="txt" scope="col">Status</th></tr></thead>
<tbody>{table_rows(summ)}</tbody></table></div>
<p class="prose">Queries missing from both published files, scored 0: {", ".join(f'<span class="mono">{esc(x)}</span>' for x in miss)}.
With the denominator over only the domains that can be found (no extra 1), the superfamily AUC
would be {fs['superfamily']['auc_possible']:.3f} for Foldseek and {mm['superfamily']['auc_possible']:.3f}
for MMseqs2; both are kept in <span class="mono">data/scope40_long.parquet</span>
(<span class="mono">sensitivity_possible</span>).</p>
<p class="refs">Cheng, H., B.-H. Kim, and N. V. Grishin. "MALISAM: A Database of Structurally
Analogous Motifs in Proteins." <i>Nucleic Acids Research</i> 36, Database issue (2008): D211–17.
<a href="https://doi.org/10.1093/nar/gkm698">https://doi.org/10.1093/nar/gkm698</a>. Other
references: <a href="methods.html#scope40">Methods</a>.</p>
<script defer src="assets/figures.js"></script>
<script defer src="assets/scope40.js"></script>
"""
    return body, dj
