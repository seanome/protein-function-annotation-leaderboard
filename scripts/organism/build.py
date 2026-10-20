"""Write the organism pages into docs/: the home page ("I work on"), one page per pre-built
species, the page that builds any other species in the browser, and the tables they read.

Called from build_site.py when data/organism/lookups.jsonl exists. Inputs:
  data/results.parquet, data/compute.parquet, data/pileup/   (the leaderboard's tables)
  data/organism/                                             (lookups, from fetch.py and lookups.py)
"""

import json
import shutil
from html import escape as esc
import subprocess
import tempfile
from pathlib import Path

import strata as strata_mod

ROOT = Path(__file__).resolve().parents[2]
PRICES_NOTE = ("c7i.large $0.08925 per hour for 2 vCPU and g5.xlarge $1.006 per GPU-hour, "
               "US East, read 2026-10-15")


def head(depth):
    up = "../" * depth
    return (f'<script defer src="{up}assets/organism.js"></script>\n'
            f'<script defer src="{up}assets/organism-page.js"></script>')


def build(data: Path, out: Path, bs) -> int:
    od = data / "organism"
    src = json.loads((data / "SOURCE.json").read_text())
    sources = json.loads((od / "sources.json").read_text())
    S = strata_mod.build(data)
    S["uniprot_release"] = sources["swissprot_per_taxon"]["release"]
    S["timetree_release"] = sources["timetree_mrca"]["release"].replace("TimeTree 5, ", "")
    S["goat_release"] = sources["goat_ebp_list"]["release"]
    S["prices_note"] = PRICES_NOTE
    dl = json.loads((od / "web/dark_lineages.json").read_text())
    for d in S["dark"]:
        d["lineage"] = dl[d["organism"]]
    S["lookup_sources"] = sources

    web = out / "data/organism"
    if web.exists():
        shutil.rmtree(web)
    (web / "names").mkdir(parents=True)
    (web / "strata.json").write_text(json.dumps(S, separators=(",", ":"), allow_nan=False))
    for f in ("reviewed.json", "mrca.json"):
        shutil.copyfile(od / "web" / f, web / f)
    for f in (od / "web/names").glob("*.json"):
        shutil.copyfile(f, web / "names" / f.name)

    recs = [json.loads(x) for x in (od / "lookups.jsonl").read_text().splitlines() if x.strip()]
    seen = {}
    for r in recs:
        if r["slug"] in seen:
            r["slug"] = f'{r["slug"]}_{r["taxid"]}'
        seen[r["slug"]] = r["taxid"]
    (web / "prebuilt.json").write_text(json.dumps({str(r["taxid"]): r["slug"] for r in recs}, separators=(",", ":")))

    with tempfile.TemporaryDirectory() as td:
        sp, lp, ap = Path(td) / "s.json", Path(td) / "l.jsonl", Path(td) / "a.jsonl"
        sp.write_text(json.dumps(S, allow_nan=False))
        lp.write_text("\n".join(json.dumps(r) for r in recs) + "\n")
        subprocess.run(["node", str(ROOT / "scripts/organism/render.cjs"), str(sp), str(lp), str(ap)], check=True)
        rendered = [json.loads(x) for x in ap.read_text().splitlines() if x.strip()]

    pages = out / "organism"
    if pages.exists():
        shutil.rmtree(pages)
    pages.mkdir()
    index_rows = []
    for rec, r in zip(recs, rendered):
        A, html_body = r["answer"], r["html"]
        name = f'{rec["common_name"]} ({rec["name"]})' if rec["common_name"] else rec["name"]
        body = (f'<nav class="crumb" aria-label="Breadcrumb"><a href="../index.html">Your organism</a> › {bs.esc(name)}</nav>\n'
                + html_body + sources_html(S, src))
        rel = f'organism/{rec["slug"]}.html'
        (out / rel).write_text(bs.layout(f"{name}: which tools to run", "organism", body, src, depth=1,
                                         head=head(1), description=f"Which annotation tools to run on {name}, "
                                         "read from measured results, with a command to copy."))
        bs.page_json(out, rel, src, {"lookups": rec, "answer": A})
        index_rows.append((name, rec["slug"], A))

    (out / "index.html").write_text(page_search(bs, src))
    bs.page_json(out, "index.html", src, {"prebuilt_species": len(recs),
                                          "pages": "organism/index.html.json lists them"})
    (pages / "live.html").write_text(page_live(bs, src))
    (pages / "index.html").write_text(page_list(bs, src, index_rows))
    bs.page_json(out, "organism/index.html", src, {"pages": [
        {"name": n, "page": f"{s}.html", "json": f"{s}.html.json"} for n, s, _ in index_rows]})
    return len(recs)


def sources_html(S, src):
    rows = []
    for k, v in sorted(S["lookup_sources"].items()):
        rows.append(f'<tr><td>{k}</td><td>{esc(v.get("release") or "")}</td><td>{v["fetched_utc"][:10]}</td>'
                    f'<td class="mono small">{esc(v["url"])}</td></tr>')
    return (f'<h2 id="sources">Where every number comes from</h2>'
            f'<p class="small">Tool scores and costs: results.parquet and compute.parquet, {src["repo"]} commit '
            f'<span class="mono">{src["commit"][:8]}</span>, experiment {strata_ex(S)}. cost_usd_est is an estimate: '
            f'CPU-hours times the AWS on-demand price ({PRICES_NOTE}), plus GPU-hours times the GPU price. '
            f'Lookups for the four facts:</p>'
            f'<div class="tablewrap"><table class="mini"><thead><tr><th>lookup</th><th>release</th><th>fetched</th><th>address</th></tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>')


def strata_ex(S):
    return S["truth"].split("experiment ")[-1]


def page_search(bs, src):
    body = """
<section class="hero" aria-labelledby="iwork">
<h1 id="iwork"><label for="species">I work on</label></h1>
<div class="combo">
<input id="species" type="text" autocomplete="off" spellcheck="false" placeholder="a species or a common name, for example Pacific oyster"
  role="combobox" aria-autocomplete="list" aria-expanded="false" aria-controls="suggest">
<ul id="suggest" role="listbox" hidden></ul>
<p id="suggest-status" class="small muted" aria-live="polite"></p>
</div>
</section>
<section class="below">
<p class="lede">Type the species you study. Its page gives three sentences and a command: which
annotation tools to run on its proteins, in what order, and how much of its proteome no tool here
will label.</p>
<p>Reference cases: <a href="organism/homo_sapiens.html">human</a>, <a href="organism/mus_musculus.html">mouse</a>,
<a href="organism/saccharomyces_cerevisiae.html">baker's yeast</a>, <a href="organism/escherichia_coli.html">E. coli</a>,
<a href="organism/arabidopsis_thaliana.html">thale cress (Arabidopsis)</a>. <a href="organism/index.html">All pre-built pages</a>.
The <a href="leaderboard.html">leaderboard</a> has every tool and setting.</p>
</section>
"""
    return bs.layout("Which annotation tools to run on your organism", "organism", body, src, depth=0,
                     head=head(0),
                     description="Type a species and get the annotation tools to run on its proteins, "
                                 "read from measured results.").replace(
        '<body data-page="organism">', '<body data-page="organism" data-organism="search" data-base="">')


def page_live(bs, src):
    body = """<nav class="crumb" aria-label="Breadcrumb"><a href="../index.html">Your organism</a> › built in this browser</nav>
<p id="progress" class="small muted" aria-live="polite">Looking up…</p>
<div id="answer"></div>
<p class="small"><a id="json-link" hidden>This answer as JSON</a></p>
<noscript><p>This page builds the answer in the browser and needs JavaScript. Pre-built pages: <a href="index.html">all species</a>.</p></noscript>
"""
    return bs.layout("Built in this browser: which tools to run", "organism", body, src, depth=1, head=head(1),
                     description="The organism page for a species with no pre-built page, built from live lookups.") \
        .replace('<body data-page="organism">', '<body data-page="organism" data-organism="live" data-base="../">')


def page_list(bs, src, rows):
    items = "".join(f'<li><a href="{s}.html">{bs.esc(n)}</a></li>' for n, s, _ in sorted(rows, key=lambda x: x[0].lower()))
    body = f"""<nav class="crumb" aria-label="Breadcrumb"><a href="../index.html">Your organism</a> › all pre-built pages</nav>
<h1>All pre-built organism pages</h1>
<p class="lede">{len(rows):,} species: every species on the Earth BioGenome Project target list that has a UniProt
reference proteome, and the five reference cases. Any other species is built in the browser from the box on
<a href="../index.html">the first page</a>.</p>
<ul class="splist">{items}</ul>"""
    return bs.layout("All pre-built organism pages", "organism", body, src, depth=1)
