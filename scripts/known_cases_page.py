"""The "Known cases" page: fourteen protein pairs whose shared function was proven by experiment.

Called from build_site.main with its layout(), so the page shares the header, fonts and footer.
Every number on the page is read from data/known_cases/cases.json and from the leaderboard's own
compute table; nothing here is typed by hand.
"""

from __future__ import annotations

import json
from pathlib import Path

# Figure numbers continue from the site's last figure (17, on the leaderboard and tool cards).
FIG_MATRIX = 18

# The four results a tool can have on a case, plus the mark for a search that has not been run.
# Each mark is a shape as well as a colour, and reuses a mark the site already has: a filled ink
# dot for found and a small dot for not found (as Figure 15 uses), the warning colour for a label
# that should not be transferred, and the grey × the whole site uses for "not run".
MARKS = {
    "found": {
        "glyph": '<circle cx="8" cy="8" r="5" fill="var(--ink)"/>',
        "words": "found: a call on the right target protein that overlaps the true region with "
                 "an IoU of 0.5 or more",
    },
    "wrong place": {
        "glyph": '<rect x="2.5" y="2.5" width="11" height="11" fill="none" stroke="var(--ink)" '
                 'stroke-width="1.8"/>',
        "words": "wrong place: a call on the right target protein, but its overlap with the true "
                 "region is an IoU below 0.5",
    },
    "wrong label": {
        "glyph": '<path d="M8,2 L14.5,13.5 L1.5,13.5 Z" fill="var(--warn)"/>',
        "words": "wrong label: the call carries the human protein's job across to a protein that "
                 "does not do that job",
    },
    "missed": {
        "glyph": '<circle cx="8" cy="8" r="1.6" fill="var(--muted)"/>',
        "words": "missed: no call on that protein above the tool's own cutoff",
    },
    "not run": {
        "glyph": '<path d="M3.5,3.5 L12.5,12.5 M12.5,3.5 L3.5,12.5" stroke="var(--muted)" '
                 'stroke-width="1.8"/>',
        "words": "not run: the search has not been run; each card says why",
    },
}


def esc(x) -> str:
    return (str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def mark_svg(result: str, size: int = 16) -> str:
    m = MARKS[result]
    return (f'<svg width="{size}" height="{size}" viewBox="0 0 16 16" aria-hidden="true">'
            f'{m["glyph"]}</svg>')


def legend() -> str:
    items = "".join(f'<span class="item">{mark_svg(k)}{esc(v["words"])}</span>'
                    for k, v in MARKS.items())
    return f'<div class="legend" aria-label="What each mark means">{items}</div>'


def case_title(c: dict) -> str:
    return f'Case {c["id"]}. {c["name"]}'


def refs_line(c: dict) -> str:
    return ", ".join(f'<a href="#ref-{n}">reference {n}</a>' for n in c["refs"])


def alignment_block(c: dict) -> str:
    a = c["alignment"]
    q, t, m = a["query"], a["target"], a["match"]
    qlab = c["queries"][0]["label"]
    tlab = c["targets"][0]["label"]
    width = max(len(qlab), len(tlab), 7)
    qs, ts = a["query_start"], a["target_start"]
    qe = qs + len(q.replace("-", "")) - 1
    te = ts + len(t.replace("-", "")) - 1
    rows = []
    step = 60
    cap = 180  # a long repeat block is cut here, and the note below says so
    cut = len(q) > cap
    if cut:
        q, t, m = q[:cap], t[:cap], m[:cap]
    for i in range(0, len(q), step):
        qi = qs + len(q[:i].replace("-", ""))
        ti = ts + len(t[:i].replace("-", ""))
        rows.append(f'{qlab:<{width}} {qi:>6} {q[i:i + step]}')
        rows.append(f'{"":<{width}} {"":>6} {m[i:i + step]}')
        rows.append(f'{tlab:<{width}} {ti:>6} {t[i:i + step]}')
        rows.append("")
    block = "\n".join(rows).rstrip()
    count = (f'{a["identical"]} of {a["columns"]} identical')
    if "identical_in_motif" in a:
        mc = a["motif_columns"]
        count = (f'{a["identical_in_motif"]} of {mc[1] - mc[0] + 1} identical in the motif, '
                 f'{a["identical"]} of {a["columns"]} identical over the whole block')
    cut_note = (f' The block is cut after the first {cap} columns; the counts are over the whole '
                f'alignment.' if cut else "")
    return (f'<pre class="aln"><code>{esc(block)}</code></pre>'
            f'<p class="alnnote">{esc(count)}.{cut_note} Residue numbers count from 1 at the '
            f'first residue of the UniProt entry. Alignment: {esc(a["source"])}. '
            f'Whole block: {esc(qlab)} {qs}–{qe}, {esc(tlab)} {ts}–{te}.</p>')


def tool_table(c: dict, tools: dict, kinds: dict, cost: dict) -> str:
    rows = []
    for tid, r in c["tools"].items():
        t = tools[tid]
        k = t["kind"]
        cc = cost.get(tid, {})
        score = "not run" if r["score"] is None else esc(r["score"])
        decoy = "not run" if r["decoy_score"] is None else esc(r["decoy_score"])
        call = "not run" if not r["call"] else esc(r["call"])
        iou = "not run" if r["iou"] is None else f'{r["iou"]:.2f}'
        inp = c["queries"][0]["input_per_tool"][tid]
        rows.append(
            f'<tr><td class="txt"><span class="kind k-{k}" aria-hidden="true"></span>'
            f'{esc(t["name"])} <span class="mono runs">{esc(r["setting"])}</span></td>'
            f'<td class="txt" title="{esc(r.get("reason", ""))}">{mark_svg(r["result"], 14)} '
            f'{esc(r["result"])}</td>'
            f'<td class="txt">{esc(inp["kind"])}</td>'
            f'<td class="num">{score}</td><td class="num">{decoy}</td>'
            f'<td class="num">{call}</td><td class="num">{iou}</td>'
            f'<td class="num"><span class="cost1k" data-cpu="{cc.get("cpu_h_per_1000", 0):.6f}" '
            f'data-gpu="{cc.get("gpu_h_per_1000", 0):.6f}">—</span></td></tr>')
    reasons = sorted({r["reason"] for r in c["tools"].values()
                      if r["result"] == "not run" and r.get("reason")})
    why = "".join(f'<p class="alnnote">Why a row reads not run: {esc(x)}.</p>' for x in reasons)
    return (
        '<div class="tablebox"><table class="lb"><thead><tr>'
        '<th class="txt">Tool and setting</th><th class="txt">Result</th>'
        '<th class="txt">Given</th><th>Its own score</th><th>Decoy score</th>'
        '<th>Call on the query</th><th>IoU</th>'
        '<th>Cost per 1,000 queries</th></tr></thead><tbody>'
        + "".join(rows) + "</tbody></table></div>" + why)


def matrix_svg(doc: dict, tools: dict) -> str:
    """Tools as rows, cases as columns grouped by box. Hover shows the score, the decoy's score
    and the coordinates; clicking scrolls to the card."""
    cases = doc["cases"]
    tool_ids = list(doc["tool_settings"])
    left, top = 210, 86
    colw, rowh = 34, 26
    gap = 16
    xs, x = {}, left
    last_box = None
    for c in cases:
        if last_box is not None and c["box"] != last_box:
            x += gap
        xs[c["id"]] = x
        x += colw
        last_box = c["box"]
    width = x + 16
    height = top + rowh * len(tool_ids) + 26
    out = [f'<svg viewBox="0 0 {width} {height}" role="img" '
           f'aria-label="Which tool finds which case, fourteen cases by nine tools">']
    # Box bands and case numbers along the top.
    for b in sorted({c["box"] for c in cases}):
        ids = [c["id"] for c in cases if c["box"] == b]
        x0 = xs[ids[0]] - 3
        x1 = xs[ids[-1]] + colw - 3
        out.append(f'<rect x="{x0}" y="{top - 20}" width="{x1 - x0}" '
                   f'height="{rowh * len(tool_ids) + 24}" fill="var(--chip)" rx="5"/>')
        out.append(f'<text x="{(x0 + x1) / 2:.0f}" y="{top - 30}" text-anchor="middle" '
                   f'font-size="12" fill="var(--muted)">Box {b}</text>')
    for c in cases:
        out.append(f'<text x="{xs[c["id"]] + colw / 2 - 3:.0f}" y="{top - 6}" '
                   f'text-anchor="middle" font-size="12">{c["id"]}</text>')
    for i, tid in enumerate(tool_ids):
        y = top + rowh * i + 17
        t = tools[tid]
        out.append(f'<text x="0" y="{y}" font-size="12">{esc(t["name"])}</text>')
        for c in cases:
            r = c["tools"][tid]
            cx = xs[c["id"]] + colw / 2 - 3
            cy = y - 4
            tip = (f'{t["name"]} on case {c["id"]}, {c["name"]}: {r["result"]}. '
                   f'Score {r["score"] or "not run"}, decoy {r["decoy_score"] or "not run"}, '
                   f'call {r["call"] or "not run"}.')
            glyph = MARKS[r["result"]]["glyph"]
            out.append(f'<g class="mk" transform="translate({cx - 8:.0f},{cy - 8:.0f})" '
                       f'data-case="{c["id"]}" tabindex="0" role="button" '
                       f'aria-label="{esc(tip)}"><title>{esc(tip)}</title>'
                       f'<rect width="16" height="16" fill="transparent"/>{glyph}</g>')
    out.append("</svg>")
    return '<div class="plotbox kcmatrix">' + "".join(out) + "</div>"


def card(c: dict, doc: dict, tools: dict, kinds: dict, cost: dict) -> str:
    box = doc["boxes"][str(c["box"])]
    eco = doc["eco"][c["eco"]]
    q = c["queries"][0]
    t = c["targets"][0]
    if c.get("elm"):
        names = f'{len(c["queries"])} pathogen proteins, listed below'
        tnames = f'{len({x["acc"] for x in c["targets"]})} human proteins, listed below'
    else:
        names = ", ".join(f'{esc(x["label"])} (<span class="mono">{esc(x["acc"])}</span>, '
                          f'{esc(x["species"])})' for x in c["queries"])
        tnames = ", ".join(f'{esc(x["label"])} (<span class="mono">{esc(x["acc"])}</span>, '
                           f'{esc(x["species"])})' for x in c["targets"])
    st = c.get("structure")
    given = (f'Foldseek is given PDB {esc(st["pdb"])} chain {esc(st["chain"])}: {esc(st["what"])}'
             if st else
             f'Foldseek is given the AlphaFold Database model of '
             f'<span class="mono">{esc(q["acc"])}</span>, since the papers cited here report no '
             f'experimental structure of this protein')
    notes = []
    if c.get("preprint"):
        notes.append("The evidence for this pair is a preprint, not a peer-reviewed paper.")
    for k in ("structure_note", "motif_note", "target_rule", "pick_rule"):
        if c.get(k):
            notes.append(c[k])
    if c.get("wrong_label"):
        notes.append('A tool that transfers "' + c["wrong_label"] + '" to this protein has made '
                     "a wrong-label call, however well it places the region.")
    if c.get("earlier_runs"):
        notes.append(c["earlier_runs"])
    note_html = "".join(f"<p class=\"muted\">{esc(n)}</p>" for n in notes)
    elm_list = ""
    if c.get("elm"):
        items = "".join(
            f'<li><span class="mono">{esc(x["elm_instance"])}</span> '
            f'{esc(x["elm_class"])} in {esc(x["label"])} '
            f'(<span class="mono">{esc(x["acc"])}</span>, {esc(x["species"])}), residues '
            f'{x["region"][0]}–{x["region"][1]}, partner '
            f'<span class="mono">{esc(x["partner"]["acc"])}</span></li>'
            for x in c["queries"])
        elm_list = f'<ol class="elm">{items}</ol>'
    return f"""
<section class="card kc" id="case-{c["id"]}">
<h2>{esc(case_title(c))}</h2>
<dl class="facts">
<dt>Searched</dt><dd>{names}</dd>
<dt>Label at stake</dt><dd>{tnames}</dd>
<dt>Box {c["box"]}</dt><dd>{esc(box["title"])}. A right call names {esc(box["right"])}.</dd>
<dt>The experiment</dt><dd>{esc(c["experiment"])} ({refs_line(c)}, evidence code
<span class="mono">{esc(c["eco"])}</span>, {esc(eco)}).</dd>
<dt>Given to Foldseek</dt><dd>{given}.</dd>
</dl>
{elm_list}
<h3>The residues</h3>
{alignment_block(c)}
<h3>What each tool does with it</h3>
{tool_table(c, tools, kinds, cost)}
{note_html}
<p class="takeaway"><b>For an annotator:</b> {esc(c["takeaway"])}</p>
</section>"""


def page(layout, src: dict, data_dir: Path, tools: dict, kinds: dict, cost: dict,
         out_dir: Path) -> None:
    doc = json.loads((data_dir / "known_cases" / "cases.json").read_text())
    n_cases = len(doc["cases"])
    n_queries = sum(len(c["queries"]) for c in doc["cases"])
    n_run = sum(1 for c in doc["cases"] for r in c["tools"].values() if r["result"] != "not run")
    boxes = "".join(
        f'<li><b>Box {b}. {esc(v["title"])}.</b> {esc(v["words"])}</li>'
        for b, v in sorted(doc["boxes"].items()))
    refs = "".join(
        f'<li id="ref-{r["n"]}">{r["text"]} '
        f'<a href="https://doi.org/{r["doi"]}">https://doi.org/{esc(r["doi"])}</a>.'
        + (f' <b>{esc(r["note"])}</b>' if r.get("note") else "") + "</li>"
        for r in doc["references"])
    cards = "".join(card(c, doc, tools, kinds, cost) for c in doc["cases"])
    n_pairs = sum(len(c["queries"]) * len(doc["tool_settings"]) for c in doc["cases"])
    if n_run:
        run_line = f"Results are in for {n_run} of the {n_pairs} searches."
    else:
        run_line = ("<b>None of these searches has been run yet.</b> Each result lands in its "
                    "own place on the cards below as the runs finish; until then every mark "
                    "reads not run and says why.")
    body = f"""<h1>Known cases</h1>
<p class="lede">Fourteen pairs of proteins whose shared job was proven in a laboratory, and what
each of the nine tools on this board does with them. This page shows illustrations, not rates.
The rates stay in Table 1.</p>
<p class="prose">Each case is one query protein searched against {esc(doc["target_db"])}. That is
the opposite direction from Table 1, which searches human proteins against one relative's
proteome.</p>
<p class="prose">Three words this page uses throughout. <b>IoU</b> (intersection over union) is
the number of residues a call and the true region share, divided by the number either one
covers. The site counts a call as correct at an IoU of 0.5 or more. An <b>E-value</b> is how
many hits as good as this one a search of a database this size would be expected to turn up by
chance, so smaller is better. A <b>decoy</b> is the same query with its residues reshuffled so
that every pair of neighbouring residues still occurs as often as before. A hit counts only
when the real query beats its decoy. Every query here has one decoy, shuffled with seed
{doc["shuffle_seed"]}.</p>
<p class="prose">{n_queries} query proteins over {n_cases} cases, each searched by
{len(doc["tool_settings"])} tools at the setting Table 1 calls that tool's best.
{run_line}</p>

<h2 id="boxes">The four boxes</h2>
<p class="prose">Each case sits in one of four boxes, and the box decides what counts as a right
call.</p>
<ol class="boxes">{boxes}</ol>

<h2 id="matrix">Which tool finds which case</h2>
{legend()}
<figure>
<figcaption><b>Figure {FIG_MATRIX}.</b> One mark per tool and case, grouped by box. Hovering a
mark shows the score, the decoy's score and the coordinates of the call; clicking it opens that
case below. Numbers along the top are the case numbers of the cards below.</figcaption>
{matrix_svg(doc, tools)}
</figure>

<h2 id="cost">What a search costs</h2>
<div class="prices"><label>Price of a core-hour, US dollars
<input id="kc-rcpu" type="number" min="0" step="0.001" value="{cost["r_cpu"]:.6f}"></label>
<label>Price of a GPU-hour, US dollars
<input id="kc-rgpu" type="number" min="0" step="0.01" value="{cost["r_gpu"]:.3f}"></label></div>
<p class="prose">The cost on every card is core-hours × the price of a core-hour, plus GPU-hours ×
the price of a GPU-hour, the same arithmetic the cost page uses. Until a case has been run, the
core-hours are the ones the Table 1 run measured for that tool against one relative's proteome.
Swiss-Prot is larger than any one proteome, so those numbers are a floor.</p>

{cards}

<h2 id="references">References</h2>
<ol class="refs">{refs}</ol>
"""
    html = layout("Known cases", "known-cases", body, src,
                  head='<script defer src="assets/known-cases.js"></script>',
                  description="Fourteen protein pairs whose shared function was proven by "
                              "experiment, and what each tool on this board does with them.")
    (out_dir / "known-cases.html").write_text(html)
    (out_dir / "data").mkdir(exist_ok=True)
    (out_dir / "data" / "known_cases.json").write_text(json.dumps(doc, indent=1,
                                                                  ensure_ascii=False) + "\n")
