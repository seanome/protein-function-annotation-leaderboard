#!/usr/bin/env python3
"""Headless check of the built site against data/results.parquet and data/compute.parquet.

1. Leaderboard: render every chip combination (truth set x dedup x relative x stratum).
   For each tool and setting row (all are in the page, the other settings hidden), the page
   must show the parquet's numbers when the parquet has a row, and "× not run" in every cell
   when it has none. The visible row of each tool must be its best setting by Fmax, and the
   sentence above the table must name the best tool by Fmax, the cheapest within 0.01 of
   it, and the best by AUPRC, all recomputed here. The every-setting view is checked too.
2. Compare two columns (Figure 1 on cost.html): render every combination for each of the
   three scores against core-hours, and every combination with no stratum chosen for five
   more axis pairs (cost against AUPRC, peak RAM against median IoU, Fmax against AUPRC,
   wall-clock time against Fmax, database build against Fmax). The plotted markers must be
   exactly the rows with a value on both axes, each marker's x and y must equal the value
   computed here from the parquet files, every other row must be listed as not plotted, and
   the line through the best trade-offs must be drawn only for a cost against a score.
3. Every page: no console errors, no broken internal links, and no horizontal page scroll
   at phone width (375 px). Every number written before "Mya" or "million years" must be a
   relative's divergence time from human in the parquet's target_mya_to_human, rounded.
4. Divergence times: each relative's chip label and the page's species table show the
   parquet's target_mya_to_human rounded to whole million years, and no number where it is
   null (E. coli).

Uses the installed Google Chrome through Playwright (channel="chrome").
"""

import argparse
import functools
import http.server
import itertools
import json
import re
import socketserver
import sys
import threading
from pathlib import Path
from urllib.parse import urljoin, urlparse

import polars as pl
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]


def serve(site: Path) -> tuple[socketserver.TCPServer, str]:
    class Quiet(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

    handler = functools.partial(Quiet, directory=str(site))
    httpd = socketserver.TCPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd, f"http://127.0.0.1:{httpd.server_address[1]}/"


def combos(meta: dict) -> list[dict]:
    rows = {r["id"]: r for r in meta["chiprows"]}
    vals = lambda rid: [c["value"] for g in rows[rid]["groups"] for c in g["chips"]]  # noqa: E731
    strata = [("all", "all")]
    for rid in ("identity", "length", "plddt", "ftype"):
        strata += [(rows[rid]["axis"], v) for v in vals(rid) if v != "all"]
    return [
        {"truth": t, "calls": c, "species": s, "axis": a, "stratum": st}
        for t, c, s, (a, st) in itertools.product(vals("truth"), vals("calls"), vals("species"), strata)
    ]


# Runs inside the page: render a batch of filter combinations and return what each shows.
SNAP_TABLE = """(batch) => batch.map((s) => {
  window.__lb.set(s);
  const rows = [...document.querySelectorAll('#lb tbody tr')].map((tr) => {
    const cells = {};
    for (const td of tr.querySelectorAll('td[data-col]')) {
      cells[td.dataset.col] = {text: td.textContent.trim(), missing: td.dataset.missing || null, none: td.dataset.none || null};
    }
    return {tool: tr.dataset.tool, setting: tr.dataset.setting, hidden: tr.hidden, cells};
  });
  const v = document.getElementById('verdict').dataset;
  return {s, rows, verdict: {best: v.best, cheap: v.cheap, auprc: v.auprc}};
})"""

SNAP_PLOT = """(batch) => batch.map((s) => {
  window.__lb.set(s);
  const marks = [...document.querySelectorAll('#plot [data-tool]')].map((m) => m.dataset.tool + '|' + m.dataset.setting);
  const notplotted = [...document.querySelectorAll('#notplotted .np')].map((m) => m.dataset.tool + '|' + m.dataset.setting);
  const st = window.__lb.state();
  return {s, marks, notplotted, plotted: window.__lb.plotted(), state: {x: st.x, y: st.y},
          frontier: !!document.querySelector('#plot polyline.frontier')};
})"""


def fmt3(x):
    return f"{x:.3f}"


# Runs inside alphabets.html: everything the page shows, with the data attributes behind it.
SNAP_ALPHABETS = """() => {
  const D = (el) => Object.assign({}, el.dataset);
  const rows = [...document.querySelectorAll('table.al tbody tr[data-alphabet]')].map((tr) => ({
    alphabet: tr.dataset.alphabet,
    classes: [...tr.querySelectorAll('td.classes .cls')].map((c) => ({
      letters: [...c.querySelectorAll('.sq')].map((s) => s.textContent).join(''),
      hp: ['h', 'p', 'c'].find((k) => c.classList.contains('hp-' + k)) || null})),
    nClasses: tr.querySelector('th .key').textContent,
    ev: [...tr.querySelectorAll('td.ev')].map((td) => {
      const L = [...td.querySelectorAll('line')], c = td.querySelector('circle');
      return {...D(td), text: td.querySelector('.v').textContent, grey: td.classList.contains('grey'), title: td.title,
              one: +L[0].getAttribute('x1'), ciLo: +L[1].getAttribute('x1'), ciHi: +L[1].getAttribute('x2'), dot: +c.getAttribute('cx')};
    }),
    kp: [...tr.querySelectorAll('td.kp')].map((td) => ({...D(td), text: td.textContent.trim(),
                                                        grey: td.classList.contains('grey'), title: td.title})),
    uk: [...tr.querySelectorAll('td.uk .kb')].map((s) => {
      const L = [...s.querySelectorAll('line')], R = [...s.querySelectorAll('rect')];
      return {...D(s), text: s.querySelector('.v').textContent, hidden: s.hidden, title: s.title,
              rangeLo: +L[0].getAttribute('x1'), rangeHi: +L[0].getAttribute('x2'),
              bars: R.map((r) => [+r.getAttribute('x'), +r.getAttribute('x') + +r.getAttribute('width')])};
    }),
  }));
  const fullHeads = [...document.querySelectorAll('table.full tbody tr')].map((tr) => tr.querySelector('th').textContent);
  const hpText = document.querySelector('#hp + p').textContent.replace(/\\s+/g, ' ');
  const evHead = [...document.querySelectorAll('th.evh')].map((th) => ({...D(th), text: th.querySelector('.key').textContent}));
  const kpHead = [...document.querySelectorAll('th.kph')].map((th) => ({...D(th), text: th.querySelector('.key').textContent}));
  const hp = [...document.querySelectorAll('table.hp116 td.sq116')].map((td) => ({...D(td), text: td.textContent}));
  const hpN = [...document.querySelectorAll('table.hp116 tbody tr')].map((tr) => ({
    alphabet: tr.querySelector('th').textContent, n: tr.querySelector('td.num').textContent}));
  const full = [...document.querySelectorAll('table.full td.full')].map((td) => ({...D(td), text: td.textContent}));
  return {rows, evHead, kpHead, hp, hpN, full, fullHeads, hpText, lead: document.getElementById('al-lead').textContent.replace(/\\s+/g, ' '),
          kbinN: document.getElementById('kbin-n').textContent};
}"""


def check_alphabets(page, base: str, data: Path, site: Path, failures: list, stats: dict) -> None:
    """alphabets.html against data/alphabets/*.parquet and PR 116's CSV, read here directly."""
    import csv
    d = data / "alphabets"
    ev = pl.read_parquet(d / "evidence.parquet")
    kp = pl.read_parquet(d / "kappa.parquet")
    uk = pl.read_parquet(d / "usable_k.parquet")
    reach = pl.read_parquet(d / "k_reach.parquet")
    cl = pl.read_parquet(d / "classes.parquet")
    with open(d / "hp_borderline_residues.csv") as f:
        hp116 = list(csv.DictReader(f))
    page.goto(base + "alphabets.html")
    snap = page.evaluate(SNAP_ALPHABETS)
    n = {"numbers": 0}
    fair = ev.filter(pl.col("mode") == "low_half")
    groups = ["Pfam seed", "SCOPe: same family", "SCOPe: same superfamily, different family"]
    alphabets = sorted(cl["alphabet"].unique())
    if sorted(r["alphabet"] for r in snap["rows"]) != alphabets:
        failures.append(f"alphabets.html: rows {[r['alphabet'] for r in snap['rows']]} are not the 19 alphabets")
    for row in snap["rows"]:
        a = row["alphabet"]
        # classes: the same residues grouped together as classes.parquet
        want = {frozenset(m) for m in cl.filter(pl.col("alphabet") == a)["class_members"].unique()}
        got = {frozenset(c["letters"]) for c in row["classes"]}
        if got != want or sum(len(c["letters"]) for c in row["classes"]) != 20:
            failures.append(f"alphabets.html {a}: classes {sorted(map(sorted, got))} != {sorted(map(sorted, want))}")
        for c in row["classes"]:
            hp = cl.filter((pl.col("alphabet") == a) & (pl.col("residue") == c["letters"][0]))["hp_class"][0]
            if c["hp"] != hp:
                failures.append(f"alphabets.html {a}: class {c['letters']} drawn as {c['hp']}, table says {hp}")
        # evidence, fair test
        if [e["group"] for e in row["ev"]] != groups:
            failures.append(f"alphabets.html {a}: evidence columns {[e['group'] for e in row['ev']]}")
        for e in row["ev"]:
            r = fair.filter((pl.col("group") == e["group"]) & (pl.col("alphabet") == a)).row(0, named=True)
            text = f"{r['ratio_to_20aa']:.3f} ({r['lo']:.3f}–{r['hi']:.3f})"
            if (e["text"] != text or float(e["ratio"]) != r["ratio_to_20aa"] or float(e["lo"]) != r["lo"]
                    or float(e["hi"]) != r["hi"] or int(e["n"]) != r["n_items"] or int(e["clusters"]) != r["n_clusters"]
                    or e["grey"] != (r["n_items"] < 30)):
                failures.append(f"alphabets.html {a} {e['group']}: page {e} != parquet {text}, n {r['n_items']}")
            n["numbers"] += 3
        # kappa
        bins = ["<20%", "20-30%", "30-40%", "40-60%", ">=60%"]
        if [k["bin"] for k in row["kp"]] != bins:
            failures.append(f"alphabets.html {a}: kappa bins {[k['bin'] for k in row['kp']]}")
        for k in row["kp"]:
            r = kp.filter((pl.col("alphabet") == a) & (pl.col("identity_bin") == k["bin"])).row(0, named=True)
            if (k["text"] != f"{r['kappa']:.2f}" or float(k["kappa"]) != r["kappa"] or int(k["n"]) != r["n_pairs"]
                    or float(k["agreement"]) != r["agreement"] or float(k["chance"]) != r["chance"]
                    or k["grey"] != (r["n_pairs"] < 30)):
                failures.append(f"alphabets.html {a} kappa {k['bin']}: page {k} != parquet {r}")
            n["numbers"] += 1
        # usable k, all three bins (two hidden)
        for s in row["uk"]:
            r = uk.filter((pl.col("alphabet") == a) & (pl.col("identity_bin") == s["bin"])).row(0, named=True)
            # recomputed from the share of pairs reached at each k: inside notebook 274's range and >= 5%
            rr = reach.filter((pl.col("alphabet") == a) & (pl.col("identity_bin") == s["bin"]))
            ok = rr.filter(pl.col("k").is_between(r["k_min"], r["k_max"]) & (pl.col("share_reached") >= 0.05))["k"]
            most = rr.filter(pl.col("k").is_between(r["k_min"], r["k_max"]) & (pl.col("share_reached") >= 0.20))["k"]
            lo, hi = ok.min(), ok.max()
            text = str(lo) if lo == hi else f"{lo}–{hi}"
            if (s["text"] != text or int(s["kLo"]) != lo or int(s["kHi"]) != hi or (lo, hi) != (r["k_use_lo"], r["k_use_hi"])
                    or (s["kMost"] or None) != (str(most.max()) if most.len() else None)
                    or int(s["kMin"]) != r["k_min"] or int(s["kMax"]) != r["k_max"] or int(s["n"]) != r["n_pairs"]
                    or s["hidden"] != (s["bin"] != "20-30%")):
                failures.append(f"alphabets.html {a} usable k {s['bin']}: page {s} != {text}, {r}")
            n["numbers"] += 2
        if [s["bin"] for s in row["uk"]] != ["20-30%", "30-40%", "40-60%"]:
            failures.append(f"alphabets.html {a}: usable k bins {[s['bin'] for s in row['uk']]}")
    for h in snap["evHead"]:
        e = fair.filter(pl.col("group") == h["group"])
        if not h["text"].startswith(f"{e['n_items'][0]:,} pair-halves, {e['n_clusters'][0]:,} "):
            failures.append(f"alphabets.html header {h}: n is not {e['n_items'][0]:,}")
    for h in snap["kpHead"]:
        if h["text"] != f"{kp.filter((pl.col('identity_bin') == h['bin']) & (pl.col('alphabet') == 'protein20'))['n_pairs'][0]:,}":
            failures.append(f"alphabets.html kappa header {h}")
    if snap["kbinN"] != f"n = {uk.filter(pl.col('identity_bin') == '20-30%')['n_pairs'][0]:,} pairs":
        failures.append(f"alphabets.html usable k header n {snap['kbinN']}")
    # PR 99's comment of 2026-10-14, Pfam seed fair test: sdm12 1.032 (1.027-1.038), hp_thomas_dill2 0.79
    by = {(r["alphabet"], e["group"]): e["text"] for r in snap["rows"] for e in r["ev"]}
    pr99 = {"sdm12": by[("sdm12", "Pfam seed")], "hp_thomas_dill2": by[("hp_thomas_dill2", "Pfam seed")]}
    if pr99["sdm12"] != "1.032 (1.027–1.038)" or not pr99["hp_thomas_dill2"].startswith("0.790 "):
        failures.append(f"alphabets.html: Pfam fair-test rows {pr99} differ from PR 99's comment")
    if "1.032" not in snap["lead"] or "1.027–1.038" not in snap["lead"] or "0.79 " not in snap["lead"]:
        failures.append(f"alphabets.html: lead sentence {snap['lead']!r}")
    above = fair.filter((pl.col("group") == "Pfam seed") & (pl.col("lo") > 1) & (pl.col("alphabet") != "protein20")).sort("rank")
    if f"({', '.join(above['alphabet'])})" not in snap["lead"]:
        failures.append(f"alphabets.html: lead sentence does not list {list(above['alphabet'])}")
    stats["alphabets_pr99_rows"] = pr99
    # PR 116's CSV, all 140 squares, and the residues-in-h column
    want = {(r["alphabet"], r["residue"]): r["class"] for r in hp116}
    got = {(s["alphabet"], s["residue"]): (s["class"], s["text"]) for s in snap["hp"]}
    if len(snap["hp"]) != 140 or len(want) != 140 or {k: v[0] for k, v in got.items()} != want \
            or any(c != t for c, t in got.values()):
        failures.append(f"alphabets.html: H/P squares differ from PR 116's CSV ({len(snap['hp'])} squares)")
    n_h = {r["alphabet"]: r["n_residues_in_h"] for r in hp116}
    if {r["alphabet"]: r["n"] for r in snap["hpN"]} != n_h:
        failures.append(f"alphabets.html: residues in h {snap['hpN']} != {n_h}")
    stats["alphabets_hp116_squares_matched"] = sum(1 for k, v in got.items() if want.get(k) == v[0])
    # every group and pair selection
    for c in snap["full"]:
        r = ev.filter((pl.col("group") == c["group"]) & (pl.col("mode") == c["mode"]) & (pl.col("alphabet") == c["alphabet"])).row(0, named=True)
        if c["text"] != f"{r['ratio_to_20aa']:.2f}" or int(c["n"]) != r["n_items"]:
            failures.append(f"alphabets.html full table {c} != {r['ratio_to_20aa']}")
        n["numbers"] += 1
    if len(snap["full"]) != ev.height:
        failures.append(f"alphabets.html: full table has {len(snap['full'])} cells, parquet {ev.height} rows")
    check_alphabets_words(snap, ev, kp, uk, reach, cl, hp116, failures, n)
    # the chips switch the usable-k bin, and the n in the column header follows
    page.click('.kbins .chip[data-kbin="40-60%"]')
    vis = page.evaluate("[...document.querySelectorAll('td.uk .kb')].filter(s => !s.hidden).map(s => s.dataset.bin)")
    if set(vis) != {"40-60%"} or len(vis) != 19:
        failures.append(f"alphabets.html: after the 40-60% chip, visible bins {set(vis)} ({len(vis)})")
    want_n = f"n = {uk.filter(pl.col('identity_bin') == '40-60%')['n_pairs'][0]:,} pairs"
    if page.evaluate("document.getElementById('kbin-n').textContent") != want_n:
        failures.append(f"alphabets.html: usable k header after the 40-60% chip is not {want_n!r}")
    # the JSON twin's data equal the parquet files
    doc = json.loads((site / "alphabets.html.json").read_text())["data"]
    for name, df in (("evidence", ev), ("kappa", kp), ("usable_k", uk)):
        cols = list(doc[name][0])
        key = lambda r: json.dumps(r, sort_keys=True)  # noqa: E731
        if len(doc[name]) != df.height or not set(cols) <= set(df.columns) \
                or sorted(map(key, doc[name])) != sorted(map(key, df.select(cols).to_dicts())):
            failures.append(f"alphabets.html.json: {name} differs from data/alphabets/{name}.parquet")
        n["numbers"] += df.height
    if doc["hp_classes_pr116"] != [{k: r[k] for k in ("alphabet", "residue", "class")} for r in hp116]:
        failures.append("alphabets.html.json: hp_classes_pr116 differs from PR 116's CSV")
    stats["alphabets_numbers_checked"] = n["numbers"]


def check_alphabets_words(snap, ev, kp, uk, reach, cl, hp116, failures, n) -> None:
    """The sentences, tooltips, labels and drawn positions of alphabets.html, rebuilt here."""
    fair = ev.filter(pl.col("mode") == "low_half")
    pf = fair.filter(pl.col("group") == "Pfam seed")
    above = pf.filter((pl.col("lo") > 1) & (pl.col("alphabet") != "protein20")).sort("rank")
    two = pf.filter(pl.col("classes") == 2).sort("rank")
    sd = pf.filter(pl.col("alphabet") == "sdm12").row(0, named=True)
    words = {0: "no alphabet keeps", 1: "one alphabet keeps", 2: "two alphabets keep", 3: "three alphabets keep",
             4: "four alphabets keep"}.get(above.height, f"{above.height} alphabets keep")
    below = ("Every two-class alphabet keeps less, with its interval below 1; the best of them"
             if (two["hi"] < 1).all() else "The best two-class alphabet")
    want = (f"In the fair test on Pfam seed pairs at 20% identity or less, {words} more evidence per aligned "
            f"position than the 20 amino acids, with the 95% interval above 1 ({', '.join(above['alphabet'])}). "
            f"sdm12 keeps {sd['ratio_to_20aa']:.3f} times as much (95% interval {sd['lo']:.3f}–{sd['hi']:.3f}). "
            f"{below}, {two['alphabet'][0]}, keeps {two['ratio_to_20aa'][0]:.2f} times as much. "
            f"n = {sd['n_items']:,} pair-halves from {sd['n_clusters']:,} Pfam families.")
    if snap["lead"].strip() != want:
        failures.append(f"alphabets.html lead sentence:\n    page {snap['lead'].strip()!r}\n    want {want!r}")
    n["numbers"] += 8
    # H/P sentence, from PR 116's CSV
    by_res = {}
    for r in hp116:
        by_res.setdefault(r["residue"], set()).add(r["class"])
    differ = [r for r in "AFILMVCGPWYDEHKNQRST" if len(by_res[r]) > 1]
    if f"agree on {20 - len(differ)} residues and differ only on {', '.join(differ[:-1])} and {differ[-1]}." not in snap["hpText"]:
        failures.append(f"alphabets.html H/P sentence {snap['hpText']!r}, differing residues {differ}")
    for row in snap["rows"]:
        a = row["alphabet"]
        nc = cl.filter(pl.col("alphabet") == a)["class_members"].n_unique()
        if row["nClasses"] != f"{nc} classes":
            failures.append(f"alphabets.html {a}: label {row['nClasses']!r}, {nc} classes in the table")
        for e in row["ev"]:
            r = fair.filter((pl.col("group") == e["group"]) & (pl.col("alphabet") == a)).row(0, named=True)
            kind = "Pfam families" if e["group"] == "Pfam seed" else "SCOPe folds"
            tip = (f"{a}: {r['ratio_to_20aa']:.3f} ({r['lo']:.3f}–{r['hi']:.3f}) times the evidence per aligned position "
                   f"of the 20 amino acids, 95% interval in brackets. Rank {r['rank']} of 19. "
                   f"n = {r['n_items']:,} pair-halves from {r['n_clusters']:,} {kind}.")
            if e["title"] != tip:
                failures.append(f"alphabets.html {a} {e['group']} tooltip {e['title']!r} != {tip!r}")
            # drawn: the dot inside its interval, and on the side of the dashed line its value is on
            ok = e["ciLo"] <= e["dot"] <= e["ciHi"] and (e["dot"] > e["one"]) == (r["ratio_to_20aa"] > 1) \
                and (e["ciLo"] > e["one"]) == (r["lo"] > 1) and (e["ciHi"] > e["one"]) == (r["hi"] > 1)
            if a != "protein20" and not ok:
                failures.append(f"alphabets.html {a} {e['group']}: drawn interval {e} does not match {r['lo']}-{r['hi']}")
            n["numbers"] += 2
        for k in row["kp"]:
            r = kp.filter((pl.col("alphabet") == a) & (pl.col("identity_bin") == k["bin"])).row(0, named=True)
            words_b = {"<20%": "20% or less", "20-30%": "20–30%", "30-40%": "30–40%", "40-60%": "40–60%", ">=60%": "over 60%"}
            tip = (f"{a}, Pfam seed pairs at {words_b[k['bin']]} identity: kappa {r['kappa']:.2f}. Same class at "
                   f"{100 * r['agreement']:.0f}% of aligned positions, {100 * r['chance']:.0f}% expected by chance. "
                   f"n = {r['n_pairs']:,} pairs.")
            if k["title"] != tip:
                failures.append(f"alphabets.html {a} kappa tooltip {k['title']!r} != {tip!r}")
            n["numbers"] += 3
        for s in row["uk"]:
            r = uk.filter((pl.col("alphabet") == a) & (pl.col("identity_bin") == s["bin"])).row(0, named=True)
            rr = reach.filter((pl.col("alphabet") == a) & (pl.col("identity_bin") == s["bin"]))
            share = lambda k: rr.filter(pl.col("k") == k)["share_reached"][0]  # noqa: E731
            lo, hi = r["k_use_lo"], r["k_use_hi"]
            rng = str(lo) if lo == hi else f"{lo}–{hi}"
            tip = (f"{a}, Pfam seed pairs at {s['bin'].replace('-', '–')} identity (n = {r['n_pairs']:,}): usable k {rng}. "
                   f"{100 * share(r['k_min']):.0f}% of pairs share a same-class run of k = {r['k_min']} or longer, "
                   f"{100 * share(hi):.0f}% one of k = {hi}. Chance-match range from notebook 274: k {r['k_min']}–{r['k_max']}.")
            if s["title"] != tip:
                failures.append(f"alphabets.html {a} usable k tooltip {s['title']!r} != {tip!r}")
            # drawn: the range line spans k_min-0.5 to k_max+0.5, so one k is (hi - lo) / (k_max - k_min + 1)
            unit = (s["rangeHi"] - s["rangeLo"]) / (r["k_max"] - r["k_min"] + 1)
            at = lambda k: s["rangeLo"] + (k - r["k_min"] + 0.5) * unit  # noqa: E731
            want_bars = [[at(lo - 0.5), at(hi + 0.5)]] + ([[at(lo - 0.5), at(r["k_most_hi"] + 0.5)]] if r["k_most_hi"] is not None else [])
            if len(s["bars"]) != len(want_bars) or any(abs(g - w) > 0.15 for gb, wb in zip(s["bars"], want_bars) for g, w in zip(gb, wb)):
                failures.append(f"alphabets.html {a} usable k {s['bin']}: bars {s['bars']} != {want_bars}")
            n["numbers"] += 4
    # row headers of the full table: group, homologs or not, pair selection, counts
    modes = {"low_half": "fair test: 20% identity or less, grouped on one half, measured on the other",
             "all": "all identities, whole pairs",
             "low_same": "20% identity or less, grouped and measured on the same positions (biased)"}
    want_heads = []
    for g in ev["group"].unique(maintain_order=True).to_list():
        for m in ("low_half", "all", "low_same"):
            e = ev.filter((pl.col("group") == g) & (pl.col("mode") == m))
            hom = "homologs" if e["homologs_by_truth"][0] else "not homologs by SCOPe"
            want_heads.append(f"{g} ({hom}){modes[m]}; {e['n_items'][0]:,} "
                              f"{'pair-halves' if m == 'low_half' else 'pairs'}, {e['n_clusters'][0]:,} clusters")
    if snap["fullHeads"] != want_heads:
        bad = [(g, w) for g, w in zip(snap["fullHeads"], want_heads) if g != w][:2]
        failures.append(f"alphabets.html full table row headers differ: {bad}")
    n["numbers"] += 2 * len(want_heads)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    ap.add_argument("--batch", type=int, default=150)
    a = ap.parse_args()

    # Expected values, from the long tables directly (not through build_site.py).
    long = pl.read_parquet(a.data / "results.parquet").filter(pl.col("experiment") == "qfo_region_midi_plus")
    truth = {"Pfam-A domain instances": "pfam", "Pfam-N domain instances": "pfamn", "Swiss-Prot features": "swissprot"}
    expect = {}
    for r in long.filter(pl.col("metric").is_in(["fmax", "auprc", "median_iou_tp"])).iter_rows(named=True):
        sp = "all" if r["target_db"] == "Pfam-A HMM library" else r["target_db"].removeprefix("QfO ")
        k = (truth[r["truth_source"]], r["dedup_or_raw"], sp, r["stratum_axis"], r["stratum_bin"], r["tool"], r["setting"])
        expect.setdefault(k, {"n_truth_instances": r["n_instances"]})[r["metric"]] = r["value"]
    for v in expect.values():
        v.setdefault("median_iou_tp", None)
    # Search cost per relative per 1,000 queries: cpu_h / n_target_dbs * 1000 / n_queries, for the
    # tool and setting's search step, shown for every relative that tool and setting searched.
    comp = pl.read_parquet(a.data / "compute.parquet").filter(pl.col("experiment") == "qfo_region_midi_plus")
    search = comp.filter((pl.col("setting") != "database build") & ~pl.col("step").str.contains("Index|Db|DB"))
    build = comp.filter((pl.col("setting") == "database build") | pl.col("step").str.contains("Index|Db|DB"))

    def db_per_db(tool, setting):
        # Build rows of this tool and setting, or of the tool family's shared build
        # ("mmseqs2" for both MMseqs2 searches): mean core-hours per database built.
        b = build.filter(pl.col("tool").is_in([tool, tool.split("_")[0]])
                         & pl.col("setting").is_in([setting, "database build"]))
        return b["cpu_h"].sum() / b["n_target_dbs"].sum() if b.height else None
    per_ts = {(r["tool"], r["setting"]): {
        "cpu_h_per_1000_queries": r["cpu_h"] / r["n_target_dbs"] * 1000 / r["n_queries"],
        "peak_rss_gb": r["peak_ram_gb"], "gpu_h_per_1000_queries": 0.0 if r["gpu_h"] is None else r["gpu_h"] / r["n_target_dbs"] * 1000 / r["n_queries"],
        "cpu_h_reported": None, "n_completed_runs": 1,
        "wall_h_per_1000_queries": r["wall_s"] / 3600 / r["n_target_dbs"] * 1000 / r["n_queries"],
        "db_cpu_h_per_db": db_per_db(r["tool"], r["setting"])} for r in search.iter_rows(named=True)}
    cexp = {(k[5], k[6], k[2]): per_ts[(k[5], k[6])] for k in expect if (k[5], k[6]) in per_ts}

    # Divergence time from human per relative, from the long table itself.
    mya_rows = long.select("target_db", "target_mya_to_human").unique()
    mya_want = {}
    for r in mya_rows.iter_rows(named=True):
        if r["target_db"] == "Pfam-A HMM library":
            continue
        sp = r["target_db"].removeprefix("QfO ")
        v = r["target_mya_to_human"]
        if sp in mya_want and mya_want[sp] != v:
            raise SystemExit(f"{sp}: more than one target_mya_to_human in the parquet")
        mya_want[sp] = v
    mya_shown = {sp: (None if v is None else int(round(v))) for sp, v in mya_want.items()}
    mya_allowed = {v for v in mya_shown.values() if v is not None}
    MYA_RE = re.compile(r"(\d[\d,]*)\s*(?:Mya|million years)")

    failures: list[str] = []
    stats = {"combinations": 0, "rows_checked": 0, "numbers_checked": 0,
             "missing_results_cells": 0, "missing_compute_cells": 0,
             "plot_renders": 0, "markers_checked": 0}
    httpd, base = serve(a.site)
    console_errors: list[str] = []
    with sync_playwright() as p:
        browser = p.chromium.launch(channel="chrome", headless=True)
        page = browser.new_page()
        page.on("console", lambda m: console_errors.append(f"{page.url}: {m.text}") if m.type == "error" else None)
        page.on("pageerror", lambda e: console_errors.append(f"{page.url}: {e}"))

        # ---- 1. leaderboard ----
        # Every row is in the page, the hidden "other settings" too, so every number is
        # checked; then the one visible row per tool and the sentence above the table.
        page.goto(base + "index.html")
        meta = page.evaluate("window.__lb.meta")
        order = {r["tool"] + "|" + r["setting"]: r["order"] for r in meta["rows"]}
        tools_in_order = list(dict.fromkeys(r["tool"] for r in meta["rows"]))
        all_combos = combos(meta)

        def cost_of(cw):
            return cw["cpu_h_per_1000_queries"] * 0.044625 + cw["gpu_h_per_1000_queries"] * 1.006

        def check_snap(snap, view):
            s = snap["s"]
            got_rows = {r["tool"] + "|" + r["setting"]: r for r in snap["rows"]}
            if set(got_rows) != set(order) or len(snap["rows"]) != len(order):
                failures.append(f"{s}: table rows {sorted(set(got_rows) ^ set(order))} differ from the tool-setting rows")
            vals = {}
            for row in snap["rows"]:
                stats["rows_checked"] += 1
                key = (s["truth"], s["calls"], s["species"], s["axis"], s["stratum"], row["tool"], row["setting"])
                want = expect.get(key)
                cells = row["cells"]
                tag = f"{key}"
                cw = cexp.get((row["tool"], row["setting"], s["species"]))
                if want is None:
                    bad = [c for c, v in cells.items() if c != "tool" and (v["missing"] != "results" or v["text"] != "× not run")]
                    if bad:
                        failures.append(f"{tag}: no parquet row but cells {bad} are not '× not run': {cells}")
                    stats["missing_results_cells"] += 1
                    continue
                shown = [cells[c]["text"] for c in ("fmax", "auprc", "iou", "n")]
                exp = [want.get("fmax"), want.get("auprc"), want.get("median_iou_tp")]
                # Three decimals on the page: each value within half of the last
                # digit (the page rounds halves up, Python's format rounds them to even).
                bad = any(abs(float(t) - e) > 0.0005 + 1e-9 for t, e in zip(shown[:2], exp[:2]))
                # A stored median IoU of 0 (no correct call longer than one residue)
                # is shown as words, never as a number or a mark: with Fmax above 0
                # some call was correct, so those correct calls are all on single
                # residues; with Fmax 0 no call was correct.
                if exp[2] is None or exp[2] == 0:
                    words = ("correct calls only on single residues" if exp[0] > 0 else "no correct calls")
                    bad |= cells["iou"]["none"] != "iou" or shown[2] != words
                    k_iou = "iou_single_residue_words" if exp[0] > 0 else "iou_no_correct_words"
                    stats[k_iou] = stats.get(k_iou, 0) + 1
                else:
                    bad |= abs(float(shown[2]) - exp[2]) > 0.0005 + 1e-9
                if bad or shown[3] != f"{want['n_truth_instances']:,}" or any(cells[c]["missing"] for c in ("fmax", "auprc", "iou", "n")):
                    failures.append(f"{tag}: page {shown} != parquet {exp + [want['n_truth_instances']]}")
                stats["numbers_checked"] += 4
                c_cost = None
                if cw is None:
                    # The search was scored here, so it ran: its cost is "not measured".
                    bad = [c for c in ("cost", "cpu", "ram") if cells[c]["missing"] != "compute" or cells[c]["text"] != "not measured"]
                    if bad:
                        failures.append(f"{tag}: no compute row but cells {[cells[c] for c in bad]}")
                    stats["missing_compute_cells"] += 3
                else:
                    c_cost = cost_of(cw)
                    cost_shown = float(cells["cost"]["text"].lstrip("$"))
                    if abs(cost_shown - c_cost) > 0.006 * c_cost + 1e-12:
                        failures.append(f"{tag}: cost {cost_shown} != {c_cost}")
                    cpu = float(cells["cpu"]["text"])
                    if abs(cpu - cw["cpu_h_per_1000_queries"]) > 0.006 * cw["cpu_h_per_1000_queries"] + 1e-9:
                        failures.append(f"{tag}: core-hours {cpu} != {cw['cpu_h_per_1000_queries']}")
                    if cells["ram"]["text"] != f"{cw['peak_rss_gb']:.1f}":
                        failures.append(f"{tag}: RAM {cells['ram']['text']} != {cw['peak_rss_gb']:.1f}")
                    stats["numbers_checked"] += 3
                vals[row["tool"] + "|" + row["setting"]] = (exp[0], exp[1], c_cost)
            # Sorted by Fmax, highest first, ties and missing values by row order.
            def best_of(keys, i):
                have = [k for k in keys if k in vals and vals[k][i] is not None]
                return min(have, key=lambda k: (-vals[k][i], order[k])) if have else None
            visible = [r["tool"] + "|" + r["setting"] for r in snap["rows"] if not r["hidden"]]
            if view == "best":
                want_vis = []
                for t in tools_in_order:
                    ks = [k for k in order if k.split("|")[0] == t]
                    want_vis.append(best_of(ks, 0) or min(ks, key=order.get))
                if sorted(visible) != sorted(want_vis):
                    failures.append(f"{s}: visible rows {sorted(set(visible) ^ set(want_vis))} are not each tool's best by Fmax")
                stats["best_rows_checked"] = stats.get("best_rows_checked", 0) + len(want_vis)
            elif sorted(visible) != sorted(order):
                failures.append(f"{s}: every-setting view hides rows")
            b = best_of(list(order), 0)
            near = [k for k in vals if b and vals[k][0] is not None and vals[k][0] >= vals[b][0] - 0.01 and vals[k][2] is not None]
            cheap = min(near, key=lambda k: (vals[k][2], order[k])) if near else None
            want_v = {"best": b or "", "cheap": cheap or "", "auprc": best_of(list(order), 1) or ""}
            if snap["verdict"] != want_v:
                failures.append(f"{s}: sentence above Table 1 names {snap['verdict']}, expected {want_v}")
            stats["sentences_checked"] = stats.get("sentences_checked", 0) + 1

        chips = {c["value"]: c["label"] for r in meta["chiprows"] if r["id"] == "species"
                 for g in r["groups"] for c in g["chips"]}
        for sp, want in mya_shown.items():
            name = meta["species"][sp]["name"]
            want_label = name if want is None else f"{name} {want}"
            if chips.get(sp) != want_label:
                failures.append(f"chip {sp}: label {chips.get(sp)!r}, parquet says {want_label!r}")
            if meta["species"][sp]["mya"] != want:
                failures.append(f"species {sp}: mya {meta['species'][sp]['mya']}, parquet says {want}")
            stats["mya_checked"] = stats.get("mya_checked", 0) + 1
        for i in range(0, len(all_combos), a.batch):
            for snap in page.evaluate(SNAP_TABLE, all_combos[i:i + a.batch]):
                stats["combinations"] += 1
                check_snap(snap, "best")
        # The every-setting view, with no subset chosen.
        every = [dict(c, view="all") for c in all_combos if c["axis"] == "all"]
        for snap in page.evaluate(SNAP_TABLE, every):
            check_snap(snap, "all")
        stats["every_setting_views"] = len(every)
        page.evaluate("window.__lb.set({view: 'best'})")

        # ---- 2. compare two columns ----
        page.goto(base + "cost.html")
        if page.evaluate("document.querySelector('#choices').open"):
            failures.append("cost.html: the filters panel starts open; the figure should come first")
        if page.evaluate("document.querySelector('#choices').compareDocumentPosition(document.querySelector('#plot')) & 4") == 0:
            failures.append("cost.html: the figure is not after the filters panel")
        score_field = {"fmax": "fmax", "auprc": "auprc", "iou": "median_iou_tp"}

        def axis_value(ax, k, cw):
            """Value of one row on one axis, straight from the parquet files, or None."""
            if ax in score_field:
                v = expect.get(k, {}).get(score_field[ax])
                return None if v is None or (ax == "iou" and v == 0) else v
            if k not in expect or cw is None:
                return None
            return {"cpu": cw["cpu_h_per_1000_queries"],
                    "cost": cw["cpu_h_per_1000_queries"] * 0.044625 + cw["gpu_h_per_1000_queries"] * 1.006,
                    "wall": cw["wall_h_per_1000_queries"], "ram": cw["peak_rss_gb"],
                    "db": cw["db_cpu_h_per_db"]}[ax]

        axis_ids = {a["id"] for a in meta["axes"]}
        no_stratum = [c for c in all_combos if c["axis"] == "all"]
        pairs = [(c, "cpu", y) for y in ("fmax", "auprc", "iou") for c in all_combos] + \
            [(c, x, y) for x, y in (("cost", "auprc"), ("ram", "iou"), ("fmax", "auprc"), ("wall", "fmax"), ("db", "fmax"))
             for c in no_stratum]
        for x, y in {(x, y) for _, x, y in pairs}:
            if x not in axis_ids or y not in axis_ids:
                failures.append(f"axis {x} or {y} not offered on cost.html")
        plot_combos = [dict(c, x=x, y=y) for c, x, y in pairs]
        stats["axis_pairs"] = sorted({f"{x}~{y}" for _, x, y in pairs})
        for i in range(0, len(plot_combos), a.batch):
            for snap in page.evaluate(SNAP_PLOT, plot_combos[i:i + a.batch]):
                stats["plot_renders"] += 1
                s = snap["s"]
                if (snap["state"]["x"], snap["state"]["y"]) != (s["x"], s["y"]):
                    failures.append(f"plot {s}: page drew {snap['state']['x']} against {snap['state']['y']}")
                want = {}
                for r in meta["rows"]:
                    k = (s["truth"], s["calls"], s["species"], s["axis"], s["stratum"], r["tool"], r["setting"])
                    cw = cexp.get((r["tool"], r["setting"], s["species"]))
                    vx, vy = axis_value(s["x"], k, cw), axis_value(s["y"], k, cw)
                    if vx is not None and vy is not None:
                        want[r["tool"] + "|" + r["setting"]] = (vx, vy)
                got = {m["tool"] + "|" + m["setting"]: (m["x"], m["y"]) for m in snap["plotted"]}
                dom = set(snap["marks"])
                if set(got) != set(want) or len(snap["plotted"]) != len(got) or dom != set(got):
                    failures.append(f"plot {s}: markers {sorted(set(got) ^ set(want))} differ from the rows with both values")
                for key in set(got) & set(want):
                    for gv, wv, name in zip(got[key], want[key], ("x", "y")):
                        if abs(gv - wv) > 1e-9 * max(1.0, abs(wv)):
                            failures.append(f"plot {s} {key}: {name} {gv} != parquet {wv}")
                    stats["marker_values_checked"] = stats.get("marker_values_checked", 0) + 2
                listed = set(snap["notplotted"])
                rows_all = {r["tool"] + "|" + r["setting"] for r in meta["rows"]}
                if listed | set(got) != rows_all or listed & set(got):
                    failures.append(f"plot {s}: not-plotted list {sorted(listed)} does not complete the markers")
                # The line: rows no other row beats on both (lower or equal x and higher or
                # equal y, one strictly), only when x is a cost and y a score.
                cost_x = s["x"] not in score_field
                front = {k for k, (x0, y0) in want.items()
                         if not any(k2 != k and x2 <= x0 and y2 >= y0 and (x2 < x0 or y2 > y0)
                                    for k2, (x2, y2) in want.items())} if cost_x else set()
                got_front = {m["tool"] + "|" + m["setting"] for m in snap["plotted"] if m["frontier"]}
                if got_front != front or snap["frontier"] != (len(front) > 1):
                    failures.append(f"plot {s}: line through {sorted(got_front)}, expected {sorted(front)}")
                stats["lines_drawn"] = stats.get("lines_drawn", 0) + int(snap["frontier"])
                stats["markers_checked"] += len(got)

        # ---- 2b. the URL hash: bad values fall back to the defaults, hashchange applies ----
        defaults = meta["defaults"]
        n_rows = len(meta["rows"])
        for bad in ("species=Zebrafish", "axis=identity&stratum=nope", "axis=nope&stratum=all",
                    "truth=pfam2", "metric=x&sort=evil&dir=up", "x=nope&y=cpu", "view=some", "r_cpu=-1&r_gpu=abc"):
            page.goto(base + "index.html#" + bad)
            got = page.evaluate("window.__lb.state()")
            shown = page.evaluate("document.querySelectorAll('#lb tbody tr').length")
            cost_cells = page.evaluate("[...document.querySelectorAll('#lb tbody td')].map(t => t.textContent).filter(t => t.startsWith('$-'))")
            shown = page.evaluate("document.querySelectorAll('#lb tbody tr').length")
            want = {k: defaults[k] for k in ("truth", "calls", "species", "axis", "stratum", "x", "y", "sort", "dir", "view", "r_cpu", "r_gpu")}
            if {k: got[k] for k in want} != want or shown != n_rows or cost_cells:
                failures.append(f"hash #{bad}: state {got}, {shown} rows, negative costs {cost_cells[:2]}")
            stats["bad_hashes_checked"] = stats.get("bad_hashes_checked", 0) + 1
        page.goto(base + "index.html")
        page.evaluate("location.hash = 'truth=swissprot&species=mouse&axis=feature_type&stratum=TRANSMEM&sort=auprc&dir=asc'")
        page.wait_for_function("window.__lb.state().species === 'mouse'", timeout=5000)
        got = page.evaluate("window.__lb.state()")
        caption = page.evaluate("document.getElementById('rankline').textContent")
        if (got["truth"], got["stratum"], got["sort"]) != ("swissprot", "TRANSMEM", "auprc"):
            failures.append(f"hashchange not applied: {got}")
        if caption != "Tools ranked by AUPRC, lowest first, best setting of each tool.":
            failures.append(f"caption after sorting by AUPRC ascending: {caption!r}")
        # Switching truth set to Pfam drops the Swiss-Prot-only feature-type stratum.
        page.evaluate("[...document.querySelectorAll('.chiprow[data-row=truth] .chip')].find(b => b.textContent === 'Pfam').click()")
        got = page.evaluate("window.__lb.state()")
        n_run = page.evaluate("[...document.querySelectorAll('#lb tbody tr')].filter(tr => !tr.querySelector('[data-missing=results]')).length")
        if (got["truth"], got["axis"], got["stratum"]) != ("pfam", "all", "all") or n_run == 0:
            failures.append(f"Swiss-Prot TRANSMEM -> Pfam kept the stratum: {got}, {n_run} rows with results")
        stats["hash_and_caption_checks"] = 3

        # ---- 2c. dark mode: plot markers take the dark kind colours of the legend ----
        # Same page (the test server answers one connection at a time): load in light mode,
        # then switch the colour scheme, which must redraw the plot without a reload.
        page.goto(base + "cost.html")
        page.emulate_media(color_scheme="dark")
        page.wait_for_timeout(200)
        dark = page
        col = dark.evaluate("""() => {
            const s = getComputedStyle(document.documentElement);
            const m = document.querySelector('#plot [data-tool="hmmer3_phmmer"]');
            const leg = [...document.querySelectorAll('.legend svg [stroke]')].map(e => e.getAttribute('stroke'));
            return {tok: s.getPropertyValue('--k-letters').trim(), marker: m && m.getAttribute('stroke'), legendHex: leg.filter(x => x.startsWith('#'))};
        }""")
        if not col["tok"] or col["marker"] != col["tok"] or col["legendHex"]:
            failures.append(f"dark mode cost.html: marker {col['marker']} vs token {col['tok']}, hex legend strokes {col['legendHex']}")
        stats["dark_mode_marker"] = col["marker"]
        page.emulate_media(color_scheme="light")

        # ---- 2d. the alphabet page ----
        if (a.site / "alphabets.html").exists():
            check_alphabets(page, base, a.data, a.site, failures, stats)

        # ---- 3. every page: links, errors, phone width ----
        pages = ["index.html", "cost.html", "coverage.html", "methods.html", "explainers.html", "tools/index.html"] + \
        (["pileup.html"] if (a.site / "pileup.html").exists() else []) + \
        (["alphabets.html"] if (a.site / "alphabets.html").exists() else []) + \
        [x for x in ("claims.html", "claims-schema.html") if (a.site / x).exists()] + \
            [f"tools/{t}.html" for t in sorted(meta["tools"])]
        files_ok = {str(q.relative_to(a.site)) for q in a.site.rglob("*") if q.is_file()}
        for width in (1280, 375):
            page.set_viewport_size({"width": width, "height": 900})
            for pg in pages:
                page.goto(base + pg)
                page.wait_for_load_state("networkidle")
                if not page.evaluate("document.body.textContent.includes(document.getElementById('source-commit').dataset.commit.slice(0,8))"):
                    failures.append(f"{pg}: source commit not printed")
                over = page.evaluate("document.documentElement.scrollWidth - window.innerWidth")
                if over > 0:
                    failures.append(f"{pg} at {width}px: page scrolls sideways by {over}px")
                if width == 1280:
                    text = page.evaluate("document.body.innerText")
                    for m in MYA_RE.finditer(text):
                        stats["mya_mentions_checked"] = stats.get("mya_mentions_checked", 0) + 1
                        if int(m.group(1).replace(",", "")) not in mya_allowed:
                            failures.append(f"{pg}: {m.group(0)!r} is no relative's divergence time "
                                            f"in the parquet ({sorted(mya_allowed)})")
                    for href in page.evaluate("[...document.querySelectorAll('a[href]')].map(a => a.getAttribute('href'))"):
                        u = urlparse(urljoin(base + pg, href))
                        if u.netloc != urlparse(base).netloc:
                            continue
                        path = u.path.lstrip("/") or "index.html"
                        if path not in files_ok:
                            failures.append(f"{pg}: broken link {href}")
        browser.close()
    httpd.shutdown()

    # ---- 4. every page has its JSON twin, naming the same source commit ----
    site_src = json.loads((a.site / "data/SOURCE.json").read_text())
    commit = site_src["commit"]
    for pg in pages:
        jp = a.site / (pg + ".json")
        if not jp.exists():
            failures.append(f"{pg}: no {pg}.json")
            continue
        doc = json.loads(jp.read_text())
        # The alphabet page's numbers come from their own commit.
        want_commit = site_src.get("alphabets_commit") if pg == "alphabets.html" else commit
        if doc["source"]["commit"] != want_commit or not doc["text"]:
            failures.append(f"{pg}.json: wrong commit or no text")
    stats["json_pages"] = len(pages)
    failures += [f"console: {e}" for e in console_errors]

    print(json.dumps(stats, indent=1))
    if failures:
        print(f"{len(failures)} FAILURES")
        for f in failures[:60]:
            print("  " + f)
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
