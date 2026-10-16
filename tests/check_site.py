#!/usr/bin/env python3
"""Headless check of the built site against data/results.parquet and data/compute.parquet.

1. Leaderboard: render every chip combination (truth set x dedup x relative x stratum).
   For each tool and setting row, the page must show the parquet's numbers when the
   parquet has a row, and "—" with "not run" when it has none. Every "—" on the page
   must be a cell marked missing, and every missing cell must be a missing row.
2. Compare two columns (Figure 1 on cost.html): render every combination for each of the
   three scores against core-hours, and every combination with no stratum chosen for five
   more axis pairs (cost against AUPRC, peak RAM against median IoU, Fmax against AUPRC,
   wall-clock time against Fmax, database build against Fmax). The plotted markers must be
   exactly the rows with a value on both axes, each marker's x and y must equal the value
   computed here from the parquet files, every other row must be listed as not plotted, and
   the line through the best trade-offs must be drawn only for a cost against a score.
   Core-hours, cost and wall-clock time are checked at the default 20,000 proteins.
2b. The cost calculator: for ten random (proteins, budget) pairs on random relatives, every
   scaled cost and the set of rows that fit the budget are recomputed from compute.parquet;
   markers are filled exactly for those rows and the caption names each tool that fits.
   cost.html.json's "calculator" block is recomputed the same way at the default inputs.
   Every tool card has a hardware line whose task counts match data/hardware.parquet.
3. Every page: no console errors, no broken internal links, and no horizontal page scroll
   at phone width (375 px).

Uses the installed Google Chrome through Playwright (channel="chrome").
"""

import argparse
import functools
import http.server
import itertools
import json
import random
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
  const rows = [...document.querySelectorAll('#lb tbody tr')].map((tr) => ({
    tool: tr.dataset.tool, setting: tr.dataset.setting,
    cells: [...tr.children].slice(2).map((td) => ({
      text: td.textContent.trim(), missing: td.dataset.missing || null, none: td.dataset.none || null,
      span: td.colSpan})),
  }));
  const dashes = (document.querySelector('#lb tbody').textContent.match(/—/g) || []).length;
  return {s, rows, dashes};
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

    failures: list[str] = []
    stats = {"combinations": 0, "rows_checked": 0, "numbers_checked": 0, "dashes": 0,
             "missing_results_cells": 0, "missing_compute_cells": 0, "missing_reported_cells": 0,
             "plot_renders": 0, "markers_checked": 0}
    httpd, base = serve(a.site)
    console_errors: list[str] = []
    with sync_playwright() as p:
        browser = p.chromium.launch(channel="chrome", headless=True)
        page = browser.new_page()
        page.on("console", lambda m: console_errors.append(f"{page.url}: {m.text}") if m.type == "error" else None)
        page.on("pageerror", lambda e: console_errors.append(f"{page.url}: {e}"))

        # ---- 1. leaderboard ----
        page.goto(base + "index.html")
        meta = page.evaluate("window.__lb.meta")
        all_combos = combos(meta)
        for i in range(0, len(all_combos), a.batch):
            for snap in page.evaluate(SNAP_TABLE, all_combos[i:i + a.batch]):
                stats["combinations"] += 1
                s = snap["s"]
                n_missing_cells = 0
                for row in snap["rows"]:
                    stats["rows_checked"] += 1
                    key = (s["truth"], s["calls"], s["species"], s["axis"], s["stratum"], row["tool"], row["setting"])
                    want = expect.get(key)
                    cells = row["cells"]
                    tag = f"{key}"
                    if want is None:
                        c0 = cells[0]
                        if not (c0["missing"] == "results" and c0["span"] == 4 and c0["text"] == "— not run"):
                            failures.append(f"{tag}: no parquet row but page shows {c0}")
                        stats["missing_results_cells"] += 1
                        rest = cells[1:]
                    else:
                        shown = [c["text"] for c in cells[:4]]
                        exp = [want.get("fmax"), want.get("auprc"), want.get("median_iou_tp")]
                        # Three decimals on the page: each value within half of the last
                        # digit (the page rounds halves up, Python's format rounds them to even).
                        bad = any(abs(float(t) - e) > 0.0005 + 1e-9 for t, e in zip(shown[:2], exp[:2]))
                        # A stored median IoU of 0 (no correct call longer than one residue)
                        # is shown as words, never as a number or a dash.
                        if exp[2] is None or exp[2] == 0:
                            bad |= cells[2]["none"] != "iou" or shown[2] != "no correct calls"
                            stats["iou_as_words"] = stats.get("iou_as_words", 0) + 1
                        else:
                            bad |= abs(float(shown[2]) - exp[2]) > 0.0005 + 1e-9
                        if bad or shown[3] != f"{want['n_truth_instances']:,}" or any(c["missing"] for c in cells[:4]):
                            failures.append(f"{tag}: page {shown} != parquet {exp + [want['n_truth_instances']]}")
                        stats["numbers_checked"] += 4
                        rest = cells[4:]
                    cw = cexp.get((row["tool"], row["setting"], s["species"]))
                    if cw is None:
                        kinds = [c["missing"] for c in rest]
                        if kinds != ["compute", "reported", "compute", "compute", "compute"]:
                            failures.append(f"{tag}: no compute row but cells {rest}")
                        stats["missing_compute_cells"] += 4
                        stats["missing_reported_cells"] += 1
                    else:
                        cpu = rest[0]["text"].split()[0]
                        if cw["n_completed_runs"] > 1 and f"{cw['n_completed_runs']} runs" not in rest[0]["text"]:
                            failures.append(f"{tag}: {cw['n_completed_runs']} runs but no range shown")
                        gpu_text = rest[3]["text"]
                        if cw["gpu_h_per_1000_queries"] > 0:
                            g = float(gpu_text.split()[0])
                            if abs(g - cw["gpu_h_per_1000_queries"]) > 0.006 * cw["gpu_h_per_1000_queries"]:
                                failures.append(f"{tag}: GPU {gpu_text} != {cw['gpu_h_per_1000_queries']}")
                        elif gpu_text != "none":
                            failures.append(f"{tag}: GPU {gpu_text} but parquet has none")
                        cost_shown = float(rest[4]["text"].lstrip("$"))
                        cost_want = cw["cpu_h_per_1000_queries"] * 0.044625 + cw["gpu_h_per_1000_queries"] * 1.006
                        if abs(cost_shown - cost_want) > 0.006 * cost_want + 1e-12:
                            failures.append(f"{tag}: cost {cost_shown} != {cost_want}")
                        ram = rest[2]["text"]
                        if abs(float(cpu) - cw["cpu_h_per_1000_queries"]) > 0.006 * cw["cpu_h_per_1000_queries"] + 1e-9:
                            failures.append(f"{tag}: CPU-h {cpu} != {cw['cpu_h_per_1000_queries']}")
                        if ram != f"{cw['peak_rss_gb']:.1f}":
                            failures.append(f"{tag}: RAM {ram} != {cw['peak_rss_gb']:.1f}")
                        if cw["cpu_h_reported"] is None and rest[1]["missing"] != "reported":
                            failures.append(f"{tag}: reported CPU-h shown though parquet has none")
                        if rest[1]["missing"]:
                            stats["missing_reported_cells"] += 1
                        stats["numbers_checked"] += 2
                    n_missing_cells += sum(1 for c in cells if c["missing"])
                stats["dashes"] += snap["dashes"]
                if snap["dashes"] != n_missing_cells:
                    failures.append(f"{s}: {snap['dashes']} dashes on the page but {n_missing_cells} missing cells")

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
            # Core-hours, cost and wall-clock time are drawn for the calculator's n proteins
            # (default 20,000): the value per 1,000 queries times n / 1000.
            f = meta["defaults"]["n"] / 1000
            return {"cpu": cw["cpu_h_per_1000_queries"] * f,
                    "cost": (cw["cpu_h_per_1000_queries"] * 0.044625 + cw["gpu_h_per_1000_queries"] * 1.006) * f,
                    "wall": cw["wall_h_per_1000_queries"] * f, "ram": cw["peak_rss_gb"],
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

        # ---- 2b. the cost calculator: ten random (proteins, budget) pairs ----
        # For each pair, on a random relative, every scaled number and the set of tool and
        # setting rows that fit is recomputed here from compute.parquet: cost for N proteins =
        # (core-hours per 1,000 x r_cpu + GPU-hours per 1,000 x r_gpu) x N / 1000; it fits
        # when that is at most the budget. Markers must be filled exactly for those rows.
        def db_gpu_per_db(tool, setting):
            b = build.filter(pl.col("tool").is_in([tool, tool.split("_")[0]])
                             & pl.col("setting").is_in([setting, "database build"]))
            g = b["gpu_h"].drop_nulls()
            return g.sum() / b["n_target_dbs"].sum() if b.height and g.len() else None

        rng = random.Random(20261016)
        species = sorted({k[2] for k in expect if k[2] != "all"})
        calc_pairs = [(int(10 ** rng.uniform(2, 5.5)), round(10 ** rng.uniform(-2, 3), 2), rng.choice(species))
                      for _ in range(10)]
        # Two more at the edge: a budget equal to one row's cost (it fits, since the rule is
        # "at most") and 5% under it (it does not). Random budgets rarely land that close.
        edge = cexp[("hmmer3_phmmer", "default", "zebrafish")]
        edge_usd = (edge["cpu_h_per_1000_queries"] * 0.044625 + edge["gpu_h_per_1000_queries"] * 1.006) * (20000 / 1000)
        calc_pairs += [(20000, edge_usd, "zebrafish"), (20000, edge_usd * 0.95, "zebrafish")]
        stats["calculator_pairs"] = [list(p) for p in calc_pairs]
        for n_prot, budget, sp in calc_pairs:
            s = {"truth": "pfam", "calls": "raw", "species": sp, "axis": "all", "stratum": "all",
                 "x": "cost", "y": "fmax", "n": n_prot, "budget": budget, "r_cpu": 0.044625, "r_gpu": 1.006}
            snap = page.evaluate("""(s) => { window.__lb.set(s);
              return {calc: window.__lb.calc(), plotted: window.__lb.plotted(),
                      fitline: document.getElementById('fitline').textContent,
                      budget: !!document.querySelector('#plot rect.overbudget'),
                      filled: [...document.querySelectorAll('#plot [data-tool]')].map((m) =>
                        [m.dataset.tool + '|' + m.dataset.setting, m.getAttribute('fill') !== 'var(--surface)'])}; }""", s)
            f = n_prot / 1000
            want_fit, rows_seen = set(), 0
            got_rows = {r["tool"] + "|" + r["setting"]: r for r in snap["calc"]["rows"]}
            for r in meta["rows"]:
                k = ("pfam", "raw", sp, "all", "all", r["tool"], r["setting"])
                if k not in expect:
                    if r["tool"] + "|" + r["setting"] in got_rows:
                        failures.append(f"calc {s}: {r['tool']} {r['setting']} has no result but is in the calculator")
                    continue
                key = r["tool"] + "|" + r["setting"]
                got = got_rows.get(key)
                cw = cexp.get((r["tool"], r["setting"], sp))
                if got is None:
                    failures.append(f"calc {s}: {key} missing from the calculator")
                    continue
                rows_seen += 1
                if cw is None:
                    if got["measured"]:
                        failures.append(f"calc {s}: {key} has no compute row but the page has a cost")
                    continue
                usd = (cw["cpu_h_per_1000_queries"] * 0.044625 + cw["gpu_h_per_1000_queries"] * 1.006) * f
                want = {"core_h": cw["cpu_h_per_1000_queries"] * f, "usd": usd,
                        "wall_h": cw["wall_h_per_1000_queries"] * f, "gpu_h": cw["gpu_h_per_1000_queries"] * f,
                        "db_core_h_per_db": cw["db_cpu_h_per_db"],
                        "db_gpu_h_per_db": db_gpu_per_db(r["tool"], r["setting"])}
                for name, wv in want.items():
                    gv = got.get(name)
                    if (wv is None) != (gv is None) or (wv is not None and abs(gv - wv) > 1e-9 * max(1.0, abs(wv))):
                        failures.append(f"calc {s}: {key} {name} {gv} != parquet {wv}")
                    stats["calculator_values_checked"] = stats.get("calculator_values_checked", 0) + 1
                if usd <= budget:
                    want_fit.add(key)
            if set(snap["calc"]["fit"]) != want_fit:
                failures.append(f"calc {s}: fits {sorted(snap['calc']['fit'])}, expected {sorted(want_fit)}")
            filled = {k for k, fl in snap["filled"] if fl}
            if filled != want_fit & {k for k, _ in snap["filled"]}:
                failures.append(f"calc {s}: filled markers {sorted(filled)} != rows that fit {sorted(want_fit)}")
            names = {meta["tools"][k.split("|")[0]]["name"] for k in want_fit}
            fit_part = snap["fitline"].split("Over the budget")[0]
            for nm in names:
                if nm not in fit_part:
                    failures.append(f"calc {s}: caption does not name {nm} among the tools that fit")
            if not snap["budget"]:
                failures.append(f"calc {s}: no shaded over-budget area on the cost axis")
            stats["calculator_rows_checked"] = stats.get("calculator_rows_checked", 0) + rows_seen

        # By hand (PR 7's number x 20): phmmer, Pfam, raw calls, zebrafish, 20,000 proteins.
        snap = page.evaluate("""() => { window.__lb.set({truth: 'pfam', calls: 'raw', species: 'zebrafish',
            axis: 'all', stratum: 'all', n: 20000, budget: 10}); return window.__lb.calc(); }""")
        ph = [r for r in snap["rows"] if r["tool"] == "hmmer3_phmmer"]
        if not ph or f"{ph[0]['core_h']:.2f}" != "8.33":
            failures.append(f"calc: phmmer at 20,000 proteins is {ph and ph[0]['core_h']} core-hours, not 8.33")
        stats["phmmer_core_h_at_20000"] = ph[0]["core_h"] if ph else None

        # ---- 3. every page: links, errors, phone width ----
        pages = ["index.html", "cost.html", "coverage.html", "methods.html", "tools/index.html"] + \
        (["pileup.html"] if (a.site / "pileup.html").exists() else []) + \
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
    commit = json.loads((a.site / "data/SOURCE.json").read_text())["commit"]
    for pg in pages:
        jp = a.site / (pg + ".json")
        if not jp.exists():
            failures.append(f"{pg}: no {pg}.json")
            continue
        doc = json.loads(jp.read_text())
        if doc["source"]["commit"] != commit or not doc["text"]:
            failures.append(f"{pg}.json: wrong commit or no text")
    stats["json_pages"] = len(pages)

    # ---- 5. the calculator's JSON twin, at the default inputs ----
    twin = json.loads((a.site / "cost.html.json").read_text())["data"]["calculator"]
    d = meta["defaults"]
    if {k: twin["inputs"][k] for k in ("n", "budget", "r_cpu", "r_gpu")} != {k: d[k] for k in ("n", "budget", "r_cpu", "r_gpu")}:
        failures.append(f"cost.html.json: calculator inputs {twin['inputs']} are not the page defaults")
    f = d["n"] / 1000
    want_fit = set()
    for r in twin["rows"]:
        cw = cexp.get((r["tool"], r["setting"], d["species"]))
        if cw is None:
            if r["measured"]:
                failures.append(f"cost.html.json: {r['tool']} {r['setting']} has a cost but no compute row")
            continue
        usd = (cw["cpu_h_per_1000_queries"] * d["r_cpu"] + cw["gpu_h_per_1000_queries"] * d["r_gpu"]) * f
        for name, wv in (("usd", usd), ("core_h", cw["cpu_h_per_1000_queries"] * f), ("wall_h", cw["wall_h_per_1000_queries"] * f)):
            if abs(r[name] - wv) > 1e-9 * max(1.0, wv):
                failures.append(f"cost.html.json: {r['tool']} {r['setting']} {name} {r[name]} != {wv}")
        if usd <= d["budget"]:
            want_fit.add(f"{r['tool']}|{r['setting']}")
    if set(twin["fit"]) != want_fit:
        failures.append(f"cost.html.json: fit {twin['fit']} != {sorted(want_fit)}")
    stats["json_twin_rows"] = len(twin["rows"])
    # The twin's rule at the budget's edge: build_site.calculator() must also keep "at most".
    sys.path.insert(0, str(ROOT / "scripts"))
    from build_site import calculator  # noqa: E402
    cj = json.loads((a.site / "cost.html.json").read_text())["data"]
    for b, fits in ((edge_usd, True), (edge_usd * 0.95, False)):
        got = calculator(cj["results"], cj["compute"], cj["rows"], dict(d, species="zebrafish", budget=b))
        if ("hmmer3_phmmer|default" in got["fit"]) != fits:
            failures.append(f"build_site.calculator(): phmmer at budget {b} fits={not fits}, expected {fits}")

    # ---- 6. hardware line on every tool card, recounted from data/hardware.parquet ----
    hwp = a.data / "hardware.parquet"
    for t in sorted(meta["tools"]):
        card = json.loads((a.site / f"tools/{t}.json").read_text())
        md = (a.site / f"tools/{t}.md").read_text()
        if "Hardware of the search" not in md or "Hardware of the database build" not in md:
            failures.append(f"tools/{t}.md: no hardware line")
        if not hwp.exists():
            continue
        hw = pl.read_parquet(hwp)
        is_build = (pl.col("setting") == "database build") | pl.col("step").str.contains("Index|Db|DB")
        q = hw.filter(pl.col("tool").is_in([t, t.split("_")[0]]) & ~is_build
                      & pl.col("setting").is_in([card["best_setting"] or "", "database build"]))
        h = card["hardware"]["search"]
        if q.height == 0:
            if h["recorded"]:
                failures.append(f"tools/{t}.json: hardware recorded but no task in hardware.parquet")
            continue
        cpus = {r["cpu_model"]: r["len"] for r in q.group_by("cpu_model").len().iter_rows(named=True)}
        if h["n_tasks"] != q.height or {c["value"]: c["n_tasks"] for c in h["cpu_models"]} != cpus:
            failures.append(f"tools/{t}.json: hardware {h} != parquet {cpus}")
        stats["hardware_cards_checked"] = stats.get("hardware_cards_checked", 0) + 1
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
