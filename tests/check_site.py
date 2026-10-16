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
   at phone width (375 px).

Uses the installed Google Chrome through Playwright (channel="chrome").
"""

import argparse
import functools
import http.server
import itertools
import json
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
        page.goto(base + "leaderboard.html")
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
            page.goto(base + "leaderboard.html#" + bad)
            got = page.evaluate("window.__lb.state()")
            shown = page.evaluate("document.querySelectorAll('#lb tbody tr').length")
            cost_cells = page.evaluate("[...document.querySelectorAll('#lb tbody td')].map(t => t.textContent).filter(t => t.startsWith('$-'))")
            shown = page.evaluate("document.querySelectorAll('#lb tbody tr').length")
            want = {k: defaults[k] for k in ("truth", "calls", "species", "axis", "stratum", "x", "y", "sort", "dir", "view", "r_cpu", "r_gpu")}
            if {k: got[k] for k in want} != want or shown != n_rows or cost_cells:
                failures.append(f"hash #{bad}: state {got}, {shown} rows, negative costs {cost_cells[:2]}")
            stats["bad_hashes_checked"] = stats.get("bad_hashes_checked", 0) + 1
        page.goto(base + "leaderboard.html")
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

        # ---- 3. every page: links, errors, phone width ----
        pages = ["index.html", "cost.html", "coverage.html", "methods.html", "explainers.html", "tools/index.html"] + \
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
