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
                        # is shown as words, never as a number or a dash: with Fmax above 0
                        # some call was correct, so those correct calls are all on single
                        # residues; with Fmax 0 no call was correct.
                        if exp[2] is None or exp[2] == 0:
                            words = ("correct calls only on single residues" if exp[0] > 0
                                     else "no correct calls")
                            bad |= cells[2]["none"] != "iou" or shown[2] != words
                            k_iou = "iou_single_residue_words" if exp[0] > 0 else "iou_no_correct_words"
                            stats[k_iou] = stats.get(k_iou, 0) + 1
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
                    "truth=pfam2", "metric=x&sort=evil&dir=up", "x=nope&y=cpu", "r_cpu=-1&r_gpu=abc"):
            page.goto(base + "index.html#" + bad)
            got = page.evaluate("window.__lb.state()")
            shown = page.evaluate("document.querySelectorAll('#lb tbody tr').length")
            cost_cells = page.evaluate("[...document.querySelectorAll('#lb tbody td')].map(t => t.textContent).filter(t => t.startsWith('$-'))")
            want = {k: defaults[k] for k in ("truth", "calls", "species", "axis", "stratum", "x", "y", "sort", "dir", "r_cpu", "r_gpu")}
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
        if caption != "Tools ranked by AUPRC, lowest first":
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
