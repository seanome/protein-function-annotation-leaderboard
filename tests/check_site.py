#!/usr/bin/env python3
"""Headless check of the built site against data/results.parquet and data/compute.parquet.

1. Leaderboard: render every chip combination (truth set x dedup x relative x stratum).
   For each tool and setting row, the page must show the parquet's numbers when the
   parquet has a row, and "—" with "not run" when it has none. Every "—" on the page
   must be a cell marked missing, and every missing cell must be a missing row.
2. Score against cost: render every combination for each of the three scores. The plotted
   markers must be exactly the rows with both a score and a measured cost.
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
  return {s, marks, notplotted: document.getElementById('notplotted').textContent};
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
    per_ts = {(r["tool"], r["setting"]): {
        "cpu_h_per_1000_queries": r["cpu_h"] / r["n_target_dbs"] * 1000 / r["n_queries"],
        "peak_rss_gb": r["peak_ram_gb"], "gpu_h_per_1000_queries": 0.0 if r["gpu_h"] is None else r["gpu_h"] / r["n_target_dbs"] * 1000 / r["n_queries"],
        "cpu_h_reported": None, "n_completed_runs": 1} for r in search.iter_rows(named=True)}
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

        # ---- 2. score against cost ----
        page.goto(base + "cost.html")
        plot_combos = [dict(c, metric=m["id"]) for c in all_combos for m in meta["metrics"]]
        field = {"fmax": "fmax", "auprc": "auprc", "iou": "median_iou_tp"}
        for i in range(0, len(plot_combos), a.batch):
            for snap in page.evaluate(SNAP_PLOT, plot_combos[i:i + a.batch]):
                stats["plot_renders"] += 1
                s = snap["s"]
                want = set()
                for r in meta["rows"]:
                    k = (s["truth"], s["calls"], s["species"], s["axis"], s["stratum"], r["tool"], r["setting"])
                    if k in expect and (r["tool"], r["setting"], s["species"]) in cexp \
                            and expect[k].get(field[s["metric"]]) not in (None, 0.0 if s["metric"] == "iou" else None):
                        want.add(r["tool"] + "|" + r["setting"])
                got = set(snap["marks"])
                stats["markers_checked"] += len(got)
                if got != want or len(snap["marks"]) != len(got):
                    failures.append(f"plot {s}: markers {sorted(got ^ want)} differ from the rows with score and cost")

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
