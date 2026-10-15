#!/usr/bin/env python3
"""Headless checks of the region pile-up page against the raw tables.

1. Three rows picked by rule from data/pileup/*.parquet (not from pileup.json): a kmerseek
   call that covers part of a domain, a domain near the protein's start (so the left flank
   runs off the protein), and a domain that runs past the protein's end. For each, the
   expected state of every column (covered, not covered, past the end) is computed here from
   the raw call and truth coordinates with code separate from scripts/pileup.py, and compared
   with the pixels drawn on the full-size canvas.
2. For every panel, the curve's body height (mean coverage inside the domain, from the
   page's own readouts()) against per-residue recall computed from the raw call intervals
   without the matrix (data/pileup_check.json).
3. Every sort order renders; under "coverage", every domain with no call inside it sits
   below every domain with one.
"""

import functools
import http.server
import json
import socketserver
import sys
import threading
from pathlib import Path

import polars as pl
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *a):
        pass


def expected_states(s, e, plen, calls):
    """'on' / 'off' / 'null' for the 300 columns, written independently of pileup.py."""
    L = e - s + 1
    out = []
    for j in range(300):
        if j < 100:
            res = s - 100 + j
        elif j < 200:
            res = s + int(((j - 100) + 0.5) * L / 100)
        else:
            res = e + 1 + (j - 200)
        if res < 1 or res > plen:
            out.append("null")
        else:
            out.append("on" if any(a <= res <= b for a, b in calls) else "off")
    return out


def main() -> int:
    reg = pl.read_parquet(ROOT / "data/pileup/regions.parquet").sort("row_id")
    calls = pl.read_parquet(ROOT / "data/pileup/calls.parquet")
    check = json.loads((ROOT / "docs/data/pileup_check.json").read_text())
    picks = []
    km = calls.filter((pl.col("tool") == "kmerseek") & (pl.col("setting") == "gbmr4_k21_lcFalse")).join(reg, on="row_id")
    r1 = km.filter((pl.col("call_start") > pl.col("true_start")) & (pl.col("call_end") < pl.col("true_end"))).sort("row_id").row(0, named=True)
    picks.append(("kmerseek", "gbmr4_k21_lcFalse", r1["row_id"], "kmerseek call inside a domain"))
    r2 = reg.filter(pl.col("true_start") < 40).sort("row_id").row(0, named=True)
    picks.append(("hmmer3_phmmer", "default", r2["row_id"], "domain near the protein's start"))
    r3 = reg.filter(pl.col("true_end") > pl.col("protein_length")).filter(pl.col("true_start") <= pl.col("protein_length")).sort("row_id").row(0, named=True)
    picks.append(("hmmer3_phmmer", "default", r3["row_id"], "domain running past the protein's end"))

    failures, report = [], []
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(ROOT / "docs")))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    errors = []
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        page = b.new_page(viewport={"width": 1280, "height": 900}, color_scheme="light")
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
        page.goto(f"http://127.0.0.1:{httpd.server_address[1]}/pileup.html")
        page.wait_for_function("window.__pileup && window.__pileup.ready", timeout=60000)
        panels = page.evaluate("window.__pileup.panels")

        # ---- 1. three rows by hand ----
        for tool, setting, row_id, why in picks:
            k = panels.index(f"{tool}|{setting}")
            page.evaluate(f"window.__pileup.setFocus({k})")
            r = reg.filter(pl.col("row_id") == row_id).row(0, named=True)
            raw = calls.filter((pl.col("tool") == tool) & (pl.col("setting") == setting) & (pl.col("row_id") == row_id))
            ivs = [(a, b_) for a, b_ in raw.select("call_start", "call_end").iter_rows()]
            want = expected_states(r["true_start"], r["true_end"], r["protein_length"], ivs)
            # Read the drawn pixels of this row: canvas x = 10 + 2*j, y = its position in the order.
            got = page.evaluate(
                """([rowIndex]) => {
                  const ord = window.__pileup.order, y = ord.indexOf(rowIndex);
                  const c = document.getElementById('pu-fheat').getContext('2d');
                  const px = c.getImageData(0, y, c.canvas.width, 1).data;
                  const out = []; for (let j = 0; j < 300; j++) { const o = (10 + 2 * j) * 4; out.push([px[o], px[o+1], px[o+2]]); }
                  return out;
                }""", [reg.with_row_index("i").filter(pl.col("row_id") == row_id)["i"][0]])
            colors = page.evaluate("""() => { const s = getComputedStyle(document.documentElement);
                return {null: s.getPropertyValue('--null').trim(), surface: s.getPropertyValue('--surface').trim()}; }""")
            hexrgb = lambda h: [int(h[i:i + 2], 16) for i in (1, 3, 5)]  # noqa: E731
            nul, off = hexrgb(colors["null"]), hexrgb(colors["surface"])
            drawn = ["null" if px == nul else "off" if px == off else "on" for px in got]
            bad = [j for j in range(300) if drawn[j] != want[j]]
            report.append({
                "why": why, "tool": tool, "setting": setting, "accession": r["accession"],
                "pfam": r["pfam_id"], "true": [r["true_start"], r["true_end"]],
                "protein_length": r["protein_length"], "calls": ivs[:5], "n_calls": len(ivs),
                "columns_on": want.count("on"), "columns_null": want.count("null"),
                "columns_disagreeing": len(bad),
            })
            if bad:
                failures.append(f"{why}: {len(bad)} columns differ, first {bad[:5]}")

        # ---- 2. body height against per-residue recall ----
        diffs = []
        for k, key in enumerate(panels):
            ro = page.evaluate(f"window.__pileup.readoutsOf({k})")
            want = check[key]["per_residue_recall_mean_over_regions"]
            diffs.append((key, ro["bodyMean"], want, ro["bodyMean"] - want))
            if abs(ro["bodyMean"] - want) > 0.01:
                failures.append(f"{key}: curve body {ro['bodyMean']:.4f} vs per-residue recall {want:.4f}")

        # ---- 3. sorts ----
        for srt in ("coverage", "identity", "length", "plddt"):
            page.evaluate(f"window.__pileup.setSort('{srt}')")
        page.evaluate("window.__pileup.setFocus(0); window.__pileup.setSort('coverage')")
        sink = page.evaluate("""() => {
            const P = window.__pileup, m = P.cov[0], ord = P.order;
            const hit = ord.map((i) => { let c = 0; for (let j = 100; j < 200; j++) if (m[i*300+j] === 1) c++; return c > 0; });
            const firstMiss = hit.indexOf(false); return firstMiss < 0 || hit.slice(firstMiss).every((h) => !h);
        }""")
        if not sink:
            failures.append("coverage sort: a covered domain sits below a missed one")
        b.close()
    httpd.shutdown()
    failures += [f"console: {e}" for e in errors]
    print(json.dumps({"hand_checked_rows": report}, indent=1, default=str))
    print("curve body vs per-residue recall (largest differences):")
    for key, got, want, d in sorted(diffs, key=lambda x: -abs(x[3]))[:6]:
        print(f"  {key:40s} curve {got:.4f}  recall {want:.4f}  diff {d:+.2e}")
    print(f"max |diff| {max(abs(d[3]) for d in diffs):.2e} over {len(diffs)} panels")
    if failures:
        print(f"{len(failures)} FAILURES")
        for f in failures:
            print("  " + f)
        return 1
    print("all pile-up checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
