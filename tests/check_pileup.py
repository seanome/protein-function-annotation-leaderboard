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

import argparse
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


def kmerseek_setting(data: Path) -> str:
    """The kmerseek setting the pile-up page draws in Figure 4, by the page's rule: highest
    Fmax against zebrafish, Pfam truth set, raw calls, all proteins (ties to the setting
    name). Read from results.parquet directly, not through build_site.py."""
    long = pl.read_parquet(data / "results.parquet")
    q = long.filter((pl.col("experiment") == "qfo_region_midi_plus") & (pl.col("tool") == "kmerseek")
                    & (pl.col("metric") == "fmax") & (pl.col("truth_source") == "Pfam-A domain instances")
                    & (pl.col("dedup_or_raw") == "raw") & (pl.col("target_db") == "QfO zebrafish")
                    & (pl.col("stratum_axis") == "all"))
    return q.sort(["value", "setting"], descending=[True, False])["setting"][0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    a = ap.parse_args()
    reg = pl.read_parquet(a.data / "pileup/regions.parquet").sort("row_id")
    calls = pl.read_parquet(a.data / "pileup/calls.parquet")
    check = json.loads((a.site / "data/pileup_check.json").read_text())
    picks = []
    km_setting = kmerseek_setting(a.data)
    km = calls.filter((pl.col("tool") == "kmerseek") & (pl.col("setting") == km_setting)).join(reg, on="row_id")
    r1 = km.filter((pl.col("call_start") > pl.col("true_start")) & (pl.col("call_end") < pl.col("true_end"))).sort("row_id").row(0, named=True)
    picks.append(("kmerseek", km_setting, r1["row_id"], "kmerseek call inside a domain"))
    r2 = reg.filter(pl.col("true_start") < 40).sort("row_id").row(0, named=True)
    picks.append(("hmmer3_phmmer", "default", r2["row_id"], "domain near the protein's start"))
    r3 = reg.filter(pl.col("true_end") > pl.col("protein_length")).filter(pl.col("true_start") <= pl.col("protein_length")).sort("row_id").row(0, named=True)
    picks.append(("hmmer3_phmmer", "default", r3["row_id"], "domain running past the protein's end"))

    failures, report = [], []
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(a.site)))
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

        # ---- 2b. the page's readouts against pileup.html.json (Python copy) ----
        pj = json.loads((a.site / "pileup.html.json").read_text())["data"]["panels"]
        n_cmp = 0
        for k, key in enumerate(panels):
            for fam, name in ((False, "all_domains"), (True, "one_domain_per_family")):
                js = page.evaluate(f"window.__pileup.readoutsOf({k}, {str(fam).lower()})")
                py = pj[key][name]
                for jk, pk in (("missedFrac", "missed_fraction"), ("bodyMean", "mean_coverage_inside"),
                               ("flankMean", "mean_coverage_flanks")):
                    n_cmp += 1
                    if abs(js[jk] - py[pk]) > 1e-9:
                        failures.append(f"{key} {name} {pk}: page {js[jk]} vs json {py[pk]}")
                for jk, pk in (("startAbs", "start_0.8_to_0.2"), ("endAbs", "end_0.8_to_0.2"),
                               ("startRel", "start_80_to_20_percent_of_plateau"),
                               ("endRel", "end_80_to_20_percent_of_plateau")):
                    n_cmp += 1
                    if js[jk] != py[pk]:
                        failures.append(f"{key} {name} {pk}: page {js[jk]!r} vs json {py[pk]!r}")
        report.append({"page_vs_json_readouts_compared": n_cmp})

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

        # ---- 4. hover: the merged list covers every call the count counts ----
        # First phmmer row (by row_id) whose only calls sit in the flank, off the domain.
        ph = calls.filter((pl.col("tool") == "hmmer3_phmmer") & (pl.col("setting") == "default")).join(reg, on="row_id")
        ov = (pl.col("call_end") >= pl.col("true_start")) & (pl.col("call_start") <= pl.col("true_end"))
        off_only = (ph.group_by("row_id").agg(ov.any().alias("any_ov"), pl.len().alias("n"),
                                              pl.col("call_start").min(), pl.col("call_end").max())
                    .filter(~pl.col("any_ov") & (pl.col("n") == 1)).sort("row_id").row(0, named=True))
        idx = reg.with_row_index("i").filter(pl.col("row_id") == off_only["row_id"])["i"][0]
        page.evaluate(f"window.__pileup.setFocus({panels.index('hmmer3_phmmer|default')})")
        hover = page.evaluate(f"window.__pileup.hoverText({idx})")
        want_iv = f"{off_only['call_start']}–{off_only['call_end']}"
        acc = reg.filter(pl.col("row_id") == off_only["row_id"])["accession"][0]
        report.append({"hover_row": acc, "row_id": off_only["row_id"], "call": want_iv,
                       "hover_lists_it": want_iv in hover})
        if "1 call " not in hover or want_iv not in hover:
            failures.append(f"hover for {acc} row {off_only['row_id']}: count and merged calls disagree")

        # ---- 5. dark mode: canvases and curves use the dark kind colours ----
        # Switch the colour scheme on the loaded page: the page must redraw without a reload.
        page.evaluate(f"window.__pileup.setFocus({panels.index('hmmer3_phmmer|default')})")
        page.emulate_media(color_scheme="dark")
        page.wait_for_timeout(500)
        got = page.evaluate("""() => {
            const tok = getComputedStyle(document.documentElement).getPropertyValue('--k-letters').trim();
            const stroke = document.querySelector('#pu-fcurve path').getAttribute('stroke');
            const hist = document.querySelector('#pu-hist path[data-tool="hmmer3_phmmer"]').getAttribute('stroke');
            return {tok, stroke, hist};
        }""")
        report.append({"dark_mode_colours": got})
        if not (got["tok"] and got["stroke"] == got["tok"] and got["hist"] == got["tok"]):
            failures.append(f"dark mode: phmmer curve {got['stroke']} / histogram {got['hist']} != token {got['tok']}")
        page.emulate_media(color_scheme="light")
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
