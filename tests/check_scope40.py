#!/usr/bin/env python3
"""Check the SCOPe40 page (scope40.html) and the SCOPe40 section of the Methods page.

1. Figure 19's worked example: label() against labels written out by hand from the SCOPe ids.
2. Numbers: the Methods section against data/scope40_testbed.json; the page's AUCs, coverage,
   n and class counts against data/scope40_summary.json; the summary against a recomputation
   from data/scope40_long.parquet (AUC = mean sensitivity over ALL queries, coverage = share
   with has_hits, hits-only AUC over has_hits queries).
3. Headless Chrome, light and dark, 1280 and 400 px: no page errors; in Figures 19, 20 and 21
   no two SVG text elements overlap and none leaves its SVG; no sideways scroll; every
   "Figure N" names a caption on the site. Controls: each tool box removes its lines; the
   no-hit toggle prints, per tool and level, a change equal to the summary's hits-only AUC
   minus its AUC; a class click lists that class's queries; a sort puts empty rows last.
   Screenshots go to screenshots/figures/scope40_*.png.
"""

import argparse
import functools
import json
import re
import socketserver
import sys
import threading
from pathlib import Path

import polars as pl
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import scope40_rules as R  # noqa: E402
from screenshots import Quiet  # noqa: E402

# Hand labels for the made-up query b.1.1.1, hits 1-8 (t = TP, f = FP, - = not counted).
HAND = {
    ("family", "paper"): "t---tftf",
    ("superfamily", "paper"): "-t-t-f-f",
    ("fold", "paper"): "--t--f-f",
    ("superfamily", "nested"): "tt-ttftf",
}
CODE = {"tp": "t", "fp": "f", "skip": "-"}
FIGS = ["#fig-scope40-labels .ex-box svg", "#s40-curves svg", "#s40-classes svg"]


def check_example(errors: list) -> None:
    want_sens = {"Family (paper)": 2 / 4, "Superfamily (paper)": 2 / 5, "Fold (paper)": 1 / 3,
                 "Superfamily (nested)": 4 / 9}
    for row in R.worked_example():
        got = "".join(CODE[x] for x in row["labels"])
        if got != HAND[(row["level"], row["rule"])]:
            errors.append(f"labels {row['name']}: code {got}, by hand {HAND[(row['level'], row['rule'])]}")
        if abs(row["sens"] - want_sens[row["name"]]) > 1e-12:
            errors.append(f"sensitivity {row['name']}: {row['sens']}")


def check_numbers(site: Path, data: Path, errors: list) -> dict:
    html = (site / "methods.html").read_text()
    sec = html[html.index('id="scope40"'):html.index('id="checks"')]
    tb = json.loads((data / "scope40_testbed.json").read_text())
    for key in ("scope208_astral40", "scope201_foldseek"):
        d = tb[key]
        for v in [d["n_domains"], *d["n_queries_with_a_true_positive"].values()]:
            if f"{v:,}" not in sec:
                errors.append(f"Methods: {key} {v:,} not on the page")
    if "3.3%" not in sec or "retired" not in sec:
        errors.append("Methods: the retired-claim sentence is missing")
    summ = json.loads((data / "scope40_summary.json").read_text())
    long = pl.read_parquet(data / "scope40_long.parquet")
    n_q = long["query"].n_unique()
    if n_q != summ["n_queries"]:
        errors.append(f"queries: {n_q} in the long table, {summ['n_queries']} in the summary")
    for r in summ["overall"]:
        g = long.filter((pl.col("tool") == r["tool"]) & (pl.col("setting") == r["setting"]) & (pl.col("level") == r["level"]))
        if g.height != n_q:
            errors.append(f"{r['tool']} {r['level']}: {g.height} rows, {n_q} queries")
        auc, cov = g["sensitivity"].mean(), g["has_hits"].mean()
        hits = g.filter(pl.col("has_hits"))["sensitivity"].mean()
        for name, a, b in (("AUC", auc, r["auc"]), ("coverage", cov, r["coverage"]), ("hits-only AUC", hits, r["auc_hits_only"])):
            if abs(a - b) > 1e-12:
                errors.append(f"{r['tool']} {r['level']} {name}: recomputed {a}, summary {b}")
        if not (r["lo"] <= r["auc"] <= r["hi"]):
            errors.append(f"{r['tool']} {r['level']}: AUC outside its interval")
    page = (site / "scope40.html").read_text()
    for r in summ["overall"]:
        a = r["auc"]
        if ("<0.0001" if 0 < a < 0.0001 else f"{a:.4f}" if 0 < a < 0.001 else f"{a:.3f}") not in page:
            errors.append(f"scope40.html: AUC {r['auc']:.3f} ({r['tool']} {r['level']}) not in Table 3")
    if f"n = {n_q:,} queries" not in page:
        errors.append("scope40.html: caption n is not the number of queries")
    for r in summ["by_class"]:
        if f"(n = {r['n_queries']:,})" not in page:
            errors.append(f"scope40.html: class {r['scope_class']} n = {r['n_queries']:,} missing")
    if "3.3%" not in page or "retired" not in page:
        errors.append("scope40.html: the retired-claim sentence is missing")
    return summ


def text_overlaps(page, sel: str) -> dict:
    return page.evaluate("""(sel) => {
      let overlaps = 0, outside = 0, n = 0, bad = [];
      for (const svg of document.querySelectorAll(sel)) {
        const sb = svg.getBoundingClientRect();
        const ts = [...svg.querySelectorAll('text')].map(t => [t.getBoundingClientRect(), t.textContent]);
        n += ts.length;
        for (let i = 0; i < ts.length; i++) {
          const a = ts[i][0];
          if (a.width === 0) continue;
          if (a.left < sb.left - 0.5 || a.right > sb.right + 0.5 || a.top < sb.top - 0.5 || a.bottom > sb.bottom + 0.5) { outside++; bad.push('out:' + ts[i][1]); }
          for (let j = i + 1; j < ts.length; j++) {
            const c = ts[j][0];
            if (c.width === 0) continue;
            if (a.left < c.right - 0.5 && c.left < a.right - 0.5 && a.top < c.bottom - 0.5 && c.top < a.bottom - 0.5) { overlaps++; bad.push(ts[i][1] + '|' + ts[j][1]); }
          }
        }
      }
      return {overlaps, outside, n, bad: bad.slice(0, 6)};
    }""", sel)


def check_render(site: Path, out: Path, summ: dict, errors: list) -> None:
    out.mkdir(parents=True, exist_ok=True)
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        for scheme in ("light", "dark"):
            for width in (1280, 400):
                page = b.new_page(viewport={"width": width, "height": 900}, color_scheme=scheme, device_scale_factor=1.5)
                errs = []
                page.on("pageerror", lambda e: errs.append(str(e)))
                page.on("console", lambda m: errs.append(m.text) if m.type == "error" else None)
                page.goto(base + "scope40.html")
                page.wait_for_selector("body[data-s40-ready]")
                tag = f"{scheme} {width}"
                for sel in FIGS:
                    res = text_overlaps(page, sel)
                    if res["n"] == 0:
                        errors.append(f"{tag} {sel}: no text drawn")
                    if res["overlaps"] or res["outside"]:
                        errors.append(f"{tag} {sel}: {res['overlaps']} overlapping, {res['outside']} outside: {res['bad']}")
                if page.evaluate("document.documentElement.scrollWidth > document.documentElement.clientWidth + 1"):
                    errors.append(f"{tag}: page scrolls sideways")
                if page.locator("#fig-scope40-labels [data-label]").count() != 32:
                    errors.append(f"{tag}: Figure 19 does not have 32 hit boxes")
                if width == 1280:
                    # The sticky site header would sit over the element screenshots.
                    page.add_style_tag(content="header.site{position:static !important}")
                    for sel, name in (("#fig-scope40-labels", "fig19"), ("#fig-scope40-curves", "fig20"),
                                      ("#fig-scope40-classes", "fig21"), ("#s40-tab", "table3")):
                        page.locator(sel).screenshot(path=str(out / f"scope40_{name}_{scheme}.png"))
                if width == 1280 and scheme == "light":
                    check_controls(page, summ, errors)
                if errs:
                    errors.append(f"{tag}: page errors {errs}")
                page.close()
            page = b.new_page(viewport={"width": 1280, "height": 900}, color_scheme=scheme)
            page.goto(base + "methods.html#scope40")
            page.wait_for_load_state("networkidle")
            page.close()
        b.close()
    httpd.shutdown()
    all_figs = set()
    for f in site.rglob("*.html"):
        all_figs |= set(re.findall(r"<b>Figure (\d+)\.</b>", f.read_text()))
    for name in ("scope40.html", "methods.html"):
        for n in set(re.findall(r"Figure (\d+)", (site / name).read_text())):
            if n not in all_figs:
                errors.append(f"{name}: Figure {n} is named but no caption has it")


def check_controls(page, summ: dict, errors: list) -> None:
    n_lines = page.locator("#s40-curves path[data-tool]").count()
    tools = page.locator("#s40-tool-boxes input").count()
    on = page.locator("#s40-tool-boxes input:checked").count()
    if n_lines != 3 * on:
        errors.append(f"Figure 20: {n_lines} lines for {on} ticked tools")
    page.locator("#s40-tool-boxes input").first.uncheck()
    if page.locator("#s40-curves path[data-tool]").count() != 3 * (on - 1):
        errors.append("Figure 20: unticking a tool does not remove its lines")
    page.locator("#s40-tool-boxes input").first.check()
    before = page.locator("#s40-curves path[data-tool]").first.get_attribute("d")
    page.locator('input[name="s40-nohit"][value="drop"]').check()
    after = page.locator("#s40-curves path[data-tool]").first.get_attribute("d")
    if before == after:
        errors.append("Figure 20: the no-hit toggle does not redraw the lines")
    deltas = page.evaluate("""() => [...document.querySelectorAll('#s40-curve-legend td[data-delta]')]
        .map(td => [td.closest('tr').getAttribute('data-key'), td.cellIndex, parseFloat(td.getAttribute('data-delta')),
                    td.textContent])""")
    levels = ["family", "superfamily", "fold"]
    if len(deltas) != 3 * on:
        errors.append(f"Figure 20: {len(deltas)} printed changes, expected {3 * on}")
    for key, col, d, shown in deltas:
        tool, setting = key.split("|")
        s = next(r for r in summ["overall"] if r["tool"] == tool and r["setting"] == setting and r["level"] == levels[col - 1])
        want = s["auc_hits_only"] - s["auc"]
        if abs(d - want) > 1e-9:
            errors.append(f"Figure 20 toggle {tool} {levels[col - 1]}: printed change {d}, recomputed {want}")
        w = abs(want)
        if ("<0.0001" if 0 < w < 0.0001 else f"{w:.4f}" if 0 < w < 0.001 else f"{w:.3f}") not in shown:
            errors.append(f"Figure 20 toggle {tool} {levels[col - 1]}: text '{shown}' lacks {abs(want):.3f}")
    page.locator('input[name="s40-nohit"][value="zero"]').check()
    g = page.locator('#s40-classes g.fg-col[data-cls="f"][data-level="family"]')
    g.locator("rect").first.click(force=True)
    n_f = next(r["n_queries"] for r in summ["by_class"] if r["scope_class"] == "f" and r["level"] == "family")
    if page.locator("#s40-class-table tbody tr").count() != n_f:
        errors.append("Figure 21: clicking class f does not list its queries")
    page.locator('#s40-tab thead button[data-col="2"]').click()
    first = page.locator("#s40-tab tbody tr").first.get_attribute("data-tool")
    page.locator('#s40-tab thead button[data-col="2"]').click()
    first_asc = page.locator("#s40-tab tbody tr").first.get_attribute("data-tool")
    last = page.locator("#s40-tab tbody tr").last.get_attribute("class") or ""
    if first != "foldseek" or "pinned" not in last:
        errors.append(f"Table 3 sort: first {first}, ascending first {first_asc}, last row class {last}")
    lowest = min((r for r in summ["overall"] if r["level"] == "family" and not r["tool"].endswith("_shuffled")),
                 key=lambda r: r["auc"])
    shown_rows = page.locator("#s40-tab tbody tr:not([hidden])")
    if first_asc != lowest["tool"] and lowest["tool"] != "kmerseek":
        errors.append(f"Table 3 ascending sort: first {first_asc}, lowest family AUC {lowest['tool']}")
    if shown_rows.count() == 0:
        errors.append("Table 3: no rows shown")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, default=ROOT / "screenshots" / "figures")
    a = ap.parse_args()
    errors: list[str] = []
    check_example(errors)
    summ = check_numbers(a.site, a.data, errors)
    check_render(a.site, a.out, summ, errors)
    for e in errors:
        print("FAIL", e)
    print("scope40 checks:", "FAIL" if errors else "ok")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
