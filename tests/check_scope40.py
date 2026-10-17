#!/usr/bin/env python3
"""Check the SCOPe40 section of the Methods page.

1. The worked example: label_hits() against labels written out by hand from the SCOPe ids
   (first number = class, first two = fold, first three = superfamily).
2. The numbers printed on the page equal data/scope40_testbed.json.
3. In headless Chrome, light and dark, 1280 px and 400 px: no page errors, no two SVG text
   elements of Figure 18 overlap, nothing leaves its SVG, the page does not scroll sideways,
   and every "Figure N" in the section names a figure on the site. Screenshots of the section
   go to screenshots/figures/scope40_*.png.
"""

import argparse
import functools
import json
import re
import socketserver
import sys
import threading
from pathlib import Path

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


def check_example(errors: list) -> None:
    for row in R.worked_example():
        got = "".join(CODE[x] for x in row["labels"])
        want = HAND[(row["level"], row["rule"])]
        if got != want:
            errors.append(f"labels {row['name']}: code {got}, by hand {want}")
    want_sens = {"Family (paper)": 2 / 4, "Superfamily (paper)": 2 / 5, "Fold (paper)": 1 / 3,
                 "Superfamily (nested)": 4 / 9}
    for row in R.worked_example():
        if abs(row["sens"] - want_sens[row["name"]]) > 1e-12:
            errors.append(f"sensitivity {row['name']}: {row['sens']} vs {want_sens[row['name']]}")


def check_numbers(site: Path, data: Path, errors: list) -> None:
    html = (site / "methods.html").read_text()
    sec = html[html.index('id="scope40"'):html.index('id="checks"')]
    tb = json.loads((data / "scope40_testbed.json").read_text())
    for key in ("scope208_astral40", "scope201_foldseek"):
        d = tb[key]
        for v in [d["n_domains"], *d["n_queries_with_a_true_positive"].values()]:
            if f"{v:,}" not in sec:
                errors.append(f"{key}: {v:,} not on the page")
    fs = next(x for x in tb["rocx"] if x["file"] == "foldseek.rocx")["levels"]
    for v in (fs["SFAM"]["mean_as_published"], fs["SFAM"]["mean_over_possible_true_positives"]):
        if f"{v:.3f}" not in sec:
            errors.append(f"Foldseek superfamily mean {v:.3f} not on the page")
    if "3.3%" not in sec or "retired" not in sec:
        errors.append("the retired-claim sentence is missing")


def check_render(site: Path, out: Path, errors: list) -> None:
    out.mkdir(parents=True, exist_ok=True)
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    all_figs = set()
    for f in site.rglob("*.html"):
        all_figs |= set(re.findall(r"<b>Figure (\d+)\.</b>", f.read_text()))
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        for scheme in ("light", "dark"):
            for width in (1280, 400):
                page = b.new_page(viewport={"width": width, "height": 900}, color_scheme=scheme,
                                  device_scale_factor=1.5)
                errs = []
                page.on("pageerror", lambda e: errs.append(str(e)))
                page.goto(base + "methods.html#scope40")
                page.wait_for_load_state("networkidle")
                if errs:
                    errors.append(f"{scheme} {width}: page errors {errs}")
                res = page.evaluate("""() => {
                  const fig = document.getElementById('fig-scope40-labels');
                  const svg = fig.querySelector('.ex-box svg');
                  const sb = svg.getBoundingClientRect();
                  const ts = [...svg.querySelectorAll('text')].map(t => t.getBoundingClientRect());
                  let overlaps = 0, outside = 0;
                  for (let i = 0; i < ts.length; i++) {
                    const a = ts[i];
                    if (a.left < sb.left - 0.5 || a.right > sb.right + 0.5 || a.top < sb.top - 0.5 || a.bottom > sb.bottom + 0.5) outside++;
                    for (let j = i + 1; j < ts.length; j++) {
                      const c = ts[j];
                      if (a.left < c.right - 0.5 && c.left < a.right - 0.5 && a.top < c.bottom - 0.5 && c.top < a.bottom - 0.5) overlaps++;
                    }
                  }
                  const doc = document.documentElement;
                  return {overlaps, outside, n: ts.length, sideways: doc.scrollWidth > doc.clientWidth + 1,
                          labels: [...svg.querySelectorAll('[data-label]')].length};
                }""")
                if res["overlaps"] or res["outside"]:
                    errors.append(f"{scheme} {width}: {res['overlaps']} overlapping, {res['outside']} outside SVG")
                if res["sideways"]:
                    errors.append(f"{scheme} {width}: page scrolls sideways")
                if res["labels"] != 8 * 4:
                    errors.append(f"{scheme} {width}: {res['labels']} hit boxes, expected 32")
                if width == 1280:
                    page.locator("#scope40").scroll_into_view_if_needed()
                    page.locator("#fig-scope40-labels").screenshot(path=str(out / f"scope40_fig18_{scheme}.png"))
                page.close()
        b.close()
    httpd.shutdown()
    sec = (site / "methods.html").read_text()
    sec = sec[sec.index('id="scope40"'):sec.index('id="checks"')]
    for n in set(re.findall(r"Figure (\d+)", sec)):
        if n not in all_figs:
            errors.append(f"Figure {n} is named but no caption has it")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, default=ROOT / "screenshots" / "figures")
    a = ap.parse_args()
    errors: list[str] = []
    check_example(errors)
    check_numbers(a.site, a.data, errors)
    check_render(a.site, a.out, errors)
    for e in errors:
        print("FAIL", e)
    print("scope40 checks:", "FAIL" if errors else "ok")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
