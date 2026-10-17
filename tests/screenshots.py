#!/usr/bin/env python3
"""Full-page screenshots of every page at 1280 and 400 px, light and dark, into screenshots/."""

import argparse
import functools
import http.server
import socketserver
import sys
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
# In the order of the site's navigation.
PAGES = {
    "1_leaderboard": "index.html",
    "2_tool_cards_index": "tools/index.html",
    "2b_tool_card_kmerseek": "tools/kmerseek.html",
    "3_compare_two_columns": "cost.html",
    "4_how_the_scores_work": "explainers.html",
    "5_region_pileup": "pileup.html",
    "6_methods": "methods.html",
    "7_claims": "claims.html",
    "8_coverage": "coverage.html",
    "9_known_cases": "known-cases.html",
}
EXTRA = {
    "3b_compare_cost_usd_auprc": "cost.html#x=cost&y=auprc",
    "3c_compare_fmax_auprc": "cost.html#x=fmax&y=auprc",
    "3d_compare_ram_iou": "cost.html#x=ram&y=iou",
}
# Every page at desktop and phone width, light and dark: <name>.png is 1280 px light,
# <name>_400.png 400 px light, and _dark is added for the dark theme.
WIDTHS = (1280, 400)
SCHEMES = ("light", "dark")


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def shoot(page, base, path, file):
    page.goto(base + path)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(800)
    if path.startswith("pileup.html"):
        page.wait_for_function("window.__pileup && window.__pileup.ready", timeout=60000)
    page.screenshot(path=str(file), full_page=True)
    print(file)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    ap.add_argument("--out", type=Path, default=ROOT / "screenshots")
    a = ap.parse_args()
    out = a.out
    out.mkdir(parents=True, exist_ok=True)
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(a.site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        for w in WIDTHS:
            for scheme in SCHEMES:
                page = b.new_page(viewport={"width": w, "height": 900 if w > 600 else 860}, color_scheme=scheme)
                suffix = ("" if w > 600 else f"_{w}") + ("" if scheme == "light" else "_dark")
                for name, path in PAGES.items():
                    shoot(page, base, path, out / f"{name}{suffix}.png")
                page.close()
        page = b.new_page(viewport={"width": 1280, "height": 900}, color_scheme="light")
        for name, path in EXTRA.items():
            shoot(page, base, path, out / f"{name}.png")
        # A definition open, as it looks after a tap or Tab.
        page.goto(base + "index.html")
        page.wait_for_load_state("networkidle")
        page.focus("button.term[aria-controls='tip-fmax-th']")
        page.wait_for_timeout(200)
        page.screenshot(path=str(out / "9_definition_fmax.png"), clip={"x": 0, "y": 0, "width": 1280, "height": 700})
        print(out / "9_definition_fmax.png")
        b.close()
    httpd.shutdown()


if __name__ == "__main__":
    main()
