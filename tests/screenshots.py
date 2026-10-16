#!/usr/bin/env python3
"""Full-page screenshots of each page at 1280 px, light theme, into screenshots/."""

import functools
import http.server
import socketserver
import sys
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
PAGES = {
    "1_leaderboard": "index.html",
    "2_compare_two_columns": "cost.html",
    "6_region_pileup": "pileup.html",
    "3_tool_card_kmerseek": "tools/kmerseek.html",
    "3b_tool_cards_index": "tools/index.html",
    "7_claims": "claims.html",
    "4_coverage": "coverage.html",
    "5_methods": "methods.html",
    "8_how_the_scores_work": "explainers.html",
}
EXTRA = {
    "2b_compare_cost_usd_auprc": "cost.html#x=cost&y=auprc",
    "2c_compare_fmax_auprc": "cost.html#x=fmax&y=auprc",
    "2d_compare_ram_iou": "cost.html#x=ram&y=iou",
}
# Also at phone width (400 px), as <name>_400.png.
PHONE = {"1_leaderboard", "2_compare_two_columns", "8_how_the_scores_work"}


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def main() -> None:
    out = ROOT / "screenshots"
    out.mkdir(exist_ok=True)
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(ROOT / "docs")))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        page = b.new_page(viewport={"width": 1280, "height": 900}, color_scheme="light")
        for name, path in PAGES.items():
            page.goto(base + path)
            page.wait_for_load_state("networkidle")
            page.wait_for_timeout(800)
            if path == "pileup.html":
                page.wait_for_function("window.__pileup && window.__pileup.ready", timeout=60000)
            page.screenshot(path=str(out / f"{name}.png"), full_page=True)
            print(out / f"{name}.png")
        # Figure 1 of cost.html with other axis pairs, and at phone width.
        for name, path in EXTRA.items():
            page.goto(base + path)
            page.wait_for_load_state("networkidle")
            page.wait_for_timeout(800)
            page.screenshot(path=str(out / f"{name}.png"), full_page=True)
            print(out / f"{name}.png")
        phone = b.new_page(viewport={"width": 400, "height": 860}, color_scheme="light")
        for name, path in PAGES.items():
            if name not in PHONE:
                continue
            phone.goto(base + path)
            phone.wait_for_load_state("networkidle")
            phone.wait_for_timeout(800)
            phone.screenshot(path=str(out / f"{name}_400.png"), full_page=True)
            print(out / f"{name}_400.png")
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
    sys.exit(main())
