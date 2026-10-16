#!/usr/bin/env python3
"""How far down the page the first piece of data starts, at phone width (400 px).

Prints, for the leaderboard, the region pile-up and the tool cards index, the distance in
CSS pixels from the top of the page to the top of the first data element (the first table
row, the first pile-up panel, the first tool card) and the header's height. A phone screen
is about 800 px tall, so a number under that means the data is on the first screen.
"""

import argparse
import functools
import http.server
import json
import socketserver
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
FIRST = {
    "leaderboard.html": "#lb tbody tr:not([hidden])",
    "pileup.html": "#pu-grid .pu-panel",
    "tools/index.html": "#tc-grid .tc",
}


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    ap.add_argument("--width", type=int, default=400)
    a = ap.parse_args()
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(a.site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    out = {}
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        page = b.new_page(viewport={"width": a.width, "height": 860})
        for pg, sel in FIRST.items():
            page.goto(base + pg)
            page.wait_for_load_state("networkidle")
            if pg == "pileup.html":
                page.wait_for_function("window.__pileup && window.__pileup.ready", timeout=60000)
            out[pg] = page.evaluate("""(sel) => {
                const e = document.querySelector(sel), h = document.querySelector('header.site');
                return {first_data_top_px: Math.round(e.getBoundingClientRect().top + scrollY),
                        header_px: Math.round(h.getBoundingClientRect().height),
                        page_px: document.documentElement.scrollHeight};
            }""", sel)
        b.close()
    httpd.shutdown()
    print(json.dumps({"width": a.width, **out}, indent=1))


if __name__ == "__main__":
    main()
