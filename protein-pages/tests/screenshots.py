#!/usr/bin/env python3
"""Full-page screenshots of the index and some protein pages, 1280 px wide, light theme.

    python tests/screenshots.py SITE_DIR OUT_DIR PROTEIN_ID[:ROW] [PROTEIN_ID[:ROW] ...]

For each protein it saves the page, and the page with the residue panel of its first box open
(the first box of row ROW, such as kmerseek or hmmscan, when one is given).
The index is also saved at phone width (375 px). Needs Playwright and Google Chrome.
"""
import functools
import http.server
import socketserver
import sys
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def main():
    site, out, ids = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
    out.mkdir(parents=True, exist_ok=True)
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = "http://127.0.0.1:%d/" % httpd.server_address[1]
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        page = b.new_page(viewport={"width": 1280, "height": 900}, color_scheme="light")
        page.goto(base + "index.html")
        page.wait_for_timeout(800)
        page.screenshot(path=str(out / "0_index.png"), full_page=True)
        for n, arg in enumerate(ids, 1):
            pid, _, row = arg.partition(":")
            page.goto(base + "proteins/%s.html" % pid)
            page.wait_for_timeout(800)
            page.screenshot(path=str(out / ("%d_%s.png" % (n, pid))), full_page=True)
            box = page.query_selector('#tracks svg.trk g.row%s .hit rect' % ('[data-k="%s"]' % row if row else ""))
            if box:
                box.scroll_into_view_if_needed()
                box.click()
                page.wait_for_timeout(1500)
                page.screenshot(path=str(out / ("%d_%s_residues.png" % (n, pid))))
        phone = b.new_page(viewport={"width": 375, "height": 812}, color_scheme="light", is_mobile=True)
        phone.goto(base + "index.html")
        phone.wait_for_timeout(800)
        phone.screenshot(path=str(out / "0_index_phone.png"), full_page=True)
        if ids:
            phone.goto(base + "proteins/%s.html" % ids[0].partition(":")[0])
            phone.wait_for_timeout(800)
            phone.screenshot(path=str(out / ("1_%s_phone.png" % ids[0].partition(":")[0])), full_page=True)
        b.close()
    httpd.shutdown()


if __name__ == "__main__":
    main()
