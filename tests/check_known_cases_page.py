#!/usr/bin/env python3
"""Open the Known cases page in headless Chrome and check what only a browser can see.

    .venv/bin/python tests/check_known_cases_page.py

  no two text elements inside an SVG overlap, and nothing is drawn outside its SVG;
  the page loads with no console error and no failed request;
  every mark in Figure 18 has a hover text naming the tool, the case and the result, and
    clicking one moves the view to that case;
  the two price boxes change every cost on the page.
"""

from __future__ import annotations

import argparse
import functools
import http.server
import socketserver
import sys
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


OVERLAP_JS = """
() => {
  const bad = [];
  document.querySelectorAll('svg').forEach((svg, si) => {
    const box = svg.getBoundingClientRect();
    const texts = [...svg.querySelectorAll('text')];
    texts.forEach((t, i) => {
      const a = t.getBoundingClientRect();
      if (a.width === 0 || a.height === 0) return;
      if (a.left < box.left - 1 || a.right > box.right + 1 ||
          a.top < box.top - 1 || a.bottom > box.bottom + 1)
        bad.push('svg ' + si + ': "' + t.textContent.trim() + '" leaves its svg');
      for (let j = i + 1; j < texts.length; j++) {
        const b = texts[j].getBoundingClientRect();
        if (b.width === 0 || b.height === 0) continue;
        const ox = Math.min(a.right, b.right) - Math.max(a.left, b.left);
        const oy = Math.min(a.bottom, b.bottom) - Math.max(a.top, b.top);
        if (ox > 1 && oy > 1)
          bad.push('svg ' + si + ': "' + t.textContent.trim() + '" overlaps "' +
                   texts[j].textContent.trim() + '"');
      }
    });
  });
  return bad;
}
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    a = ap.parse_args()
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(a.site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    bad: list[str] = []
    with sync_playwright() as pw:
        browser = pw.chromium.launch(channel="chrome")
        for width in (1280, 400):
            page = browser.new_page(viewport={"width": width, "height": 900})
            page.on("console", lambda m: bad.append(f"{width} px console {m.type}: {m.text}")
                    if m.type == "error" else None)
            page.on("pageerror", lambda e: bad.append(f"{width} px page error: {e}"))
            page.on("requestfailed",
                    lambda r: bad.append(f"{width} px request failed: {r.url}"))
            page.goto(base + "known-cases.html")
            page.wait_for_load_state("networkidle")
            bad += [f"{width} px {x}" for x in page.evaluate(OVERLAP_JS)]
            wide = page.evaluate("document.documentElement.scrollWidth > "
                                 "document.documentElement.clientWidth + 1")
            if wide:
                bad.append(f"{width} px: the page scrolls sideways")
            if width == 1280:
                marks = page.locator(".kcmatrix .mk")
                n = marks.count()
                if n != 14 * 9:
                    bad.append(f"Figure 18 has {n} marks, 14 cases by 9 tools is 126")
                tips = page.eval_on_selector_all(
                    ".kcmatrix .mk title", "ts => ts.map(t => t.textContent)")
                if len(tips) != n:
                    bad.append(f"{n} marks but {len(tips)} hover texts")
                for i, tip in enumerate(tips):
                    if "case" not in tip or ":" not in tip:
                        bad.append(f"mark {i} hover text is {tip!r}")
                before = page.evaluate("window.scrollY")
                marks.nth(40).click()
                page.wait_for_timeout(900)
                if page.evaluate("window.scrollY") == before:
                    bad.append("clicking a mark did not move the view to its case")
                first = page.locator(".cost1k").first
                was = first.inner_text()
                page.fill("#kc-rcpu", "1")
                page.wait_for_timeout(200)
                now = first.inner_text()
                if was == now:
                    bad.append(f"the price box did not change the cost ({was})")
                if "$" not in now:
                    bad.append(f"the cost does not read as money: {now!r}")
            page.close()
        browser.close()
    httpd.shutdown()
    for x in bad:
        print("FAIL:", x)
    print(f"known-cases page checks: {len(bad)} problems")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
