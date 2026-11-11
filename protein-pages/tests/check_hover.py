#!/usr/bin/env python3
"""Check that pointing at a hit box does not move the page, and that the labels of InterProScan's
predictions (disorder, coiled coil, membrane, signal) in the domain lane cover nothing.

    python tests/check_hover.py SITE_DIR PAGE [PAGE ...]

PAGE is a path inside the built site, such as proteins/FUN036660_MYO10.2.html or
protein.html#FUN005955.2. At each window width in WIDTHS the page is opened, every program
is unfolded, and the pointer is put on each hit box in turn (at most MAX_HITS, spread over
the page). Before and after each one, the top and height of the fixed header, of every
grid row and of every program row are read. A box fails when any of them changed by more
than half a pixel. The page had this fault until 2026-11-10: a line of text in the fixed
header was hidden while empty, so the header grew when a hit was pointed at.

In the domain lane, the label of each dashed prediction box must not lie over a dashed box other
than its own (inside its own is allowed when all of it fits), must not touch another label, and
must stay inside the drawing. Until 2026-11-10 a label was centred on its box whatever its width.

Prints one line per page and width, and exits 1 if any box moved anything or any label fails.
Needs Playwright, and Google Chrome or Playwright's own Chromium (python -m playwright install chromium).
pytest runs it too: tests/test_build_full.py calls check() on pages of the sites it builds.
"""
import functools
import http.server
import socketserver
import sys
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright

WIDTHS = (1050, 1280, 1400)
MAX_HITS = 150

PROBE = """async (maxHits) => {
  const wait = (ms) => new Promise((r) => setTimeout(r, ms));
  const all = document.getElementById("expand-all");
  if (all) { all.click(); await wait(400); }
  window.scrollTo(0, 300); await wait(200);
  const hits = [...document.querySelectorAll(".hit[data-i]")];
  const step = Math.max(1, Math.ceil(hits.length / maxHits));
  const boxes = () => [...document.querySelectorAll(".top, .top > *, .g, .row")].map((e) => {
    const r = e.getBoundingClientRect(); return [e, r.top, r.height]; });
  let tried = 0, moved = 0, worst = 0, example = "";
  for (let n = 0; n < hits.length; n += step) {
    const h = hits[n];
    if (!h.isConnected) continue;
    tried++;
    const before = boxes();
    h.dispatchEvent(new MouseEvent("mouseover", { bubbles: true, clientX: 300, clientY: 300 }));
    await wait(25);
    const after = new Map(boxes().map((x) => [x[0], x]));
    let d = 0, what = "";
    for (const [e, top, height] of before) {
      const a = after.get(e);
      if (!a) continue;  // the ruler row is drawn again for a kmerseek hit; its size is read through .top and .g
      const c = Math.max(Math.abs(a[1] - top), Math.abs(a[2] - height));
      if (c > d) { d = c; what = String(e.className.baseVal ?? e.className).slice(0, 30); }
    }
    if (d > 0.5) { moved++; if (d > worst) { worst = d; example = what; } }
    h.dispatchEvent(new MouseEvent("mouseout", { bubbles: true }));
    await wait(5);
  }
  // labels of the dashed prediction boxes in the domain lane
  const svg = document.querySelector(".row.lane svg"), bad = [];
  let labels = 0;
  if (svg) {
    const S = svg.getBoundingClientRect(), items = [...svg.querySelectorAll("g.dom")].filter((g) => g.querySelector("rect[stroke-dasharray]") && g.querySelector("text"))
      .map((g) => ({ box: g.querySelector("rect[stroke-dasharray]").getBoundingClientRect(), text: g.querySelector("text").getBoundingClientRect(), l: g.querySelector("text").textContent }));
    labels = items.length;
    const hit = (a, b) => a.left < b.right - 0.5 && b.left < a.right - 0.5 && a.top < b.bottom - 0.5 && b.top < a.bottom - 0.5;
    items.forEach((it, i) => {
      const inside = it.text.left >= it.box.left - 0.5 && it.text.right <= it.box.right + 0.5;
      if (hit(it.text, it.box) && !inside) bad.push(it.l + ": over its own box");
      if (it.text.left < S.left - 0.5 || it.text.right > S.right + 0.5) bad.push(it.l + ": outside the drawing");
      items.forEach((o, j) => { if (j !== i && hit(it.text, o.box)) bad.push(it.l + ": over the box of " + o.l); if (j > i && hit(it.text, o.text)) bad.push(it.l + ": over the label " + o.l); });
    });
  }
  return { hits: hits.length, tried, moved, worst, example, labels, bad };
}"""


def launch(p):
    """Google Chrome when it is installed, else the Chromium that `playwright install chromium` fetches."""
    try:
        return p.chromium.launch(channel="chrome", headless=True)
    except Exception:
        return p.chromium.launch(headless=True)


def failures(results):
    """One line for each page and width where a hit box moved the page or a label is misplaced."""
    return ["%s at %d px: %d of %d hit boxes moved the page (most %.0f px, %s); labels: %s" % (
        pg, w, r["moved"], r["tried"], r["worst"], r["example"], "; ".join(r["bad"]) or "fine")
        for pg, w, r in results if r["tried"] == 0 or r["moved"] or r["bad"]]


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def check(site, pages, widths=WIDTHS, max_hits=MAX_HITS):
    """[(page, width, result of PROBE), ...] for every page at every width."""
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = "http://127.0.0.1:%d/" % httpd.server_address[1]
    out = []
    try:
        with sync_playwright() as p:
            b = launch(p)
            for w in widths:
                page = b.new_page(viewport={"width": w, "height": 900})
                for pg in pages:
                    page.goto("about:blank")  # so a change of #protein loads the page again
                    page.goto(base + pg)
                    page.wait_for_selector(".hit[data-i]", state="attached", timeout=20000)
                    out.append((pg, w, page.evaluate(PROBE, max_hits)))
                page.close()
            b.close()
    finally:
        httpd.shutdown()
    return out


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    bad = 0
    for pg, w, r in check(Path(sys.argv[1]), sys.argv[2:]):
        ok = r["tried"] > 0 and r["moved"] == 0 and not r["bad"]
        bad += not ok
        print("%s at %d px: pointed at %d of %d hit boxes, %d moved the page%s; %d prediction labels, %d misplaced%s -> %s" % (
            pg, w, r["tried"], r["hits"], r["moved"],
            " (most: %.0f px, %s)" % (r["worst"], r["example"]) if r["moved"] else "",
            r["labels"], len(r["bad"]), " (%s)" % "; ".join(r["bad"][:3]) if r["bad"] else "", "OK" if ok else "FAIL"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
