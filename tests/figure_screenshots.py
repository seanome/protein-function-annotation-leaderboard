#!/usr/bin/env python3
"""Screenshot each of Figures 11 to 17 on its own, light and dark at 1280 px, into
screenshots/figures/, and check at 400 px that no page with one of them scrolls sideways.

A figure is found by its <figure class="fgfig" id="..."> element. The hover box is
checked once: hovering the first mark of each figure must show its data-tip text.
"""

import argparse
import functools
import socketserver
import sys
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright

sys.path.insert(0, str(Path(__file__).resolve().parent))
from screenshots import Quiet  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PAGES = ["index.html", "tools/kmerseek.html", "tools/hmmer3_phmmer.html", "methods.html", "proteins.html"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    ap.add_argument("--out", type=Path, default=ROOT / "screenshots" / "figures")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(a.site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    errors = []
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        for scheme in ("light", "dark"):
            page = b.new_page(viewport={"width": 1280, "height": 900}, color_scheme=scheme, device_scale_factor=1.5)
            for path in PAGES:
                if not (a.site / path).exists():
                    continue
                page.goto(base + path)
                page.wait_for_load_state("networkidle")
                for fig in page.query_selector_all("figure.fgfig"):
                    fid = fig.get_attribute("id")
                    stem = path.replace("/", "_").replace(".html", "")
                    f = a.out / f"{stem}_{fid}{'' if scheme == 'light' else '_dark'}.png"
                    fig.scroll_into_view_if_needed()
                    fig.screenshot(path=str(f))
                    print(f)
                    if scheme == "light":
                        mark = fig.query_selector("[data-tip]")
                        if mark is None:
                            errors.append(f"{path} {fid}: no hover marks")
                            continue
                        mark.hover(force=True)
                        page.wait_for_timeout(100)
                        shown = page.evaluate("(() => { const t = document.querySelector('.fg-tip'); "
                                              "return t && t.style.display === 'block' ? t.textContent : ''; })()")
                        # the mark under the pointer may be another one drawn on top
                        bb = mark.bounding_box()
                        top = page.evaluate(
                            "([x, y]) => { const e = document.elementFromPoint(x, y); "
                            "const m = e && e.closest('[data-tip]'); return m ? m.getAttribute('data-tip') : ''; }",
                            [bb["x"] + bb["width"] / 2, bb["y"] + bb["height"] / 2])
                        if not shown or shown != top:
                            errors.append(f"{path} {fid}: hover box shows {shown!r}")
            page.close()
        page = b.new_page(viewport={"width": 400, "height": 860})
        for path in PAGES:
            if not (a.site / path).exists():
                continue
            page.goto(base + path)
            page.wait_for_load_state("networkidle")
            sw = page.evaluate("document.documentElement.scrollWidth")
            if sw > 400:
                errors.append(f"{path}: page is {sw} px wide at a 400 px window")
        b.close()
    httpd.shutdown()
    for e in errors:
        print("FAIL", e)
    if errors:
        sys.exit(1)
    print("figure screenshots taken; hover boxes and 400 px width passed")


if __name__ == "__main__":
    main()
