#!/usr/bin/env python3
"""Drive the /contribute/truth/ form in headless Chrome and check what it would send.

1. A wrong ORCID iD (bad check digit) and an empty region must stop the form with messages.
2. Good answers must open github.com/<repo>/issues/new with a title, a body and the label,
   and the body must parse back, field for field, with scripts/truth_contributed.py.
3. Screenshots of both contribute pages at 1280 and 400 px, light and dark, go to --shots.

--url-out saves the issue URL the form built, so the same text can be filed as a test issue.
"""

import argparse
import functools
import http.server
import socketserver
import sys
import threading
import urllib.parse
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import truth_contributed as tc  # noqa: E402

GOOD = {"#tf-acc": "P01375", "#tf-start": "89", "#tf-end": "233",
        "#tf-eco": "ECO:0000269", "#tf-cite": "TEST ROW for the Prompt J check of the form; not a real result.",
        "#tf-orcid": "0000-0002-1825-0097"}


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    ap.add_argument("--shots", type=Path)
    ap.add_argument("--url-out", type=Path)
    a = ap.parse_args()
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(a.site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    bad = 0
    with sync_playwright() as p:
        browser = p.chromium.launch(channel="chrome")
        if a.shots:
            a.shots.mkdir(parents=True, exist_ok=True)
            for scheme in ("light", "dark"):
                for w in (1280, 400):
                    ctx = browser.new_context(viewport={"width": w, "height": 900}, color_scheme=scheme)
                    pg = ctx.new_page()
                    for name, path in (("contribute", "contribute/"), ("truth_form", "contribute/truth/")):
                        pg.goto(base + path)
                        pg.wait_for_load_state("networkidle")
                        sw = pg.evaluate("document.documentElement.scrollWidth")
                        if sw > w:
                            print(f"FAIL {path} at {w}px scrolls sideways: page is {sw}px wide")
                            bad += 1
                        sfx = ("" if w == 1280 else f"_{w}") + ("" if scheme == "light" else "_dark")
                        pg.screenshot(path=str(a.shots / f"{name}{sfx}.png"), full_page=True)
                    ctx.close()
        ctx = browser.new_context(viewport={"width": 1280, "height": 900})
        pg = ctx.new_page()
        pg.add_init_script("window.__opened = []; window.open = (u) => { window.__opened.push(u); return null; };")
        pg.goto(base + "contribute/truth/")
        pg.wait_for_load_state("networkidle")
        # 1. Bad answers.
        pg.fill("#tf-acc", "P01375")
        pg.fill("#tf-orcid", "0000-0002-1825-0098")
        pg.click("button.tf-go")
        msgs = pg.inner_text("#tf-problems")
        opened = pg.evaluate("window.__opened.length")
        for want in ("ORCID iD is not valid", "first residue", "Pick what the experiment showed", "evidence code"):
            if want.lower() not in msgs.lower():
                print(f"FAIL bad answers: no message about '{want}'")
                bad += 1
        if opened:
            print("FAIL the form opened GitHub with bad answers")
            bad += 1
        if a.shots:
            pg.screenshot(path=str(a.shots / "truth_form_problems.png"), full_page=True)
        print(f"bad answers: {len(msgs.splitlines())} lines of messages, GitHub opened {opened} times")
        # 2. Good answers.
        for sel, val in GOOD.items():
            pg.fill(sel, val)
        pg.check('input[name="shown"][value="structure"]')
        pg.check("#tf-cc0")
        pg.click("button.tf-go")
        urls = pg.evaluate("window.__opened")
        browser.close()
    httpd.shutdown()
    if len(urls) != 1:
        print(f"FAIL good answers opened {len(urls)} pages")
        return 1
    u = urllib.parse.urlparse(urls[0])
    q = dict(urllib.parse.parse_qsl(u.query))
    print("opened:", f"{u.scheme}://{u.netloc}{u.path}")
    print("title:", q.get("title"), "| labels:", q.get("labels"))
    if (u.netloc, u.path) != ("github.com", f"/{tc.REPO}/issues/new"):
        print("FAIL the form opened the wrong page")
        bad += 1
    fields = tc.parse_body(q.get("body", ""))
    want = {"accession": "P01375", "start": "89", "end": "233", "shown": "structure", "eco": "ECO:0000269",
            "orcid": "0000-0002-1825-0097", "citation": GOOD["#tf-cite"], "sequence": "", "organism": ""}
    for k, v in want.items():
        if fields.get(k) != v:
            print(f"FAIL body field {k}: {fields.get(k)!r}, expected {v!r}")
            bad += 1
    problems = tc.check(fields, {"length": 233})
    if problems:
        print("FAIL the build would refuse the form's own text:", problems)
        bad += 1
    print(f"body parses back to {len(fields)} fields; build checks: {problems or 'pass'}")
    if a.url_out:
        a.url_out.write_text(urls[0] + "\n")
    print("contribute checks", "FAILED" if bad else "passed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
