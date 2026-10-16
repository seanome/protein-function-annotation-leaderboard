#!/usr/bin/env python3
"""Checks for failure_modes/ and the pages built from it.

Offline (always):
1. Every failure_modes/<tool>.yaml matches the schema, and its tool is its file name.
2. scripts/failure_modes.py's check, on made-up sources: a correct entry passes; a quote one
   character off fails; a number one digit off in the finding fails; a number not in the
   finding fails; "222" does not match "2_222"; "cys2" is a name, not a number.
3. The built site: failure-modes.html has a row for every entry; every tool card has a
   "Known failure modes" section, with its entries or "None recorded yet."; the JSON twins
   (data/failure_modes.json, failure-modes.html.json) equal what the YAML gives.

With --live (needs the network): the real kmerseek.yaml with 307.65 changed to 307.66, and
99_232 to 99_233, must fail the real check against GitHub.
"""

import argparse
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import failure_modes as fm  # noqa: E402
from tools import TOOLS  # noqa: E402

failures: list[str] = []


def check(cond: bool, what: str) -> None:
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        failures.append(what)


class FakeGitHub:
    """Serves made-up source texts, keyed the way the real class fetches them."""

    def __init__(self, texts: dict):
        self.texts = texts

    def text(self, s):
        return self.texts[(s["repo"], s.get("pr"), s.get("comment_id"))]

    def pr_exists(self, repo, n):
        pass

    def commit_exists(self, repo, sha):
        pass


SRC = {"kind": "pr_comment", "repo": "x/y", "pr": 1, "comment_id": 5}
GH = FakeGitHub({("x/y", 1, 5): "Intro.\r\nThe score stops at 307.65 in 102 of 107 settings; `cys2` lost 2_222 rows."})


def entry(**kw):
    e = {"id": "e", "finding": "The score stops at 307.65 in 102 of 107 settings.", "number": "307.65",
         "evidence": [{"source": SRC, "quote": "The score stops at 307.65 in 102 of 107 settings;"}],
         "run": "a run", "pr": [{"repo": "x/y", "number": 1}], "what_to_do": "Fix it.",
         "date": "2026-10-16", "status": "open"}
    e.update(kw)
    return e


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    ap.add_argument("--live", action="store_true", help="also run the wrong-number test against GitHub")
    a = ap.parse_args()

    # 1. schema
    docs = fm.load()
    check(len(docs) > 0, f"{len(docs)} failure_modes files load and match the schema")
    for d in docs:
        for e in d["entries"]:
            for t in e.get("also_tools", []):
                check(t in TOOLS, f"{d['tool']} {e['id']}: also_tools {t} has a card")

    # 2. the check itself, offline
    check(fm.check_entry(entry(), GH) == [], "a correct entry passes")
    bad = fm.check_entry(entry(evidence=[{"source": SRC, "quote": "The score stops at 307.65 in 102 of 107 setting;"}]), GH)
    check(any("quote not in" in b for b in bad), f"a quote one character off fails: {bad}")
    bad = fm.check_entry(entry(finding="The score stops at 307.66 in 102 of 107 settings.", number="307.66"), GH)
    check(any("'307.66'" in b for b in bad), f"a number one digit off fails: {bad}")
    bad = fm.check_entry(entry(finding="The score stops at 307.65 in 102 of 108 settings."), GH)
    check(any("'108'" in b for b in bad), f"a finding number not in the quote fails: {bad}")
    bad = fm.check_entry(entry(number="309"), GH)
    check(any("not in the finding" in b for b in bad), f"a number not in the finding fails: {bad}")
    check(not fm.has_number("lost 2_222 rows", "222") and fm.has_number("lost 2_222 rows", "2_222"),
          "222 does not match inside 2_222")
    check(not fm.has_number("E 0.026", "0.02"), "0.02 does not match inside 0.026")
    check(fm.numbers("under `cys2` and PF14554, v0.4.0, E 1e-22, 99_232") == ["1e-22", "99_232"],
          "names with digits are not numbers; exponents and underscores are kept")
    bad = fm.check_entry(entry(number=None), GH)
    check(any("number_note" in b for b in bad), "a null number needs a number_note")

    # 3. the built site
    want = fm.site_json(docs, "seanome/protein-function-annotation-leaderboard",
                        {t: v["name"] for t, v in TOOLS.items()
                         if (a.site / "tools" / f"{t}.html").exists()})
    got = json.loads((a.site / "data/failure_modes.json").read_text())
    check(got == want, "docs/data/failure_modes.json equals what the YAML gives")
    twin = json.loads((a.site / "failure-modes.html.json").read_text())["data"]
    check(twin == want, "failure-modes.html.json equals what the YAML gives")
    page = (a.site / "failure-modes.html").read_text()
    for r in want["entries"]:
        check(f'id="fm-{r["tool"]}-{r["id"]}"' in page, f"failure-modes.html has a row for {r['tool']} {r['id']}")
    for card in sorted((a.site / "tools").glob("*.html")):
        if card.name == "index.html":
            continue
        t = card.stem
        h = card.read_text()
        sec = re.search(r'<section class="cl-card fm-card" id="failure-modes".*?</section>', h, re.S)
        check(sec is not None, f"{card.name}: has a Known failure modes section")
        if not sec:
            continue
        mine = [r for r in want["entries"] if t in r["shown_on_cards"]]
        n_rows = sec.group(0).count('class="fm-row"')
        if mine:
            check(n_rows == len(mine), f"{card.name}: {n_rows} rows for {len(mine)} entries")
        else:
            check("None recorded yet." in sec.group(0) and n_rows == 0, f"{card.name}: says none recorded yet")
        md = (a.site / "tools" / f"{t}.md").read_text()
        check("## Known failure modes" in md and (mine or "None recorded yet." in md), f"{t}.md: has the section")

    # 4. live: a wrong number in a real file fails the real check
    if a.live:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            shutil.copy(ROOT / "failure_modes/schema.json", tmp / "schema.json")
            y = (ROOT / "failure_modes/kmerseek.yaml").read_text()
            y2 = y.replace("stops at 307.65, which", "stops at 307.66, which").replace("6 of 99_232 rows\n      showed", "6 of 99_233 rows\n      showed")
            check(y2.count("307.66") == 1 and y2.count("99_233") == 1, "live: made the two wrong numbers")
            (tmp / "kmerseek.yaml").write_text(y2)
            n_bad = fm.check(fm.load(tmp), fm.GitHub(), set(TOOLS))
            check(n_bad == 2, f"live: both entries with a wrong number fail against GitHub ({n_bad} failing)")

    if failures:
        print(f"{len(failures)} FAILURES")
        return 1
    print("all failure-mode checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
