#!/usr/bin/env python3
"""Run the writing-style prose checker (check_prose.py) on every built page.

Fails on any cut-list word. Sentences over 30 words are listed for a person to read: the
checker joins list items and legend entries into one "sentence", so those need eyes, not a
hard failure.
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

DEFAULT = Path.home() / ".claude/skills/writing-style/scripts/check_prose.py"
# The checker needs beautifulsoup4, which the site venv does not carry.
PYTHON = os.environ.get("PROSE_PYTHON", "/Users/olga/anaconda3/envs/2025-kmerseek-analysis/bin/python3")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=Path("docs"))
    ap.add_argument("--checker", type=Path, default=Path(os.environ.get("PROSE_CHECK", DEFAULT)))
    a = ap.parse_args()
    if not a.checker.exists():
        print(f"prose checker not found at {a.checker}; set PROSE_CHECK")
        return 2
    bad = 0
    for page in sorted(a.site.rglob("*.html")):
        out = subprocess.run([PYTHON, str(a.checker), str(page)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout
        if "cut-list words" not in out:
            print(f"{page}: checker failed")
            bad += 1
            continue
        n_long = re.search(r"(\d+) sentences over", out)
        cut = re.search(r"cut-list words: (.*)", out)
        cut_s = cut.group(1).strip() if cut else "?"
        print(f"{page.relative_to(a.site)}: {n_long.group(1) if n_long else 0} long, cut-list words: {cut_s}")
        if cut_s != "none":
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
