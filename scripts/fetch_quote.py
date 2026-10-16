#!/usr/bin/env python3
"""Fetch a tool author's page and record a quote from it for claims/<tool>.yaml.

    python scripts/fetch_quote.py --text-url https://raw.githubusercontent.com/.../README.md \\
        --url https://github.com/<owner>/<repo> --quote "words copied from the page"
        # prints the author_claim block, or exits 1 when the quote is not on the page

    python scripts/fetch_quote.py --check-all          # fetch every author_claim again

A quote matches when its words, case and punctuation appear on the page in the same order,
with runs of whitespace (line breaks included) treated as one space and HTML tags removed.
The printed sha256 is what scripts/claims.py checks on every build, so a quote edited by hand
fails the build until this script has found the new words on the page.
"""

import argparse
import datetime as dt
import html
import re
import sys
import urllib.request
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from claims import normalise, quote_sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def page_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "protein-function-annotation-leaderboard"})
    with urllib.request.urlopen(req, timeout=60) as r:
        raw = r.read().decode("utf-8", errors="replace")
    if "<html" in raw.lower()[:2000]:
        raw = re.sub(r"(?is)<(script|style).*?</\1>", " ", raw)
        raw = html.unescape(re.sub(r"<[^>]+>", " ", raw))
    return normalise(raw)


def on_page(quote: str, text: str) -> bool:
    return normalise(quote) in text


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--text-url")
    ap.add_argument("--url")
    ap.add_argument("--quote")
    ap.add_argument("--check-all", action="store_true")
    a = ap.parse_args()
    if a.check_all:
        bad = 0
        for f in sorted((ROOT / "claims").glob("*.yaml")):
            for c in yaml.safe_load(f.read_text())["claims"]:
                q = c.get("author_claim")
                if not q:
                    continue
                ok = on_page(q["quote"], page_text(q["text_url"]))
                bad += not ok
                print(f"{'ok  ' if ok else 'GONE'} {f.name} {c['id']}: \"{q['quote']}\" at {q['text_url']}")
        return 1 if bad else 0
    if not (a.text_url and a.url and a.quote):
        ap.error("--text-url, --url and --quote are needed")
    if not on_page(a.quote, page_text(a.text_url)):
        print(f"not on the page: {a.quote!r} at {a.text_url}", file=sys.stderr)
        return 1
    block = {"author_claim": {"quote": normalise(a.quote), "url": a.url, "text_url": a.text_url,
                              "fetched": dt.date.today().isoformat(), "sha256": quote_sha256(a.quote)}}
    print(yaml.safe_dump(block, sort_keys=False, allow_unicode=True, width=200), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
