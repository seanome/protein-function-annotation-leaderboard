#!/usr/bin/env python3
"""Make the Claude artifact copy of the built site.

The artifact viewer wraps the first page in its own <!doctype html><head><body>, so that
page is written without them: its <title>, stylesheet links and body content only. Every
other page is published beside it unchanged, with the same relative links.
"""

import re
import shutil
import sys
from pathlib import Path


def main() -> None:
    site, out = Path(sys.argv[1]), Path(sys.argv[2])
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(site, out)
    html = (site / "index.html").read_text()
    head = re.search(r"<head>(.*?)</head>", html, re.S).group(1)
    body = re.search(r"<body[^>]*>(.*)</body>", html, re.S).group(1)
    keep = [m for m in re.findall(r"<title>.*?</title>|<link[^>]+>", head)]
    if not any(k.startswith("<title>") for k in keep):
        raise SystemExit("no <title> in index.html")
    (out / "index.html").write_text("\n".join(keep) + "\n" + body)
    print(f"artifact copy in {out}: {sum(1 for _ in out.rglob('*') if _.is_file())} files")


if __name__ == "__main__":
    main()
