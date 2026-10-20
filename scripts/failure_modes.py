#!/usr/bin/env python3
"""Check failure_modes/*.yaml against their sources, and write the site's copy.

    python scripts/failure_modes.py --check            # fetch every source from GitHub; exit 1 on any problem
    python scripts/failure_modes.py --out docs         # write docs/data/failure_modes.json (no network)

An entry is a way a tool is known to fail, with the run that showed it. Every entry must
rest on a source a reader can open: a PR body, a PR comment, or a file at a commit. --check
fetches each source and fails when

  - a quote is not in its source, character for character (Markdown included);
  - a number in the finding or in `number` is not in one of the entry's quotes;
  - a number in `run` or `what_to_do` is not anywhere in the entry's sources;
  - a named PR or commit does not exist.

A number is a run of digits, with any decimal point, comma or underscore between digits and
any exponent ("1e-22", "99_232", "307.65"), that does not follow a letter: "cys2" and
"PF14554" are names, not numbers. A number matches in a source only where it is not part of a
longer number, so "222" does not match "2_222" and "0.02" does not match "0.026".

Sources are fetched with the GitHub REST API. A token is read from GITHUB_TOKEN or GH_TOKEN,
or from `gh auth token`; the three source repos are public, so without one the check still
runs, at GitHub's lower rate limit.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

import jsonschema
import yaml

ROOT = Path(__file__).resolve().parents[1]
FM_DIR = ROOT / "failure_modes"
NUM = re.compile(r"(?<![\w.])\d+(?:[.,_]\d+)*(?:e-?\d+)?")


def numbers(text: str) -> list[str]:
    return NUM.findall(text or "")


def has_number(text: str, n: str) -> bool:
    """n appears in text, not as part of a longer number."""
    return re.search(r"(?<![\d.,_])" + re.escape(n) + r"(?![.,_]?\d)", text) is not None


def source_url(s: dict) -> str:
    """The page a reader opens for a source."""
    if s["kind"] == "pr_body":
        return f"https://github.com/{s['repo']}/pull/{s['pr']}"
    if s["kind"] == "pr_comment":
        return f"https://github.com/{s['repo']}/pull/{s['pr']}#issuecomment-{s['comment_id']}"
    return f"https://github.com/{s['repo']}/blob/{s['commit']}/{s['path']}"


def source_words(s: dict, site_repo: str) -> str:
    """How the page names a source: "PR 121 comment", "seanome/kmerseek PR 141", a file path."""
    repo = "" if s["repo"] == site_repo else f"{s['repo']} "
    if s["kind"] == "pr_body":
        return f"{repo}PR {s['pr']}"
    if s["kind"] == "pr_comment":
        return f"{repo}PR {s['pr']}, comment {s['comment_id']}"
    return f"{repo}{s['path']} at {s['commit'][:8]}"


class GitHub:
    """Fetches source texts once each."""

    def __init__(self):
        tok = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
        if not tok:
            try:
                tok = subprocess.run(["gh", "auth", "token"], capture_output=True, text=True, timeout=20).stdout.strip()
            except (OSError, subprocess.SubprocessError):
                tok = ""
        self.token = tok
        self.cache: dict[str, object] = {}

    def get(self, path: str, raw: bool = False):
        key = f"{raw}:{path}"
        if key not in self.cache:
            req = urllib.request.Request(f"https://api.github.com/{path}")
            req.add_header("Accept", "application/vnd.github.raw" if raw else "application/vnd.github+json")
            req.add_header("X-GitHub-Api-Version", "2022-11-28")
            if self.token:
                req.add_header("Authorization", f"Bearer {self.token}")
            last = None
            for _ in range(3):
                try:
                    with urllib.request.urlopen(req, timeout=60) as r:
                        body = r.read().decode("utf-8")
                    self.cache[key] = body if raw else json.loads(body)
                    break
                except urllib.error.HTTPError as e:
                    if e.code == 404:
                        raise LookupError(f"not found on GitHub: {path}") from e
                    last = e
                except (urllib.error.URLError, TimeoutError) as e:
                    last = e
            else:
                raise ConnectionError(f"could not fetch {path}: {last}")
        return self.cache[key]

    def text(self, s: dict) -> str:
        if s["kind"] == "pr_body":
            return self.get(f"repos/{s['repo']}/pulls/{s['pr']}")["body"] or ""
        if s["kind"] == "pr_comment":
            c = self.get(f"repos/{s['repo']}/issues/comments/{s['comment_id']}")
            if not c["issue_url"].endswith(f"/issues/{s['pr']}"):
                raise ValueError(f"comment {s['comment_id']} is not on PR {s['pr']} of {s['repo']}")
            return c["body"] or ""
        return self.get(f"repos/{s['repo']}/contents/{s['path']}?ref={s['commit']}", raw=True)

    def pr_exists(self, repo: str, n: int) -> None:
        self.get(f"repos/{repo}/pulls/{n}")

    def commit_exists(self, repo: str, sha: str) -> None:
        self.get(f"repos/{repo}/commits/{sha}")


def load(fm_dir: Path = FM_DIR) -> list[dict]:
    """Every failure_modes/<tool>.yaml, checked against the schema. No network."""
    schema = json.loads((fm_dir / "schema.json").read_text())
    docs = []
    for f in sorted(fm_dir.glob("*.yaml")):
        doc = yaml.safe_load(f.read_text())
        try:
            jsonschema.validate(doc, schema)
        except jsonschema.ValidationError as e:
            where = "/".join(map(str, e.absolute_path)) or "top level"
            raise ValueError(f"{f.name}: does not match failure_modes/schema.json at {where}: {e.message}") from e
        if doc["tool"] != f.stem:
            raise ValueError(f"{f.name}: tool is {doc['tool']!r}, but the file name says {f.stem!r}")
        ids = [e["id"] for e in doc["entries"]]
        if len(ids) != len(set(ids)):
            raise ValueError(f"{f.name}: entry ids repeat")
        doc["file"] = f"failure_modes/{f.name}"
        docs.append(doc)
    return docs


def check_entry(e: dict, gh) -> list[str]:
    """Problems with one entry; an empty list when it holds."""
    bad = []
    texts, quotes = [], []
    for ev in e["evidence"]:
        s = ev["source"]
        try:
            t = gh.text(s)
        except (LookupError, ValueError, ConnectionError) as x:
            bad.append(f"{source_url(s)}: {x}")
            continue
        # GitHub stores bodies with CRLF line ends when they were typed in the browser.
        t = t.replace("\r\n", "\n")
        texts.append(t)
        if ev["quote"] in t:
            quotes.append(ev["quote"])
        else:
            bad.append(f"quote not in {source_url(s)}: {ev['quote'][:90]!r}")
    q = "\n".join(quotes)
    full = "\n".join(texts)
    for n in numbers(e["finding"]):
        if not has_number(q, n):
            bad.append(f"finding number {n!r} is in no quote")
    if e["number"] is not None:
        for n in numbers(e["number"]):
            if not has_number(q, n):
                bad.append(f"number {n!r} is in no quote")
            if not has_number(e["finding"], n):
                bad.append(f"number {n!r} is not in the finding")
        if not numbers(e["number"]):
            bad.append(f"number {e['number']!r} holds no number")
    elif not e.get("number_note"):
        bad.append("number is null but number_note does not say why")
    for field in ("run", "what_to_do"):
        for n in numbers(e[field]):
            if not has_number(full, n):
                bad.append(f"{field} number {n!r} is in none of the sources")
    for p in e["pr"]:
        try:
            gh.pr_exists(p["repo"], p["number"])
        except (LookupError, ConnectionError) as x:
            bad.append(f"PR {p['repo']} {p['number']}: {x}")
    if "commit" in e:
        try:
            gh.commit_exists(e["commit"]["repo"], e["commit"]["sha"])
        except (LookupError, ConnectionError) as x:
            bad.append(f"commit {e['commit']['repo']}@{e['commit']['sha']}: {x}")
    return bad


def check(docs: list[dict], gh, known_tools: set[str]) -> int:
    n_bad = 0
    for d in docs:
        for e in d["entries"]:
            bad = check_entry(e, gh)
            for t in e.get("also_tools", []):
                if t not in known_tools:
                    bad.append(f"also_tools names {t!r}, which has no tool card")
            print(("ok   " if not bad else "FAIL ") + f"{d['file']} {e['id']}")
            for b in bad:
                print(f"       {b}")
            n_bad += bool(bad)
    print(f"{sum(len(d['entries']) for d in docs)} entries in {len(docs)} files, {n_bad} failing")
    return n_bad


def site_json(docs: list[dict], site_repo: str, tool_names: dict) -> dict:
    """The JSON twin of failure-modes.html: every entry, flat, with links resolved."""
    rows = []
    for d in docs:
        for e in d["entries"]:
            tools = [d["tool"]] + e.get("also_tools", [])
            rows.append({
                # A tool with a card keeps the card's name, so one tool has one name on the site.
                "tool": d["tool"], "tool_name": tool_names.get(d["tool"], d["name"]), "file": d["file"], "id": e["id"],
                "shown_on_cards": [t for t in tools if t in tool_names],
                "finding": e["finding"], "number": e["number"], "number_note": e.get("number_note"),
                "evidence": [{"source": source_words(ev["source"], site_repo), "url": source_url(ev["source"]),
                              "quote": ev["quote"]} for ev in e["evidence"]],
                "run": e["run"],
                "commit": ({**e["commit"], "url": f"https://github.com/{e['commit']['repo']}/commit/{e['commit']['sha']}"}
                           if "commit" in e else None),
                "pr": [{**p, "url": f"https://github.com/{p['repo']}/pull/{p['number']}"} for p in e["pr"]],
                "what_to_do": e["what_to_do"], "date": e["date"], "status": e["status"],
                "fixed": e["status"] != "open",
            })
    rows.sort(key=lambda r: (r["tool"], r["date"], r["id"]))
    return {"note": "Known failure modes, one row per entry of failure_modes/<tool>.yaml. Each quote "
                    "was checked against its source by scripts/failure_modes.py --check.",
            "schema": "failure_modes/schema.json", "entries": rows}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="fetch every source from GitHub and check it")
    ap.add_argument("--out", type=Path, default=None, help="site folder to write data/failure_modes.json into")
    ap.add_argument("--dir", type=Path, default=FM_DIR, help="folder of <tool>.yaml files")
    a = ap.parse_args()
    sys.path.insert(0, str(ROOT / "scripts"))
    from tools import TOOLS  # noqa: E402
    try:
        docs = load(a.dir)
    except ValueError as e:
        print(e)
        return 1
    rc = 0
    if a.check:
        rc = 1 if check(docs, GitHub(), set(TOOLS)) else 0
    if a.out:
        (a.out / "data").mkdir(parents=True, exist_ok=True)
        out = site_json(docs, "seanome/protein-function-annotation-leaderboard", {t: v["name"] for t, v in TOOLS.items()})
        (a.out / "data/failure_modes.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
