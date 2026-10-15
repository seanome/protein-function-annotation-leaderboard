#!/usr/bin/env python3
"""Copy results/leaderboard/ from one commit of the source repo into data/.

Uses `gh api` so the files come from the commit on GitHub, not from a local checkout
that may hold uncommitted changes. Writes data/SOURCE.json with the full commit hash,
which every page prints.
"""

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

FILES = ["results.parquet", "compute.parquet", "coverage.md", "provenance.json"]


def gh(*args: str, raw: bool = False) -> bytes:
    cmd = ["gh", "api"]
    if raw:
        cmd += ["-H", "Accept: application/vnd.github.raw"]
    return subprocess.run(cmd + list(args), check=True, capture_output=True).stdout


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--commit", required=True)
    ap.add_argument("--out", type=Path, default=Path("data"))
    a = ap.parse_args()

    commit = json.loads(gh(f"repos/{a.repo}/commits/{a.commit}"))
    sha = commit["sha"]
    a.out.mkdir(parents=True, exist_ok=True)
    checksums = {}
    for name in FILES:
        body = gh(f"repos/{a.repo}/contents/results/leaderboard/{name}?ref={sha}", raw=True)
        (a.out / name).write_bytes(body)
        checksums[name] = hashlib.sha256(body).hexdigest()

    # The source's own provenance.json lists the checksum of each table it wrote.
    expected = json.loads((a.out / "provenance.json").read_text())["outputs"]
    for name, digest in expected.items():
        if checksums[name] != digest:
            raise SystemExit(f"{name}: sha256 {checksums[name]} != {digest} in provenance.json")

    source = {
        "repo": a.repo,
        "commit": sha,
        "commit_date": commit["commit"]["committer"]["date"],
        "path": "results/leaderboard",
        "copied_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sha256": checksums,
    }
    (a.out / "SOURCE.json").write_text(json.dumps(source, indent=1) + "\n")
    print(f"copied {len(FILES)} files from {a.repo}@{sha[:8]}")


if __name__ == "__main__":
    main()
