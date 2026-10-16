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
OPTIONAL = ["inventory.md"]
PILEUP = ["pileup/regions.parquet", "pileup/calls.parquet", "pileup/provenance.json"]


def gh(*args: str, raw: bool = False) -> bytes:
    cmd = ["gh", "api"]
    if raw:
        cmd += ["-H", "Accept: application/vnd.github.raw"]
    return subprocess.run(cmd + list(args), check=True, capture_output=True).stdout


def check_listed(want: dict, checksums: dict, listing: str, prefix: str = "") -> None:
    """Every file a provenance listing names must have been copied, with the listed sha256."""
    for name, digest in want.items():
        got = checksums.get(prefix + name)
        if got is None:
            raise SystemExit(f"{listing} lists {prefix + name}, which was not copied "
                             f"(copied: {', '.join(sorted(checksums))})")
        if got != digest:
            raise SystemExit(f"{prefix + name}: sha256 {got} != {digest} in {listing}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--commit", required=True)
    ap.add_argument("--out", type=Path, default=Path("data"))
    ap.add_argument("--pileup-commit", default=None,
                    help="commit holding results/leaderboard/pileup/, when it is not --commit")
    a = ap.parse_args()

    commit = json.loads(gh(f"repos/{a.repo}/commits/{a.commit}"))
    sha = commit["sha"]
    a.out.mkdir(parents=True, exist_ok=True)
    checksums = {}
    for name in FILES:
        body = gh(f"repos/{a.repo}/contents/results/leaderboard/{name}?ref={sha}", raw=True)
        (a.out / name).write_bytes(body)
        checksums[name] = hashlib.sha256(body).hexdigest()

    for name in OPTIONAL:
        try:
            body = gh(f"repos/{a.repo}/contents/results/leaderboard/{name}?ref={sha}", raw=True)
        except subprocess.CalledProcessError:
            # Not in this commit: remove an older copy, or it would be credited to this commit.
            (a.out / name).unlink(missing_ok=True)
            continue
        (a.out / name).write_bytes(body)
        checksums[name] = hashlib.sha256(body).hexdigest()
    pileup_sha = json.loads(gh(f"repos/{a.repo}/commits/{a.pileup_commit}"))["sha"] if a.pileup_commit else sha
    for name in PILEUP:
        try:
            body = gh(f"repos/{a.repo}/contents/results/leaderboard/{name}?ref={pileup_sha}", raw=True)
        except subprocess.CalledProcessError:
            # A commit from before the pile-up tables existed: drop any older copy too.
            (a.out / name).unlink(missing_ok=True)
            continue
        (a.out / name).parent.mkdir(parents=True, exist_ok=True)
        (a.out / name).write_bytes(body)
        checksums[name] = hashlib.sha256(body).hexdigest()
    got_pileup = [n for n in PILEUP if n in checksums]
    if got_pileup and len(got_pileup) != len(PILEUP):
        raise SystemExit(f"{pileup_sha[:8]} has only some of the pile-up files ({got_pileup}); "
                         "the page needs all of them from one commit")
    if "pileup/provenance.json" in checksums:
        want = json.loads((a.out / "pileup/provenance.json").read_text())["outputs"]
        check_listed(want, checksums, "pileup/provenance.json", prefix="pileup/")

    # The source's own provenance.json lists the checksum of each table it wrote.
    expected = json.loads((a.out / "provenance.json").read_text())["outputs"]
    check_listed(expected, checksums, "provenance.json")

    source = {
        "repo": a.repo,
        "commit": sha,
        "commit_date": commit["commit"]["committer"]["date"],
        "path": "results/leaderboard",
        "copied_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sha256": checksums,
        "pileup_commit": pileup_sha if "pileup/regions.parquet" in checksums else None,
    }
    (a.out / "SOURCE.json").write_text(json.dumps(source, indent=1) + "\n")
    print(f"copied {len(FILES)} files from {a.repo}@{sha[:8]}")


if __name__ == "__main__":
    main()
