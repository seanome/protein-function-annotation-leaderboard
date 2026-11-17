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
# The no-search baseline rows (scripts/build_no_search_baseline_table.py in the source repo).
# Optional, and copied from their own commit and path: the first copy is the laptop check in
# tables/no_search_baselines_laptop_check/, until the Sherlock run's table replaces it.
NO_SEARCH = ["no_search_baselines.parquet", "no_search_compute.parquet"]
# The best-hit transfer rows (source repo PR 160, nextflow-runs/qfo-pfam-region-benchmark
# bin/annotator_tool_rows.py and scripts/annotator_g1_call_rules.py). Copied from their own
# commit into data/best_hit/; each file must match the sha256 pinned in
# data/best_hit_pins.json, because the source folder has no provenance listing of its own.
# Paths from the source repo's root; each lands in data/best_hit/ under its file name
# (cards in data/best_hit/cards/).
BEST_HIT = [f"results/leaderboard/annotator_tools/{n}" for n in (
    "rows.parquet", "compute.parquet", "call_rules.parquet",
    "cards/diamond_besthit.sensitive.json",
    "cards/diamond_besthit.ultra-sensitive.json",
    "cards/mmseqs2_besthit.s7.5.json")] + [
    f"figures/annotator_tools/{n}_best_hit_vs_qfo_relatives_midi_plus.png"
    for n in ("fmax", "median_iou")]


def best_hit_dest(path: str) -> str:
    """data/best_hit/<name>, keeping cards/ for the cards."""
    parts = path.split("/")
    return "/".join(parts[-2:]) if parts[-2] == "cards" else parts[-1]


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
    ap.add_argument("--nosearch-commit", default=None,
                    help="commit holding the no-search baseline tables")
    ap.add_argument("--nosearch-path", default="results/leaderboard",
                    help="folder of the no-search baseline tables at --nosearch-commit")
    ap.add_argument("--besthit-commit", default=None,
                    help="commit holding the best-hit transfer tables")
    ap.add_argument("--besthit-pins", type=Path, default=Path("data/best_hit_pins.json"),
                    help="sha256 each best-hit file must have")
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

    nosearch_sha = None
    if a.nosearch_commit:
        nosearch_sha = json.loads(gh(f"repos/{a.repo}/commits/{a.nosearch_commit}"))["sha"]
    for name in NO_SEARCH:
        if nosearch_sha is None:
            (a.out / name).unlink(missing_ok=True)
            continue
        try:
            body = gh(f"repos/{a.repo}/contents/{a.nosearch_path}/{name}?ref={nosearch_sha}", raw=True)
        except subprocess.CalledProcessError:
            # The compute table only exists once the Sherlock run has a trace.
            (a.out / name).unlink(missing_ok=True)
            continue
        (a.out / name).write_bytes(body)
        checksums[name] = hashlib.sha256(body).hexdigest()
    if nosearch_sha and "no_search_baselines.parquet" not in checksums:
        raise SystemExit(f"{a.nosearch_path}/no_search_baselines.parquet is not at {nosearch_sha[:8]}")

    besthit_sha = None
    if a.besthit_commit:
        besthit_sha = json.loads(gh(f"repos/{a.repo}/commits/{a.besthit_commit}"))["sha"]
        pins = json.loads(a.besthit_pins.read_text())
        if pins.get("commit") != besthit_sha:
            raise SystemExit(f"{a.besthit_pins} pins commit {pins.get('commit')}, "
                             f"not {besthit_sha}: update the pins with the commit")
        for path in BEST_HIT:
            body = gh(f"repos/{a.repo}/contents/{path}?ref={besthit_sha}", raw=True)
            name = best_hit_dest(path)
            (a.out / "best_hit" / name).parent.mkdir(parents=True, exist_ok=True)
            (a.out / "best_hit" / name).write_bytes(body)
            checksums["best_hit/" + name] = hashlib.sha256(body).hexdigest()
        check_listed(pins["sha256"], checksums, str(a.besthit_pins), prefix="best_hit/")
    else:
        import shutil
        shutil.rmtree(a.out / "best_hit", ignore_errors=True)

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
        "nosearch_commit": nosearch_sha,
        "nosearch_path": a.nosearch_path if nosearch_sha else None,
        "besthit_commit": besthit_sha,
    }
    (a.out / "SOURCE.json").write_text(json.dumps(source, indent=1) + "\n")
    print(f"copied {len(FILES)} files from {a.repo}@{sha[:8]}")


if __name__ == "__main__":
    main()
