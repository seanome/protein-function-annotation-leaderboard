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
# The alphabet page's tables (results/leaderboard/alphabets/, `make leaderboard-alphabets`),
# plus PR 116's CSV from the same commit, which the page check compares square by square.
ALPHABETS = ["alphabets/classes.parquet", "alphabets/evidence.parquet", "alphabets/kappa.parquet",
             "alphabets/usable_k.parquet", "alphabets/k_reach.parquet", "alphabets/provenance.json"]
HP116_CSV = "tables/hp_borderline_residues.csv"


def copy_alphabets(repo: str, ref: str, out: Path) -> dict:
    """Copy the alphabet tables of one commit into out/alphabets/; return their sha256 and commit."""
    commit = json.loads(gh(f"repos/{repo}/commits/{ref}"))
    sha = commit["sha"]
    checksums = {}
    for name in ALPHABETS:
        body = gh(f"repos/{repo}/contents/results/leaderboard/{name}?ref={sha}", raw=True)
        (out / name).parent.mkdir(parents=True, exist_ok=True)
        (out / name).write_bytes(body)
        checksums[name] = hashlib.sha256(body).hexdigest()
    body = gh(f"repos/{repo}/contents/{HP116_CSV}?ref={sha}", raw=True)
    (out / "alphabets/hp_borderline_residues.csv").write_bytes(body)
    checksums["alphabets/hp_borderline_residues.csv"] = hashlib.sha256(body).hexdigest()
    prov = json.loads((out / "alphabets/provenance.json").read_text())
    check_listed(prov["outputs_sha256"], checksums, "alphabets/provenance.json", prefix="alphabets/")
    want_csv = prov["inputs_sha256"][HP116_CSV]
    if checksums["alphabets/hp_borderline_residues.csv"] != want_csv:
        raise SystemExit(f"{HP116_CSV} at {sha[:8]} is not the file the tables were built from")
    return {"alphabets_commit": sha, "alphabets_commit_date": commit["commit"]["committer"]["date"],
            "alphabets_copied_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "alphabets_sha256": checksums}


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
    ap.add_argument("--alphabets-commit", default=None,
                    help="commit holding results/leaderboard/alphabets/ (the alphabet page)")
    ap.add_argument("--alphabets-only", action="store_true",
                    help="copy only the alphabet tables and add them to the existing SOURCE.json")
    a = ap.parse_args()
    if a.alphabets_only:
        src = json.loads((a.out / "SOURCE.json").read_text())
        src.update(copy_alphabets(a.repo, a.alphabets_commit or a.commit, a.out))
        (a.out / "SOURCE.json").write_text(json.dumps(src, indent=1) + "\n")
        print(f"copied the alphabet tables from {a.repo}@{src['alphabets_commit'][:8]}")
        return

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
    if a.alphabets_commit:
        source.update(copy_alphabets(a.repo, a.alphabets_commit, a.out))
    (a.out / "SOURCE.json").write_text(json.dumps(source, indent=1) + "\n")
    print(f"copied {len(FILES)} files from {a.repo}@{sha[:8]}")


if __name__ == "__main__":
    main()
