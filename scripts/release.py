#!/usr/bin/env python3
"""Freeze one version of the leaderboard: a Zenodo record with a DOI, mirrored to a Hugging
Face dataset (section 7 of docs/annotation_leaderboard_research.md, PR 147).

    python scripts/release.py --version v1 --dry-run   # builds the files, prints every upload
    ZENODO_TOKEN=... HF_TOKEN=... python scripts/release.py --version v1

Steps:
1. Build build/release/<version>/: results.parquet in the shape of schema/results_row.json
   (scripts/release_rows.py, checked by scripts/validate_rows.py), the cost table, the
   pile-up tables, every tool's resolved claims, the schema, SOURCE.json, a dataset card and
   MANIFEST.json with each file's size and sha256.
2. Zenodo: a new record, or a new version of the record named in releases/zenodo.json; upload
   the files; set the metadata; publish. Publishing gives the DOI. A published Zenodo record's
   files can never change, which is what freezes the version.
3. Hugging Face: upload the same files to the dataset repo, with the DOI in the card, and tag
   the commit with the version.
4. Write releases/<version>.json (DOI, record URL, dataset commit, file checksums) to commit.

The tokens are read from the environment and never printed. Without --dry-run the script
refuses to run when the git tree has uncommitted changes or the version was released before.
Set ZENODO_API=https://sandbox.zenodo.org/api to try the whole path on Zenodo's test server.
"""

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

ZENODO_API = os.environ.get("ZENODO_API", "https://zenodo.org/api")
HF_REPO = os.environ.get("HF_DATASET", "seanome/protein-function-annotation-leaderboard")
SITE = "https://seanome.github.io/protein-function-annotation-leaderboard/"
GITHUB = "https://github.com/seanome/protein-function-annotation-leaderboard"
# Zenodo's limits per record, from section 7 of the research doc (its references R27, R28).
ZENODO_MAX_FILES, ZENODO_MAX_BYTES = 100, 50 * 10**9


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=True).stdout.strip()


def card(version: str, date: str, manifest: dict, doi: str | None) -> str:
    """The dataset card: README.md on Hugging Face, and a file of the Zenodo record."""
    files = "\n".join(f"| `{f['name']}` | {f['bytes']:,} | {f['what']} |" for f in manifest["files"])
    cite = f"https://doi.org/{doi}" if doi else "the DOI is added when the Zenodo record is published"
    return f"""---
license: cc-by-4.0
pretty_name: Protein function annotation leaderboard {version}
configs:
- config_name: results
  data_files: results.parquet
---

# Protein function annotation leaderboard, {version} ({date})

Frozen tables behind {SITE} at commit `{manifest['site_commit'][:12]}` of {GITHUB}.
Scores come from seanome/2024-kmerseek-analysis at commit `{manifest['source']['commit'][:12]}`.
Cite this version: {cite}.

Each row of `results.parquet` is one score for one tool, setting, relative, stratum and metric,
with the cost of the run behind it. The columns are defined in `results_row.schema.json`.
Cost columns are measured from Nextflow traces; `usd_per_1000_from_rate` is those hours times
the price-list rates in the columns beside it, not a bill. Rows from runs made before tools were
submitted as containers have no container digest, 95% interval or run date.

| File | Bytes | What it holds |
| --- | ---: | --- |
{files}
"""


def build(version: str, date: str, out: Path) -> dict:
    import polars as pl
    import release_rows
    import validate_rows

    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = release_rows.build(ROOT / "data")
    rows.write_parquet(out / "results.parquet", compression="zstd", statistics=False)
    schema = json.loads((ROOT / "schema/results_row.json").read_text())
    errors, _ = validate_rows.validate(pl.read_parquet(out / "results.parquet"), schema, submitted=False)
    if errors:
        raise SystemExit("results.parquet fails schema/results_row.json:\n" + "\n".join(errors))
    what = {"results.parquet": f"{rows.height:,} score rows, schema/results_row.json version {schema['version']}"}
    copies = [("data/compute.parquet", "compute.parquet", "cost of each tool, setting and step, from the Nextflow traces"),
              ("data/pileup/regions.parquet", "pileup_regions.parquet", "truth regions on the region pile-up page"),
              ("data/pileup/calls.parquet", "pileup_calls.parquet", "tool calls on the region pile-up page"),
              ("data/SOURCE.json", "SOURCE.json", "commit and checksums of the analysis-repo tables"),
              ("schema/results_row.json", "results_row.schema.json", "columns of results.parquet")]
    for src, name, w in copies:
        shutil.copyfile(ROOT / src, out / name)
        what[name] = w
    claims_dir = ROOT / "docs/data/claims"
    for f in sorted(claims_dir.glob("*.json")):
        name = f"claims_{f.stem}.json"
        shutil.copyfile(f, out / name)
        what[name] = f"claims about {f.stem} and the verdict of each test"
    src = json.loads((ROOT / "data/SOURCE.json").read_text())
    manifest = {"version": version, "date": date, "site_commit": git("rev-parse", "HEAD"),
                "source": {k: src[k] for k in ("repo", "commit", "commit_date")},
                "files": [{"name": n, "bytes": (out / n).stat().st_size, "sha256": sha256(out / n), "what": what[n]}
                          for n in sorted(what)]}
    (out / "README.md").write_text(card(version, date, manifest, None))
    return manifest


def zenodo_metadata(version: str, date: str) -> dict:
    return {"metadata": {
        "title": f"Protein function annotation leaderboard, {version}",
        "upload_type": "dataset",
        "publication_date": date,
        "version": version,
        "creators": [{"name": "Botvinnik, Olga"}],
        "description": (f"Frozen tables behind the protein function annotation leaderboard ({SITE}): "
                        "scores per tool, setting, relative and stratum, with measured compute cost. "
                        "See README.md for the files and results_row.schema.json for the columns."),
        "license": "cc-by-4.0",
        "related_identifiers": [{"identifier": GITHUB, "relation": "isSupplementTo", "scheme": "url"}],
    }}


def http(method: str, url: str, token: str, body: bytes | None = None, ctype: str = "application/json") -> dict:
    req = urllib.request.Request(url, data=body, method=method,
                                 headers={"Authorization": f"Bearer {token}", "Content-Type": ctype})
    with urllib.request.urlopen(req, timeout=600) as r:
        txt = r.read()
    return json.loads(txt) if txt else {}


def zenodo(version: str, date: str, out: Path, files: list[str], token: str) -> dict:
    state_file = ROOT / "releases/zenodo.json"
    state = json.loads(state_file.read_text()) if state_file.exists() else {}
    if state.get("latest_record_id"):
        nv = http("POST", f"{ZENODO_API}/deposit/depositions/{state['latest_record_id']}/actions/newversion", token)
        dep = http("GET", nv["links"]["latest_draft"], token)
        for f in dep.get("files", []):  # a new version starts with the old files; drop them
            http("DELETE", f["links"]["self"], token)
    else:
        dep = http("POST", f"{ZENODO_API}/deposit/depositions", token, b"{}")
    bucket = dep["links"]["bucket"]
    for name in files:
        http("PUT", f"{bucket}/{name}", token, (out / name).read_bytes(), "application/octet-stream")
    http("PUT", dep["links"]["self"], token, json.dumps(zenodo_metadata(version, date)).encode())
    pub = http("POST", f"{dep['links']['self']}/actions/publish", token)
    state = {"concept_record_id": pub.get("conceptrecid"), "latest_record_id": pub["id"]}
    state_file.parent.mkdir(exist_ok=True)
    state_file.write_text(json.dumps(state, indent=1) + "\n")
    return {"doi": pub["doi"], "record_url": pub["links"]["html"], "record_id": pub["id"]}


def huggingface(version: str, out: Path, doi: str, token: str) -> str:
    from huggingface_hub import HfApi

    api = HfApi(token=token)
    api.create_repo(HF_REPO, repo_type="dataset", exist_ok=True)
    info = api.upload_folder(folder_path=str(out), repo_id=HF_REPO, repo_type="dataset",
                             commit_message=f"Leaderboard {version}, https://doi.org/{doi}")
    api.create_tag(HF_REPO, repo_type="dataset", tag=version, revision=info.oid,
                   tag_message=f"https://doi.org/{doi}")
    return info.oid


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", required=True, help="e.g. v1")
    ap.add_argument("--date", default=dt.date.today().isoformat())
    ap.add_argument("--out", type=Path, default=ROOT / "build/release")
    ap.add_argument("--dry-run", action="store_true", help="build the files and print every upload; send nothing")
    a = ap.parse_args()
    out = a.out / a.version
    record = ROOT / f"releases/{a.version}.json"
    if not a.dry_run:
        if git("status", "--porcelain"):
            raise SystemExit("the git tree has uncommitted changes; commit them so the release names one commit")
        if record.exists():
            raise SystemExit(f"{record.relative_to(ROOT)} exists: {a.version} was released already")
    manifest = build(a.version, a.date, out)
    files = [f["name"] for f in manifest["files"]] + ["README.md"]
    total = sum((out / n).stat().st_size for n in files)
    if len(files) > ZENODO_MAX_FILES or total > ZENODO_MAX_BYTES:
        raise SystemExit(f"{len(files)} files, {total:,} bytes: over Zenodo's {ZENODO_MAX_FILES} files or 50 GB")
    print(f"release {a.version} ({a.date}) from site commit {manifest['site_commit'][:12]}, "
          f"analysis commit {manifest['source']['commit'][:12]}")
    print(f"{len(files)} files, {total:,} bytes, in {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out}")
    for f in manifest["files"]:
        print(f"  {f['name']:32s} {f['bytes']:>12,}  sha256 {f['sha256'][:16]}")
    print(f"  {'README.md':32s} {(out / 'README.md').stat().st_size:>12,}")
    state_file = ROOT / "releases/zenodo.json"
    state = json.loads(state_file.read_text()) if state_file.exists() else {}
    zt, ht = os.environ.get("ZENODO_TOKEN"), os.environ.get("HF_TOKEN")
    if a.dry_run:
        print("\ndry run: nothing is sent. Without --dry-run it would:")
        if state.get("latest_record_id"):
            print(f"  POST {ZENODO_API}/deposit/depositions/{state['latest_record_id']}/actions/newversion"
                  " and DELETE the files carried over from the last version")
        else:
            print(f"  POST {ZENODO_API}/deposit/depositions   (a new record; releases/zenodo.json does not exist)")
        for n in files:
            print(f"  PUT  <bucket>/{n}")
        md = zenodo_metadata(a.version, a.date)["metadata"]
        print(f"  PUT  <deposition> metadata: title {md['title']!r}, version {md['version']}, "
              f"date {md['publication_date']}, license {md['license']}, type {md['upload_type']}")
        print("  POST <deposition>/actions/publish   -> the DOI")
        print(f"  Hugging Face: create dataset {HF_REPO} if missing, upload the {len(files)} files with the DOI "
              f"in README.md, tag the commit {a.version}")
        print(f"  write releases/{a.version}.json and releases/zenodo.json to commit")
        print(f"ZENODO_TOKEN is {'set' if zt else 'NOT set'}; HF_TOKEN is {'set' if ht else 'NOT set'}")
        return 0
    if not (zt and ht):
        raise SystemExit("set ZENODO_TOKEN and HF_TOKEN in the environment")
    try:  # checked before Zenodo publishes, so a missing package cannot leave half a release
        import huggingface_hub  # noqa: F401
    except ImportError:
        raise SystemExit("pip install huggingface_hub==1.33.0 (the version this script was checked against)")
    z = zenodo(a.version, a.date, out, files, zt)
    print(f"Zenodo: {z['record_url']}  DOI {z['doi']}")
    (out / "README.md").write_text(card(a.version, a.date, manifest, z["doi"]))
    oid = huggingface(a.version, out, z["doi"], ht)
    print(f"Hugging Face: https://huggingface.co/datasets/{HF_REPO}/tree/{a.version}  commit {oid}")
    record.write_text(json.dumps({**manifest, "doi": z["doi"], "zenodo": z["record_url"],
                                  "huggingface": f"https://huggingface.co/datasets/{HF_REPO}", "hf_commit": oid},
                                 indent=1) + "\n")
    print(f"wrote {record.relative_to(ROOT)}; commit it with releases/zenodo.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
