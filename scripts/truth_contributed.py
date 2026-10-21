#!/usr/bin/env python3
"""Build the truth_contributed table from the curation queue.

The queue is the repository's GitHub issues labelled truth-submission (the /contribute/truth/
form and .github/ISSUE_TEMPLATE/truth.yml both write them). The curator's decision on each one
is a line in data/curation/decisions.tsv, added with

    scripts/truth_contributed.py decide --issue 12 --accept BINDING --reason "..." --curator <ORCID or GitHub name>
    scripts/truth_contributed.py decide --issue 13 --reject --reason "only restates UniProt FT BINDING" --curator ...

and mirrored on the issue as a comment with the same reason (docs/curation.md). The build

    scripts/truth_contributed.py build [--issues issues.json]

reads every accepted issue, checks it again, and writes

    data/truth_contributed.tsv               one row per accepted region, provenance = wet_lab
    data/truth_contributed/<query set>.parquet  the rows on that set's query proteins, in the
                                             scorer's truth layout, so scripts/score_submission.py
                                             grades every submission on them too
    docs/data/contributed_counts.json        rows per organism, for the organism page

Without --issues it reads the issues with `gh api`. Each protein's organism and length come
from the UniProt REST API, kept in data/curation/uniprot_cache.json so a rebuild is the same
offline.
"""

import argparse
import datetime as dt
import json
import re
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parent.parent
REPO = "seanome/protein-function-annotation-leaderboard"
LABEL = "truth-submission"
SHOWN = ["binding", "activity", "localisation", "structure"]
# The Swiss-Prot feature type the curator files an accepted row under. Search tools are
# graded on it by copying the relative's Swiss-Prot features of this type, so it must be one
# of the types the Swiss-Prot truth set already uses (scorer/build_swissprot_truth.py).
FT_TYPES = ["ACT_SITE", "BINDING", "COILED", "DNA_BIND", "DOMAIN", "INTRAMEM", "MOTIF",
            "REGION", "REPEAT", "SITE", "TRANSMEM", "ZN_FING"]
ACC = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
DECISIONS = ROOT / "data/curation/decisions.tsv"
CACHE = ROOT / "data/curation/uniprot_cache.json"
DECISION_COLS = ["issue", "decision", "swissprot_type", "reason", "curator", "decided_on"]
HEADINGS = {
    "UniProt accession": "accession", "Sequence": "sequence", "Organism": "organism",
    "First residue": "start", "Last residue": "end", "What was shown": "shown",
    "ECO evidence code": "eco", "Citation or note": "citation", "ORCID iD": "orcid",
    "Licence": "cc0", "Source": "newexp",
}
TABLE_COLS = ["issue", "issue_url", "accession", "sequence", "organism", "taxon_id", "start",
              "end", "shown", "swissprot_type", "eco", "citation", "orcid", "provenance",
              "curator", "decided_on", "reason"]


def orcid_ok(s: str) -> bool:
    """ISO 7064 MOD 11-2 check character, as ORCID defines it."""
    if not re.fullmatch(r"\d{4}-\d{4}-\d{4}-\d{3}[\dX]", s):
        return False
    d = s.replace("-", "")
    total = 0
    for c in d[:15]:
        total = (total + int(c)) * 2
    r = (12 - total % 11) % 11
    return d[15] == ("X" if r == 10 else str(r))


def parse_body(body: str) -> dict:
    """The '### heading' sections an issue form (or the site's form) writes."""
    out = {}
    for block in re.split(r"^### ", body or "", flags=re.M)[1:]:
        head, _, text = block.partition("\n")
        key = HEADINGS.get(head.strip())
        if key:
            text = text.strip()
            out[key] = "" if text == "_No response_" else text
    return out


def check(f: dict, uniprot: dict | None) -> list[str]:
    p = []
    acc, seq = f.get("accession", "").upper(), re.sub(r"\s+", "", f.get("sequence", "")).upper()
    if not acc and not seq:
        p.append("no accession and no sequence")
    if acc and not ACC.match(acc):
        p.append(f"'{acc}' is not a UniProt accession")
    if acc and uniprot is None:
        p.append(f"UniProt has no entry {acc}")
    if seq and not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYUOBZXJ*]+", seq):
        p.append("sequence holds characters that are not amino-acid codes")
    if not acc and not f.get("organism"):
        p.append("a sequence without an organism")
    try:
        s, e = int(f.get("start", "")), int(f.get("end", ""))
    except ValueError:
        return p + ["first or last residue is not a whole number"]
    length = (uniprot or {}).get("length") if acc else (len(seq) if seq else None)
    if s < 1 or e < s:
        p.append(f"region {s}-{e} is not a region (positions count from 1, end >= start)")
    if length and e > length:
        p.append(f"region ends at {e}, past the protein's {length} residues")
    if f.get("shown") not in SHOWN:
        p.append(f"'what was shown' must be one of {', '.join(SHOWN)}")
    if not re.fullmatch(r"ECO:\d{7}", f.get("eco", "").upper()):
        p.append("evidence code is not ECO:nnnnnnn")
    if not f.get("citation"):
        p.append("no citation or note")
    if not orcid_ok(f.get("orcid", "").upper()):
        p.append("ORCID iD fails its check digit")
    if "[x]" not in f.get("cc0", "").lower():
        p.append("CC0 box not ticked")
    if "[x]" not in f.get("newexp", "").lower():
        p.append("'comes from an experiment' box not ticked")
    return p


def uniprot_entry(acc: str, cache: dict) -> dict | None:
    if acc in cache:
        return cache[acc]
    url = (f"https://rest.uniprot.org/uniprotkb/{acc}?format=json"
           "&fields=accession,organism_name,organism_id,length")
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            j = json.load(r)
        entry = {"organism": j["organism"]["scientificName"], "taxon_id": j["organism"]["taxonId"],
                 "length": j["sequence"]["length"], "fetched": dt.date.today().isoformat()}
    except urllib.error.HTTPError as e:
        print(f"UniProt {acc}: {e}", file=sys.stderr)
        if e.code != 404:  # a real "no such entry" is cached; other failures are tried again next run
            return None
        entry = None
    except Exception as e:  # noqa: BLE001  network error or bad reply: unknown now, not cached
        print(f"UniProt {acc}: {e}", file=sys.stderr)
        return None
    cache[acc] = entry
    return entry


def read_issues(path: Path | None) -> list[dict]:
    if path:
        return json.loads(path.read_text())
    out = subprocess.run(["gh", "api", "--paginate", f"repos/{REPO}/issues?labels={LABEL}&state=all&per_page=100"],
                         check=True, capture_output=True, text=True).stdout
    # --paginate writes one JSON array per page, back to back.
    dec, pos, issues = json.JSONDecoder(), 0, []
    while pos < len(out):
        if out[pos].isspace():
            pos += 1
            continue
        page, pos = dec.raw_decode(out, pos)
        issues.extend(page)
    return issues


def read_decisions() -> pl.DataFrame:
    if not DECISIONS.exists():
        return pl.DataFrame(schema={c: pl.String for c in DECISION_COLS})
    return pl.read_csv(DECISIONS, separator="\t", infer_schema=False)


def decide(a) -> int:
    if a.accept and a.accept not in FT_TYPES:
        raise SystemExit(f"--accept takes a Swiss-Prot feature type: {', '.join(FT_TYPES)}")
    if not a.reason.strip():
        raise SystemExit("--reason is required: every decision carries one")
    d = read_decisions()
    if str(a.issue) in d["issue"].to_list():
        raise SystemExit(f"issue {a.issue} already has a decision in {DECISIONS}; edit that line instead")
    row = pl.DataFrame([{"issue": str(a.issue), "decision": "accepted" if a.accept else "rejected",
                         "swissprot_type": a.accept or "", "reason": a.reason.replace("\t", " ").strip(),
                         "curator": a.curator, "decided_on": a.date or dt.date.today().isoformat()}])
    DECISIONS.parent.mkdir(parents=True, exist_ok=True)
    pl.concat([d, row]).sort(pl.col("issue").cast(pl.Int64)).write_csv(DECISIONS, separator="\t")
    print(f"issue {a.issue}: {row['decision'][0]} -> {DECISIONS}")
    return 0


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths, acc = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            acc = name.split("|")[1] if "|" in name else name
            lengths[acc] = 0
        elif acc:
            lengths[acc] += len(line.strip())
    return lengths


def build(a) -> int:
    issues = {str(i["number"]): i for i in read_issues(a.issues) if "pull_request" not in i}
    decisions = read_decisions().filter(pl.col("decision") == "accepted")
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    rows, bad = [], 0
    for d in decisions.iter_rows(named=True):
        issue = issues.get(d["issue"])
        if issue is None:
            print(f"issue {d['issue']}: accepted in decisions.tsv but not in the queue", file=sys.stderr)
            bad += 1
            continue
        f = parse_body(issue["body"])
        acc = f.get("accession", "").upper()
        up = uniprot_entry(acc, cache) if acc else None
        problems = check(f, up)
        if problems:
            print(f"issue {d['issue']}: accepted but fails the checks: {'; '.join(problems)}", file=sys.stderr)
            bad += 1
            continue
        org = f.get("organism", "")
        rows.append({
            "issue": int(d["issue"]), "issue_url": issue["html_url"], "accession": acc,
            "sequence": re.sub(r"\s+", "", f.get("sequence", "")).upper(),
            "organism": up["organism"] if up else org,
            "taxon_id": up["taxon_id"] if up else (int(org) if org.isdigit() else None),
            "start": int(f["start"]), "end": int(f["end"]), "shown": f["shown"],
            "swissprot_type": d["swissprot_type"], "eco": f["eco"].upper(), "citation": f["citation"],
            "orcid": f["orcid"].upper(), "provenance": "wet_lab", "curator": d["curator"],
            "decided_on": d["decided_on"], "reason": d["reason"],
        })
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(dict(sorted(cache.items())), indent=1) + "\n")
    schema = {"issue": pl.Int64, "issue_url": pl.String, "accession": pl.String, "sequence": pl.String,
              "organism": pl.String, "taxon_id": pl.Int64, "start": pl.Int64, "end": pl.Int64,
              "shown": pl.String, "swissprot_type": pl.String, "eco": pl.String, "citation": pl.String,
              "orcid": pl.String, "provenance": pl.String, "curator": pl.String,
              "decided_on": pl.String, "reason": pl.String}
    table = pl.DataFrame(rows, schema=schema).sort("issue")
    # The fourteen known cases are the table's first entries. They come from published
    # experiments rather than from the issue queue, so they are kept in their own file and read
    # back in at every build; their issue number is 0 and their issue_url is the paper's DOI.
    seed = a.data / "truth_known_cases.tsv"
    if seed.exists():
        table = pl.concat([pl.read_csv(seed, separator="\t", schema=schema), table])
    out = a.data / "truth_contributed.tsv"
    table.write_csv(out, separator="\t")

    # Per query set: the rows on its query proteins, in the scorer's truth layout. The
    # scorer's interval convention is the Swiss-Prot truth's: start as given, end + 1 for a
    # one-residue feature (is_point), as build_swissprot_truth widens point features.
    tdir = a.data / "truth_contributed"
    tdir.mkdir(exist_ok=True)
    for qs in sorted(p for p in (a.data / "query_sets").iterdir() if p.is_dir()):
        lengths = fasta_lengths(qs / "queries.fasta")
        sub = table.filter(pl.col("accession").is_in(list(lengths)))
        target = tdir / f"{qs.name}.parquet"
        if sub.height == 0:
            target.unlink(missing_ok=True)
            continue
        sp = pl.read_parquet(qs / "truth/swissprot/human_truth.parquet")
        point = pl.col("start") == pl.col("end")
        sub.select(
            pl.col("accession"), pl.col("swissprot_type").alias("pfam_id"),
            pl.col("start").cast(sp["domain_start"].dtype).alias("domain_start"),
            pl.when(point).then(pl.col("end") + 1).otherwise(pl.col("end")).cast(sp["domain_end"].dtype).alias("domain_end"),
            pl.col("accession").replace_strict(lengths, return_dtype=sp["protein_length"].dtype).alias("protein_length"),
            point.alias("is_point"),
        ).write_parquet(target)
        print(f"{target}: {sub.height} rows")

    counts = (table.group_by("taxon_id", "organism").agg(pl.len().alias("n_rows"))
              .sort("n_rows", "organism", descending=[True, False]))
    payload = {
        "note": "Accepted wet-lab rows in data/truth_contributed.tsv, per organism. Keyed by NCBI taxon ID; "
                "rows given as a sequence with an organism name and no taxon ID are under 'by_name'.",
        "n_rows": table.height,
        "by_taxon": {str(r["taxon_id"]): {"organism": r["organism"], "n_rows": r["n_rows"]}
                     for r in counts.iter_rows(named=True) if r["taxon_id"] is not None},
        "by_name": {r["organism"]: r["n_rows"] for r in counts.iter_rows(named=True) if r["taxon_id"] is None},
    }
    a.counts.parent.mkdir(parents=True, exist_ok=True)
    a.counts.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"{out}: {table.height} accepted rows; {bad} accepted issues failed the checks")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--issues", type=Path, help="issues as JSON (GitHub API layout); default: read with gh api")
    b.add_argument("--data", type=Path, default=ROOT / "data")
    b.add_argument("--counts", type=Path, default=ROOT / "docs/data/contributed_counts.json")
    d = sub.add_parser("decide")
    d.add_argument("--issue", type=int, required=True)
    g = d.add_mutually_exclusive_group(required=True)
    g.add_argument("--accept", metavar="SWISSPROT_TYPE")
    g.add_argument("--reject", action="store_true")
    d.add_argument("--reason", required=True)
    d.add_argument("--curator", required=True)
    d.add_argument("--date")
    for p in (b, d):
        p.add_argument("--curation", type=Path, default=ROOT / "data/curation",
                       help="folder holding decisions.tsv and uniprot_cache.json")
    a = ap.parse_args()
    global DECISIONS, CACHE
    DECISIONS, CACHE = a.curation / "decisions.tsv", a.curation / "uniprot_cache.json"
    return build(a) if a.cmd == "build" else decide(a)


if __name__ == "__main__":
    sys.exit(main())
