"""Write the fourteen known cases as the first entries of the wet-lab truth table.

    .venv/bin/python scripts/known_cases_truth.py

Reads data/known_cases/cases.json and writes data/truth_known_cases.tsv in the schema of
data/truth_contributed.tsv, which scripts/truth_contributed.py reads back in at every build.
One row per case: the stretch of the query protein the cited experiment covers, the kind of
evidence, and the paper. The issue number is 0 and the issue_url is the paper's DOI, because
these rows come from published experiments and not from the curation queue.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

COLS = ["issue", "issue_url", "accession", "sequence", "organism", "taxon_id", "start", "end",
        "shown", "swissprot_type", "eco", "citation", "orcid", "provenance", "curator",
        "decided_on", "reason"]

# Per case: what the experiment shows (the form's four kinds), the Swiss-Prot feature type a
# curator would file it under, and the residues it covers. "whole" means the whole protein,
# because the experiment is about the protein and not about one stretch of it.
PER_CASE = {
    1: ("activity", "DOMAIN", "whole"),
    2: ("structure", "DOMAIN", "whole"),
    3: ("localisation", "DOMAIN", "whole"),
    4: ("structure", "DOMAIN", "whole"),
    5: ("structure", "DOMAIN", "whole"),
    6: ("structure", "DOMAIN", "whole"),
    7: ("structure", "DOMAIN", "whole"),
    8: ("activity", "DOMAIN", "whole"),
    9: ("binding", "REGION", "whole"),
    10: ("structure", "REGION", "whole"),
    11: ("structure", "REGION", "whole"),
    12: ("structure", "REPEAT", "whole"),
    13: ("binding", "MOTIF", "region"),
    14: ("binding", "MOTIF", "region"),
}
DECIDED_ON = "2026-10-17"
CURATOR = "known-cases page, seanome/protein-function-annotation-leaderboard"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, default=Path("data/known_cases/cases.json"))
    ap.add_argument("--out", type=Path, default=Path("data/truth_known_cases.tsv"))
    a = ap.parse_args()
    doc = json.loads(a.cases.read_text())
    refs = {r["n"]: r for r in doc["references"]}
    rows = []
    for c in doc["cases"]:
        shown, ft, where = PER_CASE[c["id"]]
        q = c["queries"][0]
        if where == "region":
            start, end = c.get("query_region") or q["region"]
        else:
            start, end = 1, q["length"]
        ref = refs[c["refs"][0]]
        extra = ""
        if len(c["queries"]) > 1:
            extra = (" The case also covers "
                     + ", ".join(x["acc"] for x in c["queries"][1:]) + ".")
        rows.append({
            "issue": 0,
            "issue_url": f'https://doi.org/{ref["doi"]}',
            "accession": q["acc"],
            "sequence": "",
            "organism": q["species"],
            "taxon_id": "",
            "start": start,
            "end": end,
            "shown": shown,
            "swissprot_type": ft,
            "eco": c["eco"],
            "citation": ref["text"].replace("<i>", "").replace("</i>", ""),
            "orcid": "",
            "provenance": "wet_lab",
            "curator": CURATOR,
            "decided_on": DECIDED_ON,
            "reason": f'Known case {c["id"]}: {c["name"]}. {c["experiment"]}{extra}',
        })
    lines = ["\t".join(COLS)]
    for r in rows:
        lines.append("\t".join(str(r[c]).replace("\t", " ") for c in COLS))
    a.out.write_text("\n".join(lines) + "\n")
    print(f"wrote {a.out}: {len(rows)} rows")


if __name__ == "__main__":
    main()
