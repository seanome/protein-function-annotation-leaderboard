"""Pick case 14: ten pathogen motifs from the Eukaryotic Linear Motif resource.

The rule, fixed before any search was run: of the ELM instances whose own protein is a virus or
a bacterium (by its UniProt lineage), which ELM records as experimentally shown ("true
positive"), and whose recorded binding partner is a human protein at the same motif
coordinates, take the ten with the lowest ELM instance identifiers.

    .venv/bin/python scripts/known_cases_elm.py --dir data/known_cases

Inputs, all kept in data/known_cases/ so the pick can be repeated without the network:
    elm_instances.tsv     ELM's instance table   (http://elm.eu.org/instances.tsv?q=*)
    elm_interactions.tsv  ELM's interaction table (http://elm.eu.org/interactions/as_tsv)
    elm_lineage.json      the UniProt lineage of every protein in the instance table
Output: elm_picked.json, which scripts/known_cases.py reads.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

UNIPROT = "https://rest.uniprot.org/uniprotkb/search?"


def lineages(accs: list[str], cache_file: Path, offline: bool) -> dict:
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
    todo = [a for a in accs if a not in cache]
    if todo and offline:
        raise SystemExit(f"{len(todo)} accessions are not in {cache_file} and --offline was given")
    for i in range(0, len(todo), 80):
        chunk = todo[i:i + 80]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        url = UNIPROT + urllib.parse.urlencode(
            {"query": q, "fields": "accession,organism_id,lineage", "format": "tsv", "size": 500})
        for line in urllib.request.urlopen(url, timeout=120).read().decode().splitlines()[1:]:
            f = (line.split("\t") + ["", ""])[:3]
            cache[f[0]] = {"taxon": f[1], "lineage": f[2]}
        time.sleep(0.3)
    cache_file.write_text(json.dumps(dict(sorted(cache.items())), indent=0) + "\n")
    return cache


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=Path("data/known_cases"))
    ap.add_argument("--offline", action="store_true")
    a = ap.parse_args()
    rows = list(csv.DictReader(
        (l for l in open(a.dir / "elm_instances.tsv") if not l.startswith("#")), delimiter="\t"))
    inter = list(csv.DictReader(open(a.dir / "elm_interactions.tsv"), delimiter="\t"))
    lin = lineages(sorted({r["Primary_Acc"] for r in rows}), a.dir / "elm_lineage.json", a.offline)

    by_acc: dict[str, list[dict]] = {}
    for i in inter:
        by_acc.setdefault(i["interactorElm"], []).append(i)

    def pathogen(acc: str) -> str:
        """Virus, bacterium, or neither. The lineage begins at "cellular organisms" for a
        bacterium, so every step of it is looked at, not only the first."""
        parts = [p.strip().split(" (")[0] for p in (lin.get(acc, {}).get("lineage") or "").split(",")]
        if "Viruses" in parts:
            return "virus"
        if "Bacteria" in parts:
            return "bacterium"
        return ""

    picked = []
    for r in sorted(rows, key=lambda r: r["Accession"]):
        if r["InstanceLogic"] != "true positive":
            continue
        kind = pathogen(r["Primary_Acc"])
        if not kind:
            continue
        hits = [i for i in by_acc.get(r["Primary_Acc"], [])
                if i["taxonomyDomain"].startswith("9606(")
                and i["Elm"] == r["ELMIdentifier"]
                and i["StartElm"] == r["Start"] and i["StopElm"] == r["End"]]
        if not hits:
            continue
        picked.append({**r, "superkingdom": kind, "taxon": lin[r["Primary_Acc"]]["taxon"],
                       "partner": hits[0]})
        if len(picked) == 10:
            break
    (a.dir / "elm_picked.json").write_text(json.dumps(picked, indent=1) + "\n")
    for p in picked:
        print(p["Accession"], p["ELMIdentifier"], p["Primary_Acc"], p["Start"], p["End"], "->",
              p["partner"]["interactorDomain"], "|", p["superkingdom"], "|", p["Organism"])
    print(f"wrote {a.dir/'elm_picked.json'}: {len(picked)} instances")


if __name__ == "__main__":
    main()
