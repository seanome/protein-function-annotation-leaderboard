"""The residues of the four named proteins on the why page, from the UniProt entries that
scripts/why/fetch.py downloaded (Prompt E). Each file's sha256 must match results/why/sources.json.

    python3 -I scripts/why/four_sequences.py --raw ~/data/why-numbers/2026_03
Writes data/why/four_sequences.json.
"""

import argparse
import hashlib
import json
import pathlib

FOUR = ["Q25460", "P24856", "P05484", "P42212"]
KEEP = ("Signal", "Propeptide", "Peptide", "Chain", "Repeat")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=pathlib.Path)
    ap.add_argument("--sources", type=pathlib.Path, default=pathlib.Path("results/why/sources.json"))
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("data/why/four_sequences.json"))
    a = ap.parse_args()
    sources = json.loads(a.sources.read_text())
    out = {}
    for acc in FOUR:
        src = sources[f"uniprot_{acc}"]
        p = a.raw / src["file"]
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
        if sha != src["sha256"]:
            raise SystemExit(f"{p}: sha256 {sha} is not the {src['sha256']} recorded in {a.sources}")
        d = json.loads(p.read_text())
        out[acc] = {
            "accession": acc,
            "name": d["proteinDescription"]["recommendedName"]["fullName"]["value"],
            "organism": d["organism"]["scientificName"],
            "common_name": d["organism"].get("commonName"),
            "length": d["sequence"]["length"],
            "sequence": d["sequence"]["value"],
            "features": [{"type": f["type"], "start": f["location"]["start"]["value"],
                          "end": f["location"]["end"]["value"], "description": f.get("description", "")}
                         for f in d.get("features", []) if f["type"] in KEEP],
            "source": {"url": src["url"], "sha256": sha, "release": src["release"], "fetched_utc": src["fetched_utc"]},
        }
        assert len(out[acc]["sequence"]) == out[acc]["length"]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print({k: v["length"] for k, v in out.items()})


if __name__ == "__main__":
    main()
