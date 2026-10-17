"""Read the output of scripts/known_cases_run.sbatch into data/known_cases/cases.json.

Each tool writes one table of hits. This script keeps, for every case and tool, the best hit on
that case's target protein for the real query and for its decoy, works out what the hit means
(found, wrong place, wrong label, missed), and writes those into the cases file the page is
built from. A tool with no output file keeps its "not run" row and its reason.

    .venv/bin/python scripts/known_cases_results.py --run $SCRATCH/known-cases

The query names in the run come from data/known_cases/queries.fasta:
    case<N>|<accession>|real   or   case<N>|<accession>|decoy
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

IOU_FOUND = 0.5  # the same cutoff Table 1 uses
# tool -> (file name written by the run, how to read it)
FILES = {
    "hmmer3_phmmer": ("phmmer.domtbl", "domtbl"),
    "hmmer3_jackhmmer": ("jackhmmer.domtbl", "domtbl"),
    "hmmscan": ("hmmscan.domtbl", "domtbl"),
    "mmseqs2_seqseq": ("mmseqs2_seqseq.tsv", "blasttab"),
    "mmseqs2_iterative": ("mmseqs2_iterative.tsv", "blasttab"),
    "foldseek": ("foldseek.tsv", "blasttab"),
    "prostt5": ("prostt5.tsv", "blasttab"),
    "hhblits": ("hhblits.tsv", "blasttab"),
    "kmerseek": ("kmerseek.tsv", "kmerseek"),
}
ACC = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")


def iou(a: tuple[int, int], b: tuple[int, int]) -> float:
    """Overlap of two closed residue ranges divided by the residues either one covers."""
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    share = max(0, hi - lo + 1)
    union = (a[1] - a[0] + 1) + (b[1] - b[0] + 1) - share
    return share / union if union else 0.0


def parse_query(name: str) -> tuple[int, str, bool] | None:
    f = name.split("|")
    if len(f) < 3 or not f[0].startswith("case"):
        return None
    return int(f[0][4:]), f[1], f[2].startswith("decoy")


def target_acc(name: str) -> str:
    """The accession in a target name, whatever header style the database carries."""
    m = ACC.search(name.replace("|", " "))
    return m.group(0) if m else name


def read_hits(path: Path, kind: str) -> list[dict]:
    """One row per hit: query name, target accession, query start and end, score, E-value."""
    out = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t") if kind != "domtbl" else line.split()
        if kind == "domtbl":
            # HMMER's --domtblout: target name 0, query name 3, i-Evalue 12, domain score 13,
            # hmm from/to 15/16, ali from/to 17/18. hmmscan swaps query and target.
            if len(f) < 23:
                continue
            out.append({"query": f[3], "target": target_acc(f[0]),
                        "qstart": int(f[17]), "qend": int(f[18]),
                        "score": float(f[13]), "evalue": float(f[12]), "score_kind": "bits"})
        elif kind == "blasttab":
            if len(f) < 8:
                continue
            out.append({"query": f[0], "target": target_acc(f[1]),
                        "qstart": int(f[2]), "qend": int(f[3]),
                        "score": float(f[6]), "evalue": float(f[7]), "score_kind": "bits"})
        else:  # kmerseek: query, target, qstart, qend, score, evalue
            if len(f) < 6:
                continue
            out.append({"query": f[0], "target": target_acc(f[1]),
                        "qstart": int(f[2]), "qend": int(f[3]),
                        "score": float(f[4]), "evalue": float(f[5]),
                        "score_kind": "region enrichment"})
    return out


def best(hits: list[dict]) -> dict | None:
    return max(hits, key=lambda h: h["score"]) if hits else None


def classify(case: dict, hit: dict | None, decoy: dict | None) -> dict:
    """What the hit means for this case, by the rules of the case's box."""
    if hit is None:
        return {"result": "missed", "reason": "no hit on the target protein above the tool's "
                                              "own cutoff"}
    if decoy is not None and decoy["score"] >= hit["score"]:
        return {"result": "missed",
                "reason": "the shuffled copy of the query hits the same protein as well as the "
                          "real query does, so the hit is not counted"}
    truth = case.get("query_region")
    if truth:
        overlap = iou((hit["qstart"], hit["qend"]), tuple(truth))
        if overlap < IOU_FOUND:
            return {"result": "wrong place", "iou": overlap,
                    "reason": "the call is on the right protein but off the proven region"}
        return {"result": "found", "iou": overlap}
    # Box 3 has no family truth: naming the human protein's job is a wrong label.
    if case["box"] == 3:
        return {"result": "wrong label",
                "reason": "there is no shared ancestor here, so a call that names the human "
                          "protein's family carries a job it cannot inherit"}
    return {"result": "found"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True, help="the run's output directory")
    ap.add_argument("--cases", type=Path, default=Path("data/known_cases/cases.json"))
    a = ap.parse_args()
    doc = json.loads(a.cases.read_text())
    by_id = {c["id"]: c for c in doc["cases"]}
    filled = 0
    for tool, (fname, kind) in FILES.items():
        path = a.run / fname
        if not path.exists():
            print(f"{tool}: {path} is not there, rows stay 'not run'")
            continue
        hits = read_hits(path, kind)
        for c in doc["cases"]:
            targets = {t["acc"] for t in c["targets"]}
            mine = [h for h in hits
                    if (p := parse_query(h["query"])) and p[0] == c["id"] and h["target"] in targets]
            real = best([h for h in mine if not parse_query(h["query"])[2]])
            decoy = best([h for h in mine if parse_query(h["query"])[2]])
            verdict = classify(c, real, decoy)
            row = c["tools"][tool]
            row.update(verdict)
            row["score"] = None if real is None else f'{real["evalue"]:.2g} (E-value)'
            row["decoy_score"] = None if decoy is None else f'{decoy["evalue"]:.2g} (E-value)'
            row["call"] = None if real is None else f'{real["qstart"]}–{real["qend"]}'
            if real is not None:
                filled += 1
        print(f"{tool}: {len(hits)} hits read from {path}")
    a.cases.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
    print(f"wrote {a.cases}: {filled} case-and-tool rows now have a hit")


if __name__ == "__main__":
    main()
