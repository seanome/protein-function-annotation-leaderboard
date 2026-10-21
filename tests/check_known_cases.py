#!/usr/bin/env python3
"""Checks on the Known cases page and the data it is built from.

    .venv/bin/python tests/check_known_cases.py

It checks, without any network:
  the fourteen cases are the fixed list, in order, with the boxes the page claims;
  every tool and setting is the one the leaderboard tables call that tool's best;
  every aligned-residue block matches the stored sequences, residue for residue;
  every "found" mark has an IoU of at least 0.5 recomputed from the stored coordinates, and
    every hit beats its decoy;
  every decoy keeps the query's residue counts and dipeptide counts;
  every "Figure N" in the page's prose has a caption, and the numbers do not repeat one the
    rest of the site already uses;
  the truth table's first fourteen rows are these cases.
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

EXPECTED = [
    (1, 1, "P41958", "P10415"), (2, 1, "Q85295", "P10415"), (3, 1, "Q9Q8H3", "Q08722"),
    (4, 1, "Q57816", "P07437"), (5, 1, "Q9WZ57", "P60709"), (6, 1, "A0A0H3JV59", "P62314"),
    (7, 2, "Q01220", "P10415"), (8, 2, "P0DTS9", "O00206"), (9, 3, "H7C7N8", "Q08722"),
    (10, 3, "O52623", "Q13009"), (11, 3, "P74873", "Q07960"), (12, 3, "P24856", "A0A481T0N3"),
    (13, 4, "P03495", "P68431"), (14, 4, "P03254", "Q13363"),
]
IOU_FOUND = 0.5


def iou(a, b) -> float:
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    share = max(0, hi - lo + 1)
    union = (a[1] - a[0] + 1) + (b[1] - b[0] + 1) - share
    return share / union if union else 0.0


def main() -> int:
    bad = []
    doc = json.loads((ROOT / "data/known_cases/cases.json").read_text())
    seqs = json.loads((ROOT / "data/known_cases/sequences.json").read_text())
    cases = doc["cases"]

    if len(cases) != len(EXPECTED):
        bad.append(f"{len(cases)} cases, the list has {len(EXPECTED)}")
    for c, (cid, box, q, t) in zip(cases, EXPECTED):
        if (c["id"], c["box"], c["queries"][0]["acc"], c["targets"][0]["acc"]) != (cid, box, q, t):
            bad.append(f'case {c["id"]} is {c["queries"][0]["acc"]} against '
                       f'{c["targets"][0]["acc"]} in box {c["box"]}, the list says {q} against '
                       f'{t} in box {box}')

    # The settings are Table 1's best per tool, recomputed from the tables.
    try:
        import polars as pl  # noqa: F401
        import build_site as b
        res, comp, _src, _prov, _cov = b.load(ROOT / "data")
        for tool, setting in doc["tool_settings"].items():
            best = b.card_numbers(tool, res, comp).get("best_setting")
            if best != setting:
                bad.append(f"{tool}: the page uses {setting}, the tables' best is {best}")
    except Exception as e:  # the tables are large; say so rather than pass quietly
        bad.append(f"could not check the settings against the tables: {e}")

    for c in cases:
        a = c["alignment"]
        for side, acc in (("query", c["queries"][0]["acc"]), ("target", c["targets"][0]["acc"])):
            row = a[side].replace("-", "")
            start = a[f"{side}_start"]
            real = seqs[acc]["sequence"][start - 1:start - 1 + len(row)]
            if row != real:
                bad.append(f'case {c["id"]}: the {side} line does not match {acc} from residue '
                           f"{start}")
        line = "".join("|" if x == y and x != "-" else " " for x, y in zip(a["query"], a["target"]))
        if line != a["match"] or line.count("|") != a["identical"]:
            bad.append(f'case {c["id"]}: the match line or the identical count is not the one '
                       "the two rows give")

        for q in c["queries"]:
            s = seqs[q["acc"]]["sequence"]
            d = q["decoy"]["sequence"]
            if d == s:
                bad.append(f'case {c["id"]} {q["acc"]}: the decoy is the query itself')
            if Counter(s) != Counter(d):
                bad.append(f'case {c["id"]} {q["acc"]}: the decoy has different residue counts')
            if (Counter(x + y for x, y in zip(s, s[1:]))
                    != Counter(x + y for x, y in zip(d, d[1:]))):
                bad.append(f'case {c["id"]} {q["acc"]}: the decoy has different dipeptide counts')

        for tool, r in c["tools"].items():
            if r["result"] == "found":
                if r.get("iou") is None or r["iou"] < IOU_FOUND:
                    bad.append(f'case {c["id"]} {tool}: marked found with IoU {r.get("iou")}')
                call, truth = r.get("call"), c.get("query_region")
                if call and truth:
                    lo, hi = (int(x) for x in re.split(r"[–-]", call))
                    got = iou((lo, hi), tuple(truth))
                    if abs(got - r["iou"]) > 1e-6:
                        bad.append(f'case {c["id"]} {tool}: stored IoU {r["iou"]}, recomputed '
                                   f"{got:.3f}")
            if r["result"] != "not run" and r.get("score") and r.get("decoy_score"):
                real_e = float(re.split(r"[ (]", r["score"])[0])
                decoy_e = float(re.split(r"[ (]", r["decoy_score"])[0])
                if real_e >= decoy_e and r["result"] != "missed":
                    bad.append(f'case {c["id"]} {tool}: the decoy scores as well as the real '
                               "query, so the hit should not be counted")

    page = (ROOT / "docs/known-cases.html").read_text()
    used = set(re.findall(r"Figure (\d+)", page))
    captions = set(re.findall(r"<b>Figure (\d+)\.</b>", page))
    if used - captions:
        bad.append(f"the page names Figure {sorted(used - captions)} with no caption")
    others = set()
    for p in sorted((ROOT / "docs").rglob("*.html")):
        if p.name != "known-cases.html":
            others |= set(re.findall(r"<b>Figure (\d+)\.</b>", p.read_text()))
    if captions & others:
        bad.append(f"figure numbers {sorted(captions & others)} are used elsewhere on the site")

    truth = (ROOT / "data/truth_known_cases.tsv").read_text().splitlines()
    if len(truth) != 15:
        bad.append(f"the truth seed has {len(truth) - 1} rows, the page has 14 cases")
    for line, c in zip(truth[1:], cases):
        if line.split("\t")[2] != c["queries"][0]["acc"]:
            bad.append(f'the truth seed row {line.split(chr(9))[2]} is not case {c["id"]}')

    for b_ in bad:
        print("FAIL:", b_)
    print(f"known-cases checks: {len(cases)} cases, {len(bad)} problems")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
