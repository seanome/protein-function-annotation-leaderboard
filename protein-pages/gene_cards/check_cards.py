"""Check gene cards against the evidence they were written from, and compare two models' name
verdicts on the same proteins.

  python gene_cards/check_cards.py gene_cards/out claude-haiku-5-5 claude-opus-5-5

Checks, per card:
  numbers   every percentage and E-value written in the summary or verdict is on the evidence
            sheet, allowing for rounding (a number not on the sheet was made up or mis-copied)
  dark      a protein whose only evidence is FANTASIA is not given "high" or "medium" confidence
  verdict   a conflict gets a verdict other than not_applicable, and a non-conflict gets
            not_applicable
Writes out/check.tsv (one row per protein) and prints the counts.
"""
import collections
import csv
import json
import re
import sys
from pathlib import Path

NUM = re.compile(r"(\d+(?:\.\d+)?)\s*%|(\d(?:\.\d+)?e-\d+)", re.I)


def numbers(text):
    out = set()
    for pct, e in NUM.findall(text):
        out.add(("pct", float(pct)) if pct else ("e", float(e)))
    return out


PROMPT_PCT = {20, 25, 40, 50}  # identity and coverage levels prompt.md itself names


def on_sheet(n, sheet):
    """A number counts as on the sheet when it is a sheet value rounded (27% for 27.2%; 4e-163
    for 4.29e-163), 100 minus a sheet percentage (the unmatched 87% of a 13% match), or a level
    prompt.md names."""
    kind, v = n
    if kind == "pct":
        vals = [float(x) for x in re.findall(r"(\d+(?:\.\d+)?)\s*%", sheet)]
        return v in PROMPT_PCT or any(abs(v - x) <= 0.51 or abs(v - (100 - x)) <= 0.51 for x in vals)
    m, ex = ("%e" % v).split("e")
    for x in re.findall(r"\d(?:\.\d+)?e-\d+", sheet, re.I):
        xm, xe = ("%e" % float(x)).split("e")
        if int(ex) == int(xe) and abs(float(m) - float(xm)) < 1:  # rounded or cut to one figure
            return True
    return False


def main(out, *models):
    out = Path(out)
    ev = {e["protein_id"]: e for e in map(json.loads, open(out / "evidence.jsonl"))}
    cards = {m: {c["protein_id"]: c for c in map(json.loads, open(out / ("cards.%s.jsonl" % m)))} for m in models}
    rows, counts = [], collections.Counter()
    for pid, e in ev.items():
        conflict = e["kind"].startswith("conflict")
        row = {"protein_id": pid, "kind": e["kind"]}
        for m in models:
            c = cards[m].get(pid)
            if not c:
                if conflict or m == models[0]:
                    counts[m, "missing"] += 1
                continue
            bad = sorted("%g%s" % (v, "%" if k == "pct" else "") for k, v in numbers(c["summary"] + " " + c["name_verdict"]["reason"]) if not on_sheet((k, v), e["evidence"]))
            dark_bad = e["kind"] == "fantasia_only" and c["confidence"] in ("high", "medium")
            v = c["name_verdict"]["choice"]
            verdict_bad = (v == "not_applicable") == conflict
            counts[m, "cards"] += 1
            counts[m, "numbers_not_on_sheet"] += bool(bad)
            counts[m, "dark_overconfident"] += dark_bad
            counts[m, "verdict_wrong_kind"] += verdict_bad
            row.update({m + ":name": c["card_name"], m + ":confidence": c["confidence"], m + ":verdict": v,
                        m + ":recommended": c["name_verdict"]["recommended_name"],
                        m + ":numbers_not_on_sheet": " ".join(bad), m + ":dark_overconfident": int(dark_bad),
                        m + ":verdict_wrong_kind": int(verdict_bad)})
        rows.append(row)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out / "check.tsv", "w") as fh:
        w = csv.DictWriter(fh, keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for k in sorted(counts):
        print("%s\t%s\t%d" % (k[0], k[1], counts[k]))
    if len(models) == 2:
        a, b = models
        both = [r for r in rows if r["kind"].startswith("conflict") and a + ":verdict" in r and b + ":verdict" in r]
        agree = sum(r[a + ":verdict"] == r[b + ":verdict"] for r in both)
        print("verdicts on %d conflicts: %s and %s agree on %d" % (len(both), a, b, agree))
        for (x, y), n in sorted(collections.Counter((r[a + ":verdict"], r[b + ":verdict"]) for r in both).items()):
            print("  %s -> %s: %d" % (x, y, n))


if __name__ == "__main__":
    main(*sys.argv[1:])
