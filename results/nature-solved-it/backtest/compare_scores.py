"""Compare the checking agent's blind re-score with blind_scores.csv."""
import csv
from pathlib import Path
HERE = Path(__file__).resolve().parent
MINE = HERE / "rescore.csv"
mine = {r["blind_id"]: r for r in csv.DictReader(MINE.open())}
theirs = {r["blind_id"]: r for r in csv.DictReader((HERE / "blind_scores.csv").open())}
key = {r["blind_id"]: r for r in csv.DictReader((HERE / "key.csv").open())}
pts = ["T%d" % i for i in range(1, 7)]
agree_pt = {p: 0 for p in pts}
total_match = 0
print("id    entry  programme                     mine theirs  differing points")
for bid in mine:
    m, t = mine[bid], theirs[bid]
    diff = [p for p in pts if int(m[p]) != int(t[p])]
    for p in pts:
        if int(m[p]) == int(t[p]):
            agree_pt[p] += 1
    mt, tt = sum(int(m[p]) for p in pts), int(t["T"])
    if mt == tt:
        total_match += 1
    print("%-5s %-6s %-30s %-4d %-6d %s" % (bid, key[bid]["entry_id"],
          key[bid]["drug_or_programme"][:30], mt, tt,
          ", ".join("%s mine %s theirs %s" % (p, m[p], t[p]) for p in diff) or "none"))
print("\ntotals matching exactly: %d / %d" % (total_match, len(mine)))
print("per-point agreement:", {p: "%d/10" % agree_pt[p] for p in pts},
      "| all points: %d/60" % sum(agree_pt.values()))
