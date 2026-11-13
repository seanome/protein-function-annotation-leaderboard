"""Second-agent check: compare the recount (check_raw.tsv) with verdict.json and
unreadable_table.tsv, one row per proteome, written to check_rows.csv.

Usage: python3 -I check_compare.py
"""
import csv
import json

raw = {r["upid"]: r for r in csv.DictReader(open("check_raw.tsv"), delimiter="\t")}
ver = {p["upid"]: p for p in json.load(open("verdict.json"))["proteomes"]}
tab = {r["upid"]: r for r in csv.DictReader(open("run_2026-10-20/unreadable_table.tsv"),
                                            delimiter="\t")}
assert set(raw) == set(ver) == set(tab), (set(raw) ^ set(ver), set(raw) ^ set(tab))

cols = ["group", "species", "upid", "n_proteins_recount", "n_unlabelled_recount",
        "unlabelled_pct_recount", "n_unlabelled_verdict", "unlabelled_pct_verdict",
        "n_unlabelled_table", "unlabelled_pct_table", "n_uncharacterized_recount",
        "n_uncharacterized_verdict", "share_matches"]
n_match = 0
with open("check_rows.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(cols)
    for upid in sorted(ver, key=lambda u: (ver[u]["group"], ver[u]["species"])):
        r, v, t = raw[upid], ver[upid], tab[upid]
        pct = 100.0 * int(r["n_unlabelled"]) / int(r["n_proteins"])
        ok = (int(r["n_proteins"]) == v["n_proteins"] == int(t["n_proteins"])
              and int(r["n_unlabelled"]) == v["n_unlabelled"] == int(t["n_unlabelled"])
              and round(pct, 2) == v["unlabelled_pct"] == float(t["unlabelled_pct"])
              and int(r["n_uncharacterized_name"]) == v["n_uncharacterized_name"])
        n_match += ok
        w.writerow([v["group"], v["species"], upid, r["n_proteins"], r["n_unlabelled"],
                    "%.4f" % pct, v["n_unlabelled"], v["unlabelled_pct"],
                    t["n_unlabelled"], t["unlabelled_pct"], r["n_uncharacterized_name"],
                    v["n_uncharacterized_name"], "yes" if ok else "NO"])
print("proteomes matching on every count:", n_match, "of", len(ver))
