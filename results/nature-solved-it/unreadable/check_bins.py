"""Second-agent check: each control's bin against its TimeTree time (cache and live today),
and the per-bin counts of controls against atlas species.

Usage: python3 -I check_bins.py   (from this folder; needs network for the live read)
"""
import collections
import csv
import json
import urllib.request

tt = json.load(open("timetree_cache.json"))
controls = list(csv.DictReader(open("controls.csv")))
atlas = list(csv.DictReader(open("atlas_proteomes.csv")))
verdict = {p["upid"]: p for p in json.load(open("verdict.json"))["proteomes"]}

print("species\ttaxid\tcsv_my\tcache_my\tlive_my\tbin_csv\tbin_live\tverdict_my")
for c in controls:
    tax = c["proteome_taxid"]
    live = None
    try:
        with urllib.request.urlopen(
                "https://timetree.temple.edu/api/pairwise/9606/" + tax, timeout=60) as r:
            # TimeTree answers CSV text, one header line and one row
            rows = list(csv.DictReader(r.read().decode().splitlines()))
            live = float(rows[0]["adjusted_age"])
    except Exception as e:  # report, do not hide
        live = "error: %s" % e
    csv_my = float(c["divergence_my"])
    cache_my = tt[tax]["adjusted_age"]
    bin_live = int(float(live) // 100) if isinstance(live, (int, float)) else "NA"
    print("\t".join(str(x) for x in [c["species"], tax, csv_my, cache_my, live,
                                     int(csv_my // 100), bin_live,
                                     verdict[c["upid"]]["divergence_my"]]))

print()
a_bins = collections.Counter(int(float(verdict[r["upid"]]["divergence_my"]) // 100)
                             for r in atlas if r["upid"] in verdict
                             and verdict[r["upid"]]["group"] == "atlas species")
c_bins = collections.Counter(int(float(c["divergence_my"]) // 100) for c in controls)
print("bin\tatlas\tcontrols")
for b in sorted(set(a_bins) | set(c_bins)):
    print("%d\t%d\t%d" % (b, a_bins[b], c_bins[b]))
