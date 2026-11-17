#!/bin/bash
# Fetch the structure Foldseek is given for each known case: the experimental structure the
# cited paper reports, or the AlphaFold Database model when there is none.
#
#   bash scripts/known_cases_structures.sh data/known_cases/cases.json <out dir>
#
# One file per query protein, named case<N>_<accession>.pdb (PDB entry) or .cif (AlphaFold).

set -uo pipefail
CASES=${1:?cases.json}
OUT=${2:?output directory}
# The AlphaFold Database answers 403 to a user agent that names curl or this project, and 200
# to a browser-style one, so the string below starts with Mozilla and then says who is asking.
UA=${KC_UA:-"Mozilla/5.0 (compatible; known-cases; +https://github.com/seanome/protein-function-annotation-leaderboard)"}
mkdir -p "$OUT"

python3 - "$CASES" <<'PY' | while IFS=$'\t' read -r name kind id chain; do
import json
import sys

doc = json.load(open(sys.argv[1]))
for c in doc["cases"]:
    for q in c["queries"]:
        st = q.get("structure") or c.get("structure")
        if st:
            print(f'case{c["id"]}_{q["acc"]}\tpdb\t{st["pdb"]}\t{st["chain"]}')
        else:
            print(f'case{c["id"]}_{q["acc"]}\tafdb\t{q["acc"]}\t-')
PY
  if [ "$kind" = pdb ]; then
    url="https://files.rcsb.org/download/${id}.pdb"
    out="$OUT/${name}.pdb"
  else
    # Ask the AlphaFold Database which file it has: the model version moves (v4, v6, ...), so
    # a version written into the URL by hand stops working.
    url=$(curl -fsS -A "$UA" "https://alphafold.ebi.ac.uk/api/prediction/${id}" 2>/dev/null |
          python3 -c 'import json,sys
try:
    d = json.load(sys.stdin)
except Exception:
    d = []
print(d[0]["cifUrl"] if d else "")')
    out="$OUT/${name}.cif"
    if [ -z "$url" ]; then
      echo "MISSING $name: the AlphaFold Database has no model for $id; the card stays \"not run\"" >&2
      continue
    fi
  fi
  if [ -s "$out" ]; then echo "have $out"; continue; fi
  if curl -fsS -A "$UA" -o "$out" "$url"; then
    echo "fetched $out from $url (chain $chain)"
  else
    echo "MISSING $name: $url did not download; the card stays \"not run\"" >&2
    rm -f "$out"
  fi
done
