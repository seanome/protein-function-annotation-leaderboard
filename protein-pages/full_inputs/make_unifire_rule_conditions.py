"""Write unifire_rule_conditions.tsv: which signatures trigger each rule UniFIRE fired on this proteome.

    python make_unifire_rule_conditions.py unifire_predictions.tsv.gz | gzip > unifire_rule_conditions.tsv.gz

UniFIRE's output says which rule fired, not where. A rule fires on InterPro, Pfam, PROSITE and other
signatures, and InterProScan gives those signatures' positions, so the rule's conditions say where on
the protein it points. Conditions come from UniProt's REST API (rest.uniprot.org/arba and /unirule,
read page by page). One row per signature condition: rule, the condition set it belongs to (a rule
fires when every condition of one set holds), and the accepted signatures, comma-separated. Negative
conditions and taxon, length and other non-signature conditions are left out.
"""
import csv
import gzip
import json
import re
import sys
import time
import urllib.request

SIG = re.compile(r"^(InterPro|Pfam|PROSITE|PRINTS|SMART|PANTHER|NCBIfam|TIGRFAMs?|HAMAP|CDD|Gene3D|CATH|SUPFAM|SUPERFAMILY|PIRSF|SFLD|FunFam)\b", re.I)


def pages(url):
    while url:
        for attempt in range(6):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "protein-pages (olga@seanome.org)"}), timeout=120) as r:
                    link = r.headers.get("Link") or ""
                    data = json.load(r)
                break
            except Exception as e:
                print("retry", e, file=sys.stderr)
                time.sleep(5 * (attempt + 1))
        else:
            sys.exit("gave up on " + url)
        yield from data["results"]
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else None


def condition_sets(rule):
    out = list(rule.get("mainRule", {}).get("conditionSets", []))
    for case in rule.get("otherRules", []) + rule.get("caseRules", []) + rule.get("samFeatureSets", []) + rule.get("positionFeatureSets", []):
        out += case.get("conditionSets", [])
    return out


def main(predictions):
    want = {r["rule"] for r in csv.DictReader(gzip.open(predictions, "rt"), delimiter="\t")}
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["rule", "set", "type", "signatures"])
    seen, types = set(), {}
    for base in ("arba", "unirule"):
        for rule in pages("https://rest.uniprot.org/%s/search?query=*&size=500&format=json" % base):
            rid = rule["uniRuleId"]
            if rid not in want or rid in seen:
                continue
            seen.add(rid)
            for i, cs in enumerate(condition_sets(rule)):
                for c in cs.get("conditions", []):
                    types[c.get("type")] = types.get(c.get("type"), 0) + 1
                    if c.get("isNegative") or not SIG.match(c.get("type") or ""):
                        continue
                    w.writerow([rid, i, c["type"], ",".join(v["value"] for v in c.get("conditionValues", []))])
    print("%d of %d rules found; condition types: %s" % (len(seen), len(want), sorted(types.items(), key=lambda x: -x[1])), file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1])
