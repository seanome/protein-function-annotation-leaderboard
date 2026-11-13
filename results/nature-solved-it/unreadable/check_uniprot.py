"""Second-agent check: five proteome ids against rest.uniprot.org today.

Compares UniProt's geneCount (proteins in the one-protein-per-gene FASTA) and taxon with
n_proteins, proteome_taxid and species in verdict.json.

Usage: python3 -I check_uniprot.py
"""
import json
import urllib.request

PICK = ["UP000005640", "UP000694906", "UP000186922", "UP001619887", "UP001710656"]
verdict = {p["upid"]: p for p in json.load(open("verdict.json"))["proteomes"]}

print("upid\tspecies_verdict\tspecies_uniprot\ttaxid_verdict\ttaxid_uniprot\t"
      "n_proteins_verdict\tgeneCount_uniprot\tproteinCount_uniprot\tproteomeType\trelease")
for upid in PICK:
    url = "https://rest.uniprot.org/proteomes/%s?format=json" % upid
    with urllib.request.urlopen(url, timeout=60) as r:
        release = r.headers.get("x-uniprot-release")
        d = json.load(r)
    v = verdict[upid]
    print("\t".join(str(x) for x in [
        upid, v["species"], d["taxonomy"]["scientificName"], v["proteome_taxid"],
        d["taxonomy"]["taxonId"], v["n_proteins"], d.get("geneCount"),
        d.get("proteinCount"), d.get("proteomeType"), release]))
