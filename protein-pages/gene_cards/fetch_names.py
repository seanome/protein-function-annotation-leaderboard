"""Fetch the names the gene cards need and the annotation tables leave out:
the UniProt protein name, gene and organism of every Swiss-Prot match and every FANTASIA
reference protein, and the name of every GO term. Writes two files into full_inputs/:

  uniprot_names.tsv.gz  accession, protein_name, gene, organism
  go_names.tsv.gz       go_id, name, namespace

Run: python gene_cards/fetch_names.py full_inputs
"""
import csv
import gzip
import sys
import time
import urllib.request
from pathlib import Path

UNIPROT = "https://rest.uniprot.org/uniprotkb/accessions?accessions={}&fields=accession,protein_name,gene_primary,organism_name&format=tsv&size=500"
GO_OBO = "https://purl.obolibrary.org/obo/go/go-basic.obo"


def get(url, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "protein-pages/fetch_names (olgabot)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode()
        except Exception as e:  # network errors; UniProt sometimes answers 5xx
            if i == tries - 1:
                raise
            print("retry", i + 1, e, file=sys.stderr)
            time.sleep(5 * (i + 1))


def main(folder):
    folder = Path(folder)
    rows = list(csv.DictReader(gzip.open(folder / "botryllus_protein_annotations_v2.tsv.gz", "rt"), delimiter="\t"))
    accs = sorted({r["sp_acc"] for r in rows if r["sp_acc"]} | {r["fantasia_uniprot"] for r in rows if r["fantasia_uniprot"]})
    go = sorted({t for r in rows for t in r["go_terms"].split(";") if t}
                | {t.split("|")[0] for r in rows for t in r["fantasia_go_terms"].split(";") if t})

    out, seen = [], set()
    for i in range(0, len(accs), 400):
        text = get(UNIPROT.format(",".join(accs[i:i + 400])))
        for line in text.splitlines()[1:]:
            f = line.split("\t")
            if f[0] not in seen:
                seen.add(f[0])
                out.append(f[:4])
        print(len(seen), "of", len(accs), file=sys.stderr)
    with gzip.open(folder / "uniprot_names.tsv.gz", "wt") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "protein_name", "gene", "organism"])
        w.writerows(sorted(out))
    missing = [a for a in accs if a not in seen]
    print("UniProt names:", len(seen), "found,", len(missing), "not found", missing[:10], file=sys.stderr)

    names, cur = {}, {}
    for line in get(GO_OBO).splitlines() + ["[Term]"]:
        if line.startswith("[Term]") or line.startswith("[Typedef]"):
            if "id" in cur:
                names[cur["id"]] = (cur.get("name", ""), cur.get("namespace", ""))
                for a in cur.get("alt", []):
                    names.setdefault(a, names[cur["id"]])
            cur = {}
        elif line.startswith("id: "):
            cur["id"] = line[4:]
        elif line.startswith("name: "):
            cur["name"] = line[6:]
        elif line.startswith("namespace: "):
            cur["namespace"] = line[11:]
        elif line.startswith("alt_id: "):
            cur.setdefault("alt", []).append(line[8:])
    with gzip.open(folder / "go_names.tsv.gz", "wt") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["go_id", "name", "namespace"])
        for g in go:
            if g in names:
                w.writerow([g, *names[g]])
    print("GO names:", sum(g in names for g in go), "of", len(go), file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1])
