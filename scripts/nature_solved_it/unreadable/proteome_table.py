"""One row per species: its UniProt reference proteome and its divergence time from human.

Model organisms (PLAN.md group 2):

    python3 -I scripts/nature_solved_it/unreadable/proteome_table.py --models \
        --out results/nature-solved-it/unreadable/model_organisms.csv

Atlas species (group 1), from the tolerance atlas CSV of Prompt 1, column 2 (species):

    python3 -I scripts/nature_solved_it/unreadable/proteome_table.py \
        --atlas <atlas.csv> --species-column species \
        --out results/nature-solved-it/unreadable/atlas_proteomes.csv

A species with no reference proteome keeps its row, with an empty `upid` and a note, so the
count of atlas species without one can be reported. Standard library only.
"""

import argparse
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import today  # noqa: E402

# PLAN.md group 2. Taxon ids are those of the reference proteome (yeast is strain S288C).
MODELS = (
    ("Homo sapiens", 9606),
    ("Mus musculus", 10090),
    ("Danio rerio", 7955),
    ("Drosophila melanogaster", 7227),
    ("Caenorhabditis elegans", 6239),
    ("Saccharomyces cerevisiae", 559292),
    ("Arabidopsis thaliana", 3702),
)

COLUMNS = [
    "group", "species", "upid", "proteome_taxid", "uniprot_organism", "uniprot_release",
    "readme_canonical_proteins", "uniprot_entries", "superregnum", "busco_lineage", "busco_lineage_all",
    "uniprot_busco", "divergence_my", "divergence_precomputed_my", "timetree_url",
    "timetree_read_utc", "fasta_url", "proteome_query_url", "read_utc", "note",
]


def row_for(group, species, prot_rows, release, url, readme, tt, note="", prefer_taxid=None):
    row = dict.fromkeys(COLUMNS, "")
    row.update(group=group, species=species, uniprot_release=release, proteome_query_url=url, read_utc=today.now_utc())
    if not prot_rows:
        row["note"] = "; ".join(x for x in (note, "no UniProt reference proteome") if x)
        return row
    if len(prot_rows) > 1:
        # Several reference proteomes under one species (strains, subspecies): take the one
        # filed under the species' own taxon id, else the most canonical proteins in the
        # README, ties broken by proteome id.
        prot_rows = sorted(prot_rows, key=lambda r: (
            str(r["Organism Id"]) != str(prefer_taxid),
            -int(readme.get(r["Proteome Id"], {}).get("#(1)", 0)),
            r["Proteome Id"],
        ))
        note = "; ".join(x for x in (note, f"chosen from {len(prot_rows)} reference proteomes") if x)
    p = prot_rows[0]
    upid, taxid = p["Proteome Id"], int(p["Organism Id"])
    rm = readme.get(upid, {})
    t = tt.lookup(taxid)
    division = {"eukaryota": "Eukaryota", "bacteria": "Bacteria", "archaea": "Archaea", "viruses": "Viruses"}.get(rm.get("SUPERREGNUM", ""), "")
    row.update(
        upid=upid,
        proteome_taxid=taxid,
        uniprot_organism=p["Organism"],
        readme_canonical_proteins=rm.get("#(1)", ""),
        uniprot_entries=p.get("Protein count", ""),
        superregnum=rm.get("SUPERREGNUM", ""),
        busco_lineage=today.busco_lineage(p["Taxonomic lineage"]),
        busco_lineage_all=today.BUSCO_ALL,
        uniprot_busco=p.get("BUSCO", ""),
        divergence_my="" if t["adjusted_age"] is None else t["adjusted_age"],
        divergence_precomputed_my="" if t.get("precomputed_age") is None else t["precomputed_age"],
        timetree_url=t["url"],
        timetree_read_utc=t["read_utc"],
        fasta_url=f"{today.FTP}/reference_proteomes/{division}/{upid}/{upid}_{taxid}.fasta.gz" if division else "",
        note="; ".join(x for x in (note, "" if rm else "not in the reference proteomes README", t.get("note", "")) if x),
    )
    return row


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--models", action="store_true")
    g.add_argument("--atlas", type=pathlib.Path, help="tolerance atlas CSV (Prompt 1)")
    ap.add_argument("--species-column", default="species")
    ap.add_argument("--out", required=True, type=pathlib.Path)
    args = ap.parse_args()

    out_dir = args.out.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    readme, readme_date = today.readme_table(out_dir / "reference_proteomes_README.txt.gz")
    tt = today.TimeTree(out_dir / "timetree_cache.json")
    rows = []
    if args.models:
        for species, taxid in MODELS:
            prot, release, url = today.reference_proteomes(f"organism_id:{taxid}")
            rows.append(row_for("model organism", species, prot, release, url, readme, tt))
    else:
        with open(args.atlas, newline="") as fh:
            names = [r[args.species_column].strip() for r in csv.DictReader(fh)]
        # One row per species; "none found" rows and blanks are not species.
        seen = []
        for n in names:
            if n and n.lower() != "none found" and n not in seen:
                seen.append(n)
        for species in seen:
            taxid, note = today.species_taxid(species)
            if taxid is None:
                rows.append(row_for("atlas species", species, [], "", "", readme, tt, note=note))
                continue
            prot, release, url = today.reference_proteomes(f"taxonomy_id:{taxid}")
            rows.append(row_for("atlas species", species, prot, release, url, readme, tt, note=note, prefer_taxid=taxid))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{args.out}: {len(rows)} rows, {sum(bool(r['upid']) for r in rows)} with a reference proteome; README Last-Modified {readme_date}")


if __name__ == "__main__":
    main()
