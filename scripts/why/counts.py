"""Steps 1-3 of the "Why this project" numbers: Swiss-Prot by organism, PDB by organism,
species counts by kingdom. Every number goes into results/why/numbers.parquet with its
unit and the sources.json key of the file it was computed from.

    python3 -I scripts/why/counts.py --raw ~/data/why-numbers/2026_03 metazoan-names out.txt
    python3 -I scripts/why/counts.py --raw ~/data/why-numbers/2026_03 counts
"""

import argparse
import ast
import collections
import gzip
import json
import pathlib
import sys

import polars as pl

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from numbers_table import Numbers  # noqa: E402
from taxa import COL_GROUP_TAXA, Taxonomy  # noqa: E402

OUT = pathlib.Path("results/why")
# The deck's numbers as stated on slide 6 of seanome/presentations 2026-09-30-kmerseek-status.
# deck_numbers() recomputes them from the file the deck names.
DECK_STATED = {"entries": 62253, "organisms": 5918, "top5_pct": 26.4}


def deck_numbers(path):
    """The deck's method: cellular life, PDB structure counts summed by organism name."""
    d = pl.read_parquet(path).filter(pl.col("superdomain") == "Cellular Life", pl.col("pdb_structures_count") > 0)
    c = d.group_by("organism").agg(pl.col("pdb_structures_count").sum())
    total = c["pdb_structures_count"].sum()
    top = c.sort("pdb_structures_count", descending=True).head(5)["pdb_structures_count"].sum()
    return {"entries": total, "organisms": c.height, "top5_pct": round(100 * top / total, 2)}


def read_sprot_dat(path):
    """One row per Swiss-Prot entry: primary accession, OX taxon, number of PDB entries."""
    rows, acc, taxid, n_pdb = [], None, None, 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            tag = line[:2]
            if tag == "AC" and acc is None:
                acc = line[5:].split(";")[0].strip()
            elif tag == "OX":
                taxid = int(line.split("NCBI_TaxID=")[1].split(";")[0].split()[0])
            elif tag == "DR" and line[5:9] == "PDB;":
                n_pdb += 1
            elif tag == "//":
                rows.append((acc, taxid, n_pdb))
                acc, taxid, n_pdb = None, None, 0
    return pl.DataFrame(rows, schema=["accession", "taxid", "n_pdb"], orient="row")


def annotate(df, tax):
    taxids = df["taxid"].unique().to_list()
    info = pl.DataFrame(
        {
            "taxid": taxids,
            "species_taxid": [tax.species(t) for t in taxids],
            "superkingdom": [tax.superkingdom(t) for t in taxids],
            "col_group": [tax.col_group(t) for t in taxids],
            "metazoan": [tax.is_metazoan(t) for t in taxids],
            "organism": [tax.name.get(tax.resolve(t)) for t in taxids],
        },
        schema_overrides={"species_taxid": pl.Int64},
    )
    df = df.join(info, on="taxid", how="left")
    return df.with_columns(
        species_name=pl.col("species_taxid").map_elements(lambda t: tax.name.get(t), return_dtype=pl.String)
    )


def top_share(counts, key, value, n=5):
    c = counts.group_by(key).agg(pl.col(value).sum()).sort([value, key], descending=[True, False])
    total = c[value].sum()
    top = c.head(n)
    return c, total, top, 100 * top[value].sum() / total


def load_worms(path):
    """Species name -> habitat flags from the first exact, accepted WoRMS match."""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            row = json.loads(line)
            exact = [m for m in row["matches"] if m.get("match_type") == "exact"]
            exact.sort(key=lambda m: m.get("status") != "accepted")
            if exact:
                m = exact[0]
                out[row["query"]] = {k: m.get(k) for k in ("AphiaID", "status", "isMarine", "isBrackish", "isFreshwater", "isTerrestrial")}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=pathlib.Path)
    sub = ap.add_subparsers(dest="cmd", required=True)
    n = sub.add_parser("metazoan-names")
    n.add_argument("out")
    sub.add_parser("counts")
    args = ap.parse_args()
    raw = args.raw

    tax = Taxonomy(raw / "taxdump.tar.gz")
    for group, taxa in COL_GROUP_TAXA.items():
        print(group, [tax.name.get(t) for t in taxa])
    sp = annotate(read_sprot_dat(raw / "uniprot_sprot.dat.gz"), tax)

    if args.cmd == "metazoan-names":
        names = sorted(sp.filter(pl.col("metazoan"), pl.col("species_name").is_not_null())["species_name"].unique())
        pathlib.Path(args.out).write_text("\n".join(names) + "\n")
        print(len(names), "metazoan species names")
        return

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    N = Numbers()
    src_sp = "swissprot_dat"
    rel = "UniProt 2026_03 (02-Sep-2026)"

    # ---- Step 1: Swiss-Prot by organism -------------------------------------------------
    N.add("sp_entries", 1, "Swiss-Prot reviewed entries", sp.height, "entries", src_sp, rel)
    sp = sp.with_columns(one=pl.lit(1))
    for key, label in (("taxid", "organism (UniProt OX taxon)"), ("species_taxid", "species (NCBI rank species)")):
        c, total, top, share = top_share(sp.filter(pl.col(key).is_not_null()), key, "one")
        unit = "organisms" if key == "taxid" else "species"
        N.add(f"sp_{unit}", 1, f"{unit} with at least one Swiss-Prot entry", c.height, unit, src_sp, rel, note=f"grouped by {label}")
        N.add(f"sp_top5_share_by_{unit}", 1, f"share of Swiss-Prot entries held by the five largest {unit}", round(share, 2), "% of entries", src_sp, rel,
              note="top five: " + "; ".join(f"{tax.name.get(t)} {v:,}" for t, v in top.iter_rows()))
        c.rename({"one": "entries"}).with_columns(
            name=pl.col(key).map_elements(lambda t: tax.name.get(t), return_dtype=pl.String)
        ).write_parquet(OUT / "tables" / f"swissprot_entries_per_{key}.parquet")
    no_species = sp.filter(pl.col("species_taxid").is_null()).height
    N.add("sp_entries_no_species", 1, "Swiss-Prot entries whose taxon has no species above it (left out of species counts)", no_species, "entries", src_sp, rel)

    # Marine animals: metazoan species with a WoRMS isMarine = 1 record.
    worms = load_worms(raw / "worms_matches.jsonl.gz")
    met = sp.filter(pl.col("metazoan"))
    met_species = met["species_name"].drop_nulls().unique().to_list()
    matched = [s for s in met_species if s in worms]
    marine = {s for s in matched if worms[s]["isMarine"] == 1}
    marine_only = {s for s in marine if not worms[s]["isFreshwater"] and not worms[s]["isTerrestrial"]}
    src_w = "worms_matches"
    N.add("sp_metazoa_entries", 1, "Swiss-Prot entries from animals (NCBI Metazoa)", met.height, "entries", src_sp, rel)
    N.add("sp_metazoa_species", 1, "animal species with at least one Swiss-Prot entry", len(met_species), "species", src_sp, rel)
    N.add("worms_matched_species", 1, "of those, species names with an exact WoRMS match", len(matched), "species", src_w, "WoRMS REST API")
    N.add("worms_marine_species", 1, "of those, species WoRMS flags as marine (isMarine = 1)", len(marine), "species", src_w, "WoRMS REST API")
    n_marine = met.filter(pl.col("species_name").is_in(list(marine))).height
    n_marine_only = met.filter(pl.col("species_name").is_in(list(marine_only))).height
    N.add("sp_marine_animal_entries", 1, "Swiss-Prot entries from animal species WoRMS flags as marine", n_marine, "entries", src_w, "WoRMS REST API + UniProt 2026_03")
    N.add("sp_marine_animal_share", 1, "share of all Swiss-Prot entries from animal species WoRMS flags as marine", round(100 * n_marine / sp.height, 2), "% of entries", src_w, "WoRMS REST API + UniProt 2026_03",
          note="marine = isMarine 1; a species may also live in brackish or fresh water")
    N.add("sp_marine_only_animal_share", 1, "same, counting only species WoRMS flags as marine and not freshwater or terrestrial", round(100 * n_marine_only / sp.height, 2), "% of entries", src_w, "WoRMS REST API + UniProt 2026_03")
    pl.DataFrame([{"species_name": s, **worms.get(s, {})} for s in met_species]).write_parquet(OUT / "tables" / "metazoan_species_worms.parquet")

    # ---- Step 2: experimental structures ------------------------------------------------
    sp_pdb = sp.filter(pl.col("n_pdb") > 0)
    N.add("sp_pdb_entries", 2, "Swiss-Prot entries citing at least one PDB entry", sp_pdb.height, "entries", src_sp, rel)
    for key, unit in (("taxid", "organisms"), ("species_taxid", "species")):
        c, total, top, share = top_share(sp_pdb.filter(pl.col(key).is_not_null()), key, "one")
        N.add(f"sp_pdb_{unit}", 2, f"{unit} with a Swiss-Prot entry citing a PDB entry", c.height, unit, src_sp, rel)
        N.add(f"sp_pdb_top5_share_by_{unit}", 2, f"share of those entries held by the five largest {unit}", round(share, 2), "% of entries", src_sp, rel,
              note="top five: " + "; ".join(f"{tax.name.get(t)} {v:,}" for t, v in top.iter_rows()))

    # All UniProtKB (Swiss-Prot + TrEMBL), cellular life: the deck's scope.
    kb = pl.read_csv(raw / "uniprotkb_with_pdb.tsv.gz", separator="\t", infer_schema_length=0).rename(
        {"Entry": "accession", "Reviewed": "reviewed", "Organism (ID)": "taxid", "Organism": "organism_name", "PDB": "pdb"}
    ).with_columns(pl.col("taxid").cast(pl.Int64), one=pl.lit(1))
    kb = annotate(kb.select("accession", "reviewed", "taxid", "organism_name", "one"), tax)
    src_kb = "uniprotkb_with_pdb"
    kb_rel = json.loads(pathlib.Path("results/why/sources.json").read_text())[src_kb]["release"]
    N.add("kb_pdb_entries_all", 2, "UniProtKB entries (Swiss-Prot and TrEMBL) citing a PDB entry, all organisms", kb.height, "entries", src_kb, kb_rel)
    cell = kb.filter(pl.col("superkingdom").is_in(["Bacteria", "Archaea", "Eukaryota"]))
    N.add("kb_pdb_entries", 2, "same, cellular life only (the deck's scope)", cell.height, "entries", src_kb, kb_rel)
    for key, unit, how in (
        ("organism_name", "organism names", "the deck's method: summed by organism name"),
        ("taxid", "organisms", "summed by UniProt OX taxon"),
        ("species_taxid", "species", "summed by NCBI species"),
    ):
        c, total, top, share = top_share(cell.filter(pl.col(key).is_not_null()), key, "one")
        names = [r[0] if key == "organism_name" else tax.name.get(r[0]) for r in top.iter_rows()]
        N.add(f"kb_pdb_{key}", 2, f"{unit} with a UniProtKB entry citing a PDB entry, cellular life", c.height, unit, src_kb, kb_rel, note=how)
        N.add(f"kb_pdb_top5_share_by_{key}", 2, f"share held by the five largest {unit}, cellular life", round(share, 2), "% of entries", src_kb, kb_rel,
              note=how + "; top five: " + "; ".join(f"{nm} {v:,}" for nm, (_, v) in zip(names, top.iter_rows())))
        c.rename({"one": "entries"}).write_parquet(OUT / "tables" / f"uniprotkb_pdb_entries_per_{key}.parquet")
    DECK = deck_numbers(raw / "wealth_inequality_uniprot_combined_sprot_trembl_taxonomy_summary.parquet")
    for what, v in DECK.items():
        N.add(f"deck_{what}_recomputed", 2, f"deck value recomputed from its source table: {what} (slide states {DECK_STATED[what]:,})", v,
              "% of entries" if what == "top5_pct" else ("entries" if what == "entries" else "organism names"), "deck_parquet",
              "UniProt as parsed 2025-03-11 in seanome/2025-biodata-wealth-inequality")
    deck_cmp = [("entries", cell.height, DECK["entries"])]
    c_name = cell.group_by("organism_name").agg(pl.len())
    deck_cmp.append(("organisms", c_name.height, DECK["organisms"]))
    deck_cmp.append(("top5_pct", round(top_share(cell, "organism_name", "one")[3], 2), DECK["top5_pct"]))
    for what, now, then in deck_cmp:
        N.add(f"deck_{what}_now_minus_then", 2, f"current minus deck value: {what} (deck {then:,})", round(now - then, 2),
              "% points" if what == "top5_pct" else ("entries" if what == "entries" else "organism names"), src_kb, kb_rel,
              note="deck value from seanome/2025-biodata-wealth-inequality uniprot_combined_sprot_trembl_taxonomy_summary.parquet (UniProt as parsed 2025-03-11)")

    # SIFTS, the PDB side: UniProt accessions with at least one PDB chain mapped.
    sifts = pl.read_csv(raw / "sifts_uniprot_pdb.tsv.gz", separator="\t", comment_prefix="#", infer_schema_length=0)
    sifts_acc = set(sifts[sifts.columns[0]].to_list())
    kb_acc = set(kb["accession"].to_list())
    src_s = "sifts_uniprot_pdb"
    N.add("sifts_accessions", 2, "UniProt accessions with a PDB chain mapped by SIFTS", len(sifts_acc), "accessions", src_s, "SIFTS flat file")
    N.add("sifts_and_uniprot", 2, "of those, also in UniProt's own PDB cross-references", len(sifts_acc & kb_acc), "accessions", src_s, "SIFTS + UniProt 2026_03")

    # ---- Step 3: species by kingdom -----------------------------------------------------
    col = json.loads((raw / "col_accepted_species.json").read_text())
    col_rel = col["dataset"]["alias"] + " " + col["dataset"]["version"]
    sp_species = sp.filter(pl.col("species_taxid").is_not_null()).unique("species_taxid")
    per_group = sp_species.group_by("col_group").agg(pl.len().alias("n"))
    per_group = dict(per_group.iter_rows())
    rows3 = []
    for c in col["counts"]:
        g = c["col_group"]
        if g == "all":
            continue
        described = c["accepted_species"]
        with_sp = per_group.get(g, 0)
        N.add(f"col_species_{g}", 3, f"described species, {g} (Catalogue of Life accepted species)", described, "species", "col_accepted_species", col_rel)
        N.add(f"sp_species_{g}", 3, f"species with at least one Swiss-Prot entry, {g}", with_sp, "species", src_sp, rel, note="NCBI taxa mapped to COL groups, see taxa.py")
        N.add(f"sp_species_share_{g}", 3, f"species with a Swiss-Prot entry per described species, {g}", round(100 * with_sp / described, 3), "% of described species", src_sp, rel)
        rows3.append({"col_group": g, "described_species": described, "species_with_swissprot_entry": with_sp})
    # Larsen et al. 2017 (https://doi.org/10.1086/693564), Table 1 scenario 1: their
    # ESTIMATED total species, as the deck uses them, read from the notebook that tables them.
    nb = json.loads((raw / "wealth_inequality_03_make_predicted_n_proteins.ipynb").read_text())
    cell = next("".join(c["source"]) for c in nb["cells"] if "larsen2017 = pd.Series(" in "".join(c["source"]))
    larsen = ast.literal_eval(cell[cell.index("{", cell.index("pd.Series(")) : cell.index("}") + 1])
    larsen_to_col = {"Animals": "Animalia", "Plants": "Plantae", "Fungi": "Fungi", "Protists": "Chromista + Protozoa", "Bacteria": "Bacteria + Archaea"}
    for group, est in larsen.items():
        N.add(f"larsen2017_estimated_species_{group}", 3, f"ESTIMATED total species, {group} (Larsen et al. 2017, Table 1, scenario 1)", est, "species (estimate)",
              "larsen2017_notebook", "Larsen et al. 2017, Q. Rev. Biol. 92:229-265, doi 10.1086/693564",
              note=f"an estimate, not a count; compare with COL group(s) {larsen_to_col[group]}")
        rows3.append({"col_group": f"Larsen 2017: {group}", "larsen2017_estimated_species": est})
    pl.DataFrame(rows3, infer_schema_length=None).write_parquet(OUT / "tables" / "species_by_kingdom.parquet")

    N.write(OUT / "numbers.parquet", replace_steps={1, 2, 3})
    print(N.frame().select("id", "value", "unit").to_pandas().to_string())


if __name__ == "__main__":
    main()
