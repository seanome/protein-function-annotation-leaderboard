#!/usr/bin/env python3
"""The four card lookups for every pre-built species, and the tables the browser needs to do
the same lookups for any other species.

    python3 -I scripts/organism/lookups.py --raw ~/data/organism-page relatives
        writes data/organism/relatives.parquet (nearest reviewed relative of each pre-built
        species); run `fetch.py --step pairs` after it
    python3 -I scripts/organism/lookups.py --raw ~/data/organism-page tables
        writes data/organism/lookups.jsonl (one record per pre-built species) and the browser
        tables under data/organism/web/

The nearest-relative rule here and in docs/assets/organism.js is the same, and
tests/check_organism.py compares the two on live pages:
  1. Candidates are species with at least 100 reviewed Swiss-Prot entries, other than the
     species itself.
  2. Walk up the species' NCBI lineage; the first ancestor that has a candidate under it is
     the shared ancestor.
  3. Under it, the candidate with most reviewed entries wins (ties: lower taxon id).
"""

import argparse
import json
import pathlib
import re
import sys

import polars as pl

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import fetch  # noqa: E402

OUT = pathlib.Path("data/organism")
WEB = OUT / "web"
RANKS = ("order", "class")


def load_tax(raw):
    parent, n = fetch.lineage_table(raw)
    rank = dict(zip(n["taxid"].to_list(), n["rank"].to_list()))
    name = dict(zip(n["taxid"].to_list(), n["name"].to_list()))
    common = {t: c for t, c in zip(n["taxid"].to_list(), n["common_name"].to_list()) if c}
    return parent, rank, name, common


def top2_under(sp, parent):
    """node -> the two reviewed species under it with most entries, as [[taxid, entries], ...]."""
    under = {}
    for t, e in sp.iter_rows():
        for a in fetch.ancestors(t, parent)[1:]:
            lst = under.setdefault(a, [])
            lst.append((t, e))
    out = {}
    for a, lst in under.items():
        lst.sort(key=lambda x: (-x[1], x[0]))
        out[a] = [list(x) for x in lst[:2]]
    return out


def nearest(t, lineage, under):
    for a in lineage[1:]:
        for cand, e in under.get(a, []):
            if cand != t:
                return a, cand, e
    return None, None, None


def reviewed(min_n=fetch.MIN_REVIEWED):
    return pl.read_parquet(OUT / "swissprot_per_species.parquet").filter(pl.col("entries") >= min_n) \
        .select("species_taxid", "entries")


def cmd_relatives(raw):
    parent, rank, name, common = load_tax(raw)
    under = top2_under(reviewed(), parent)
    pre = fetch.prebuilt_species(raw)
    rows = []
    for t in pre["species_taxid"].to_list():
        lin = fetch.ancestors(t, parent)
        a, cand, e = nearest(t, lin, under)
        rows.append({"species_taxid": t, "shared_ancestor": a, "relative_taxid": cand, "relative_entries": e})
    pl.DataFrame(rows).write_parquet(OUT / "relatives.parquet")
    print(f"relatives: {len(rows)} pre-built species")


def slug(name):
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def lineage_list(t, parent, rank, name):
    return [{"taxid": a, "rank": rank.get(a), "name": name.get(a)} for a in fetch.ancestors(t, parent)]


def cmd_tables(raw):
    parent, rank, name, common = load_tax(raw)
    allsp = pl.read_parquet(OUT / "swissprot_per_species.parquet")
    entries = dict(allsp.iter_rows())
    rev = reviewed()
    under = top2_under(rev, parent)
    # Reviewed entries in every order and class that holds any reviewed species.
    rank_entries = {}
    for t, e in allsp.iter_rows():
        for a in fetch.ancestors(t, parent):
            if rank.get(a) in RANKS:
                rank_entries[a] = rank_entries.get(a, 0) + e
    mrca = pl.read_parquet(OUT / "timetree_mrca.parquet")
    pairs = {(a, b): m for a, b, m in pl.read_parquet(OUT / "timetree_pairs.parquet").iter_rows()} \
        if (OUT / "timetree_pairs.parquet").exists() else {}
    ages = {r["taxid"]: r for r in mrca.iter_rows(named=True)}
    ebp = dict(pl.read_parquet(OUT / "ebp_species.parquet").iter_rows())
    pre = fetch.prebuilt_species(raw)
    af = {r["upid"]: r for r in pl.read_parquet(OUT / "alphafold_counts.parquet").iter_rows(named=True)} \
        if (OUT / "alphafold_counts.parquet").exists() else {}
    src = json.loads((OUT / "sources.json").read_text())

    recs = []
    for p in pre.iter_rows(named=True):
        t = p["species_taxid"]
        lin = fetch.ancestors(t, parent)
        a, cand, e = nearest(t, lin, under)
        rec = {
            "taxid": t, "name": name.get(t), "common_name": common.get(t), "slug": slug(name.get(t)),
            "lineage": lineage_list(t, parent, rank, name),
            "self_entries": entries.get(t, 0),
            "ranks": {r: next(({"taxid": x, "name": name.get(x), "entries": rank_entries.get(x, 0)}
                               for x in lin if rank.get(x) == r), None) for r in RANKS},
            "relative": None if cand is None else {
                "taxid": cand, "name": name.get(cand), "common_name": common.get(cand), "entries": e,
                "shared_ancestor": {"taxid": a, "name": name.get(a), "rank": rank.get(a)}},
            "proteome": {"upid": p["upid"], "taxid": p["taxid"], "organism": p["organism"],
                         "protein_count": p["protein_count"]},
            "ebp": {"on_list": t in ebp, "projects": ebp.get(t)},
        }
        mya = None
        if cand is not None:
            pa = pairs.get((t, cand))
            if pa is not None:
                mya = {"value": pa, "method": "pairwise", "between": [t, cand]}
            else:
                ag = ages.get(a)
                if ag and ag["mrca_mya"] is not None:
                    mya = {"value": ag["mrca_mya"], "method": "node" if ag["method"] == "node" else "node_pair",
                           "node": a, "pair": [ag["pair_a"], ag["pair_b"]] if ag["method"] == "pair" else None}
        rec["mya"] = mya
        c = af.get(p["upid"])
        if c:
            bands = {k: v for k, v in c.items() if k.startswith("plddt_")}
            total = sum(bands.values())
            rec["alphafold"] = {
                "proteome_entries": c["proteome_entries"], "with_model": c["proteome_entries_with_afdb"],
                "share_with_model": c["proteome_entries_with_afdb"] / c["proteome_entries"] if c["proteome_entries"] else None,
                "afdb_models": c["afdb_models"], "bands5": bands, "models_in_bands": total}
        else:
            rec["alphafold"] = None
        recs.append(rec)
    with open(OUT / "lookups.jsonl", "w") as fh:
        for r in recs:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")

    # Browser tables.
    WEB.mkdir(parents=True, exist_ok=True)
    (WEB / "reviewed.json").write_text(json.dumps({
        "min_entries": fetch.MIN_REVIEWED, "release": src["swissprot_per_taxon"]["release"],
        "under": {str(k): v for k, v in sorted(under.items())},
        "names": {str(t): [name.get(t), common.get(t)] for t in rev["species_taxid"].to_list()},
        "rank_entries": {str(k): v for k, v in sorted(rank_entries.items())},
        "self_entries": {str(t): e for t, e in allsp.iter_rows()},
    }, separators=(",", ":")))
    (WEB / "mrca.json").write_text(json.dumps({
        "release": src["timetree_mrca"]["release"],
        "ages": {str(r["taxid"]): [r["mrca_mya"], r["method"], r["pair_a"], r["pair_b"]]
                 for r in mrca.iter_rows(named=True) if r["mrca_mya"] is not None},
        "names": {str(r["taxid"]): name.get(r["taxid"]) for r in mrca.iter_rows(named=True)},
    }, separators=(",", ":")))
    # Lineages of the organisms the dark-protein share was measured on.
    dark = {"human": 9606, "yeast": 4932, "worm": 6239, "borrelia": 139, "botryllus": 30301}
    (WEB / "dark_lineages.json").write_text(json.dumps(
        {k: fetch.ancestors(t, parent) for k, t in dark.items()}, separators=(",", ":")))
    names_index(raw, recs, name, rank, parent)
    print(f"tables: {len(recs)} pre-built species; {len(under)} nodes with a reviewed species; "
          f"{len(rank_entries)} orders and classes")


def shard_key(word):
    w = re.sub(r"[^a-z]", "", word.lower())[:2]
    return w if len(w) == 2 else "_"


def names_index(raw, recs, name, rank, parent):
    """Autocomplete shards: data/organism/web/names/<two letters>.json. Each entry is
    [label, taxon id, page slug or ""], filed under the first two letters of its scientific
    name, and of each word of a common name, so "oyster" finds "Pacific oyster". Species in
    the index: pre-built ones, the Earth BioGenome list, species with a UniProt reference
    proteome, and species with a Swiss-Prot entry. Any other species is found through NCBI's
    own search, from the page."""
    pre = {r["taxid"]: r["slug"] for r in recs}
    sp = set(pre)
    sp |= set(pl.read_parquet(OUT / "ebp_species.parquet")["taxid"].to_list())
    sp |= set(pl.read_parquet(OUT / "reference_proteomes.parquet")["species_taxid"].drop_nulls().to_list())
    sp |= set(pl.read_parquet(OUT / "swissprot_per_species.parquet")["species_taxid"].to_list())
    sp = {t for t in sp if rank.get(t) == "species" and name.get(t)}
    other = pl.read_parquet(raw / "taxonomy_other_names.parquet").filter(pl.col("taxid").is_in(list(sp)))
    shards: dict[str, list] = {}

    def add(key, entry):
        shards.setdefault(key, []).append(entry)

    for t in sorted(sp):
        sci = name[t]
        add(shard_key(sci), [sci, t, pre.get(t, "")])
    for t, nm, cls in other.iter_rows():
        sci = name[t]
        if cls == "synonym" or cls == "equivalent name":
            if t in pre:  # old names only for pre-built species, to keep the index small
                add(shard_key(nm), [f"{nm} (now {sci})", t, pre[t]])
            continue
        label = f"{nm} ({sci})"
        for w in sorted({shard_key(x) for x in nm.split() if len(re.sub(r'[^a-z]', '', x.lower())) >= 2}):
            add(w, [label, t, pre.get(t, "")])
    d = WEB / "names"
    d.mkdir(parents=True, exist_ok=True)
    for old in d.glob("*.json"):
        old.unlink()
    for k, v in shards.items():
        v.sort(key=lambda e: (e[2] == "", len(e[0]), e[0]))
        (d / f"{k}.json").write_text(json.dumps(v, ensure_ascii=False, separators=(",", ":")))
    print(f"names: {len(sp)} species, {sum(len(v) for v in shards.values())} entries in {len(shards)} shards")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=pathlib.Path, required=True)
    ap.add_argument("what", choices=["relatives", "tables"])
    a = ap.parse_args()
    {"relatives": cmd_relatives, "tables": cmd_tables}[a.what](a.raw)


if __name__ == "__main__":
    main()
