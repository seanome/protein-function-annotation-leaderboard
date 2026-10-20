#!/usr/bin/env python3
"""Download every lookup the organism pages read, and record where each came from.

Lookups only: no search tool runs. Each download is written under --raw (outside the repo);
the small tables the pages read are written to data/organism/, and one row per source goes
into data/organism/sources.json with the URL, the time it was fetched, the release the
server reported, and the sha256 of the raw file.

Run one step at a time, in this order:
    python3 -I scripts/organism/fetch.py --raw ~/data/organism-page --step taxonomy
    python3 -I scripts/organism/fetch.py --raw ~/data/organism-page --step swissprot
    python3 -I scripts/organism/fetch.py --raw ~/data/organism-page --step proteomes
    python3 -I scripts/organism/fetch.py --raw ~/data/organism-page --step ebp
    python3 -I scripts/organism/fetch.py --raw ~/data/organism-page --step timetree
    python3 -I scripts/organism/fetch.py --raw ~/data/organism-page --step alphafold

Steps:
  taxonomy   NCBI taxdump: names and parent of every taxon
  swissprot  reviewed UniProtKB entries per taxon (UniProt REST)
  proteomes  UniProt reference proteomes: accession, taxon, protein count
  ebp        Earth BioGenome Project target list (GoaT, attribute long_list), species only
  timetree   TimeTree age of every ancestor of a species with >= 100 reviewed entries
  pairs      TimeTree pairwise age of each pre-built species and its nearest reviewed
             relative (run scripts/organism/lookups.py first)
  alphafold  per pre-built species: proteome entries with an AlphaFold DB cross-reference
             (UniProt) and AlphaFold DB model counts in 5-point pLDDT bands (AlphaFold DB)
"""

import argparse
import concurrent.futures as cf
import datetime as dt
import gzip
import hashlib
import io
import json
import pathlib
import tarfile
import time
import urllib.parse
import urllib.request

import polars as pl

OUT = pathlib.Path("data/organism")
UA = {"User-Agent": "seanome-organism-page/1.0 (olgabot; protein-function-annotation-leaderboard)"}
TAXDUMP = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
UNIPROT = "https://rest.uniprot.org"
GOAT = "https://goat.genomehubs.org/api/v2"
TIMETREE = "https://timetree.org/api"
AFDB = "https://alphafold.ebi.ac.uk/api/search"
# The five reference species every build pre-builds, whatever the target list says.
REFERENCE = {9606: "human", 10090: "mouse", 4932: "yeast", 562: "E. coli", 3702: "Arabidopsis"}
MIN_REVIEWED = 100
PLDDT_STEP = 5


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def get(url, tries=5, timeout=300):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return dict(r.headers), r.read()
        except Exception as e:  # network errors and 5xx: wait and retry
            if i == tries - 1 or getattr(e, "code", None) in (400, 404):
                raise
            print(f"retry {i + 1} {url}: {e}")
            time.sleep(5 * (i + 1))


def record(name, raw_path, url, release, **extra):
    p = OUT / "sources.json"
    rows = json.loads(p.read_text()) if p.exists() else {}
    rows[name] = {"file": raw_path.name, "url": url, "bytes": raw_path.stat().st_size,
                  "sha256": sha256(raw_path), "fetched_utc": now(), "release": release, **extra}
    OUT.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(rows, indent=1, sort_keys=True) + "\n")


# --------------------------------------------------------------------------- taxonomy

def step_taxonomy(raw):
    f = raw / "taxdump.tar.gz"
    hdr, body = get(TAXDUMP, timeout=1800)
    f.write_bytes(body)
    release = f'NCBI Taxonomy taxdump, Last-Modified {hdr.get("Last-Modified")}'
    with tarfile.open(f) as tf:
        nodes = tf.extractfile("nodes.dmp").read().decode()
        names = tf.extractfile("names.dmp").read().decode()
        merged = tf.extractfile("merged.dmp").read().decode()
    old_new = [tuple(int(x) for x in line.rstrip("\t|").split("\t|\t")[:2]) for line in merged.splitlines()]
    pl.DataFrame(old_new, schema=["old", "new"], orient="row").write_parquet(raw / "taxonomy_merged.parquet")
    tid, parent, rank = [], [], []
    for line in nodes.splitlines():
        c = line.split("\t|\t")
        tid.append(int(c[0]))
        parent.append(int(c[1]))
        rank.append(c[2])
    sci, common = {}, {}
    for line in names.splitlines():
        c = line.rstrip("\t|").split("\t|\t")
        t, name, cls = int(c[0]), c[1], c[3]
        if cls == "scientific name":
            sci[t] = name
        elif cls in ("genbank common name", "common name"):
            # genbank common name first; another common name only if there is none
            if cls == "genbank common name" or t not in common:
                common[t] = name
    nodes_df = pl.DataFrame({"taxid": tid, "parent": parent, "rank": rank}).with_columns(
        pl.col("taxid").replace_strict(sci, default=None, return_dtype=pl.Utf8).alias("name"),
        pl.col("taxid").replace_strict(common, default=None, return_dtype=pl.Utf8).alias("common_name"),
    )
    # Every common name too, for the autocomplete (a species can have several).
    allnames = []
    for line in names.splitlines():
        c = line.rstrip("\t|").split("\t|\t")
        if c[3] in ("genbank common name", "common name", "synonym", "equivalent name"):
            allnames.append((int(c[0]), c[1], c[3]))
    alt = pl.DataFrame(allnames, schema=["taxid", "name", "name_class"], orient="row")
    OUT.mkdir(parents=True, exist_ok=True)
    (raw / "taxonomy_nodes.parquet").unlink(missing_ok=True)
    nodes_df.write_parquet(raw / "taxonomy_nodes.parquet")
    alt.write_parquet(raw / "taxonomy_other_names.parquet")
    record("ncbi_taxdump", f, TAXDUMP, release, rows=nodes_df.height)
    print(f"taxonomy: {nodes_df.height} taxa, {len(common)} with a common name")


def lineage_table(raw):
    """taxid -> list of ancestors, nearest first (the taxon itself first)."""
    n = pl.read_parquet(raw / "taxonomy_nodes.parquet")
    parent = dict(zip(n["taxid"].to_list(), n["parent"].to_list()))
    m = pl.read_parquet(raw / "taxonomy_merged.parquet")
    MERGED.update(zip(m["old"].to_list(), m["new"].to_list()))
    return parent, n


# Taxon ids NCBI has merged into another id: old -> new.
MERGED: dict[int, int] = {}


def ancestors(t, parent):
    """The taxon and its ancestors, nearest first. A merged id is followed to its new id; an
    id NCBI no longer has gives an empty list."""
    t = MERGED.get(t, t)
    if t not in parent:
        return []
    out = [t]
    while parent.get(t, 1) != t and t != 1:
        t = parent[t]
        out.append(t)
    return out


def species_of(t, parent, rank):
    """The species a taxon belongs to (itself if it is one), or None above species level."""
    for a in ancestors(t, parent):
        if rank.get(a) == "species":
            return a
    return None


# --------------------------------------------------------------------------- Swiss-Prot

def step_swissprot(raw):
    url = f"{UNIPROT}/uniprotkb/stream?query=reviewed:true&fields=organism_id&format=tsv&compressed=true"
    hdr, body = get(url, timeout=1800)
    f = raw / "swissprot_organism_id.tsv.gz"
    f.write_bytes(body)
    release = f'UniProt {hdr.get("X-UniProt-Release")} ({hdr.get("X-UniProt-Release-Date")})'
    t = pl.read_csv(io.BytesIO(gzip.decompress(body)), separator="\t")
    t = t.rename({t.columns[0]: "taxid"}).group_by("taxid").len("entries")
    parent, n = lineage_table(raw)
    rank = dict(zip(n["taxid"].to_list(), n["rank"].to_list()))
    t = t.with_columns(pl.col("taxid").map_elements(lambda x: species_of(x, parent, rank),
                                                    return_dtype=pl.Int64).alias("species_taxid"))
    no_sp = int(t.filter(pl.col("species_taxid").is_null())["entries"].sum())
    sp = t.drop_nulls("species_taxid").group_by("species_taxid").agg(pl.col("entries").sum()).sort("species_taxid")
    sp.write_parquet(OUT / "swissprot_per_species.parquet")
    t.select("taxid", "entries").sort("taxid").write_parquet(OUT / "swissprot_per_taxid.parquet")
    record("swissprot_per_taxon", f, url, release, entries=int(t["entries"].sum()),
           entries_without_species=no_sp, species=sp.height)
    print(f"swissprot: {int(t['entries'].sum())} entries, {sp.height} species, {no_sp} entries above species level")


# --------------------------------------------------------------------------- proteomes

def step_proteomes(raw):
    url = (f"{UNIPROT}/proteomes/stream?query=reference:true"
           "&fields=upid,organism_id,organism,protein_count,busco&format=tsv")
    hdr, body = get(url, timeout=1800)
    f = raw / "reference_proteomes.tsv"
    f.write_bytes(body)
    release = f'UniProt {hdr.get("X-UniProt-Release")} ({hdr.get("X-UniProt-Release-Date")})'
    t = pl.read_csv(f, separator="\t", infer_schema_length=0)
    t = t.rename({"Proteome Id": "upid", "Organism Id": "taxid", "Organism": "organism",
                  "Protein count": "protein_count", "BUSCO": "busco"}).with_columns(
        pl.col("taxid").cast(pl.Int64), pl.col("protein_count").cast(pl.Int64))
    parent, n = lineage_table(raw)
    rank = dict(zip(n["taxid"].to_list(), n["rank"].to_list()))
    t = t.with_columns(pl.col("taxid").map_elements(lambda x: species_of(x, parent, rank),
                                                    return_dtype=pl.Int64).alias("species_taxid"))
    # Which kingdom folder of the UniProt FTP site holds the FASTA file.
    sk = {2759: "Eukaryota", 2: "Bacteria", 2157: "Archaea", 10239: "Viruses"}

    def folder(x):
        for a in ancestors(x, parent):
            if a in sk:
                return sk[a]
        return None
    t = t.with_columns(pl.col("taxid").map_elements(folder, return_dtype=pl.Utf8).alias("ftp_folder"))
    t.sort("upid").write_parquet(OUT / "reference_proteomes.parquet")
    record("uniprot_reference_proteomes", f, url, release, rows=t.height)
    print(f"proteomes: {t.height} reference proteomes")


# --------------------------------------------------------------------------- EBP list

def goat_count(q):
    u = f"{GOAT}/count?query={urllib.parse.quote(q)}&result=taxon&taxonomy=ncbi"
    return json.loads(get(u)[1])["count"]


def goat_search(q, size):
    u = (f"{GOAT}/search?query={urllib.parse.quote(q)}&result=taxon&taxonomy=ncbi"
         f"&fields=long_list&size={size}&includeEstimates=false")
    hdr, body = get(u, timeout=900)
    return json.loads(body)


def step_ebp(raw):
    """GoaT holds the target lists of the Earth BioGenome Project and its affiliated projects
    in the attribute long_list. A species is on the list when long_list has any value. GoaT's
    search returns at most 50,000 rows a call, so the tree is split into clades below that."""
    parent, n = lineage_table(raw)
    children = {}
    for t, p in zip(n["taxid"].to_list(), n["parent"].to_list()):
        if t != p:
            children.setdefault(p, []).append(t)
    base = "long_list AND tax_rank(species)"
    rows, todo, calls, no_ncbi = {}, [1], [], []
    while todo:
        c = todo.pop()
        q = f"tax_tree({c}) AND {base}"
        k = goat_count(q)
        if k == 0:
            continue
        if k > 40000:
            todo.extend(children.get(c, []))
            continue
        d = goat_search(q, max(k, 1))
        calls.append({"clade": c, "count": k})
        for r in d["results"]:
            raw_id = r["result"]["taxon_id"]
            if not str(raw_id).isdigit():
                # GoaT also lists species that NCBI has no taxon id for yet (ids such as
                # "pyProAcut"); they cannot be matched to UniProt, so they are counted, not kept.
                no_ncbi.append(raw_id)
                continue
            tid = int(raw_id)
            ll = r["result"].get("fields", {}).get("long_list", {}).get("value")
            rows[tid] = ",".join(ll) if isinstance(ll, list) else ll
        print(f"ebp: clade {c}: {k} species, total {len(rows)}")
    total = goat_count(base)
    f = raw / "goat_ebp_species.json"
    f.write_text(json.dumps({"calls": calls, "species": rows}))
    t = pl.DataFrame({"taxid": list(rows), "long_list": list(rows.values())}).sort("taxid")
    t.write_parquet(OUT / "ebp_species.parquet")
    hdr, body = get(f"{GOAT}/search?query={urllib.parse.quote('tax_eq(9606)')}&result=taxon&taxonomy=ncbi&size=1")
    index = json.loads(body)["results"][0]["index"]
    record("goat_ebp_list", f, f"{GOAT}/search?query=tax_tree(<clade>) AND {base}", f"GoaT index {index}",
           species=t.height, goat_total=total, species_without_ncbi_id=len(set(no_ncbi)))
    if t.height + len(set(no_ncbi)) != total:
        # Seen 2026-10-16: 17 of 160,083 missing, GoaT-only species filed directly under a
        # clade the download split into its NCBI children. Reported, not fatal.
        print(f"ebp: WARNING GoaT says {total} species, the split download holds {t.height} + {len(set(no_ncbi))}")
    print(f"ebp: {t.height} species")


# --------------------------------------------------------------------------- TimeTree

def mrca_age(taxid):
    try:
        hdr, body = get(f"{TIMETREE}/mrca/id/{taxid}", tries=3, timeout=60)
    except Exception:  # TimeTree answers HTTP 500 for some taxa it cannot place
        return None
    lines = body.decode().strip().splitlines()
    if len(lines) < 2 or not lines[0].startswith("taxon_id"):
        return None
    v = lines[1].split(",")
    try:
        age = float(v[1])
    except (ValueError, IndexError):
        return None
    return age if age > 0 else None


def pairwise_age(a, b):
    """TimeTree's age for the split between two species, or None when TimeTree has none."""
    try:
        hdr, body = get(f"{TIMETREE}/pairwise/{a}/{b}", tries=3, timeout=60)
    except Exception:  # TimeTree answers HTTP 500 for some pairs it cannot place
        return None
    lines = body.decode().strip().splitlines()
    if len(lines) < 2 or not lines[0].startswith("taxon_a_id"):
        return None
    head, v = lines[0].split(","), lines[1].split(",")
    try:
        # TimeTree's adjusted age when it has one, else its median of the studies
        # (precomputed_age). Neurospora crassa and N. intermedia: adjusted 0, precomputed 2.97.
        # Both 0 (E. coli and E. fergusonii) means TimeTree has no estimate.
        adj, pre = float(v[head.index("adjusted_age")] or 0), float(v[head.index("precomputed_age")] or 0)
        age = adj if adj > 0 else pre
    except (ValueError, IndexError):
        return None
    return age if age > 0 else None


def step_timetree(raw):
    parent, n = lineage_table(raw)
    sp = pl.read_parquet(OUT / "swissprot_per_species.parquet").filter(pl.col("entries") >= MIN_REVIEWED)
    nodes = set()
    for t in sp["species_taxid"].to_list():
        nodes.update(ancestors(t, parent)[1:])
    nodes.discard(1)
    nodes = sorted(nodes)
    print(f"timetree: {len(nodes)} ancestor nodes of {sp.height} species with >= {MIN_REVIEWED} reviewed entries")
    with cf.ThreadPoolExecutor(6) as ex:
        ages = dict(zip(nodes, ex.map(mrca_age, nodes)))
    # A node TimeTree has no age for: the pairwise age of two reviewed species that split
    # there (the largest species under each of its two largest child branches).
    entries = dict(zip(sp["species_taxid"].to_list(), sp["entries"].to_list()))
    lin = {t: ancestors(t, parent) for t in entries}
    pair_for = {}
    for node in [x for x in nodes if ages[x] is None]:
        by_child = {}
        for t, a in lin.items():
            if node in a[1:]:
                child = a[a.index(node) - 1]
                if child not in by_child or entries[t] > entries[by_child[child]]:
                    by_child[child] = t
        best = sorted(by_child.values(), key=lambda t: (-entries[t], t))[:2]
        if len(best) == 2:
            pair_for[node] = tuple(best)
    with cf.ThreadPoolExecutor(6) as ex:
        pages = dict(zip(pair_for, ex.map(lambda ab: pairwise_age(*ab), pair_for.values())))
    hdr, body = get(f"{TIMETREE}/taxon/9606")
    version = json.loads(body).get("version")
    method = ["node" if ages[x] is not None else ("pair" if pages.get(x) else None) for x in nodes]
    t = pl.DataFrame({
        "taxid": nodes,
        "mrca_mya": [ages[x] if ages[x] is not None else pages.get(x) for x in nodes],
        "method": method,
        "pair_a": [pair_for[x][0] if m == "pair" else None for x, m in zip(nodes, method)],
        "pair_b": [pair_for[x][1] if m == "pair" else None for x, m in zip(nodes, method)],
    })
    f = raw / "timetree_mrca.json"
    f.write_text(json.dumps({"node": {str(k): v for k, v in ages.items()},
                             "pair": {str(k): [list(pair_for[k]), v] for k, v in pages.items()}}))
    t.write_parquet(OUT / "timetree_mrca.parquet")
    record("timetree_mrca", f, f"{TIMETREE}/mrca/id/<NCBI taxon id>",
           f"TimeTree 5, database version {version}", nodes=len(nodes),
           nodes_with_age=int(t["mrca_mya"].is_not_null().sum()))
    print(f"timetree: {int(t['mrca_mya'].is_not_null().sum())} of {len(nodes)} nodes have an age; version {version}")


# --------------------------------------------------------------------------- AlphaFold

def prebuilt_species(raw):
    """Species pre-built as static pages: the five reference species and every species on the
    Earth BioGenome target list that has a UniProt reference proteome."""
    ebp = set(pl.read_parquet(OUT / "ebp_species.parquet")["taxid"].to_list())
    prot = pl.read_parquet(OUT / "reference_proteomes.parquet").drop_nulls("species_taxid")
    keep = prot.filter(pl.col("species_taxid").is_in(list(ebp | set(REFERENCE))))
    return choose_proteome(keep)


def choose_proteome(prot):
    """One reference proteome per species: the one whose organism has most reviewed entries
    (E. coli K-12 rather than another strain), then most proteins, then lowest accession.
    docs/assets/organism.js applies the same rule to a species looked up in the browser."""
    per_taxid = pl.read_parquet(OUT / "swissprot_per_taxid.parquet").rename({"entries": "reviewed_entries"})
    return (prot.join(per_taxid, on="taxid", how="left").with_columns(pl.col("reviewed_entries").fill_null(0))
            .sort(["species_taxid", "reviewed_entries", "protein_count", "upid"], descending=[False, True, True, False])
            .unique("species_taxid", keep="first").sort("species_taxid"))


def afdb_count(q):
    """Number of AlphaFold DB models matching a search. The API answers HTTP 404 when nothing
    matches (checked 2026-10-16 with taxId:999999999), so 404 is a count of 0."""
    try:
        return json.loads(get(f"{AFDB}?q={urllib.parse.quote(q)}&type=main&rows=0")[1])["numFound"]
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return 0
        raise


def af_counts(row):
    upid, taxid = row["upid"], row["taxid"]
    out = {"upid": upid}
    h, _ = get(f"{UNIPROT}/uniprotkb/search?query={urllib.parse.quote(f'proteome:{upid} AND database:alphafolddb')}"
               "&size=0&fields=accession", timeout=120)
    out["proteome_entries_with_afdb"] = int(h["X-Total-Results"])
    h, _ = get(f"{UNIPROT}/uniprotkb/search?query={urllib.parse.quote(f'proteome:{upid}')}&size=0&fields=accession",
               timeout=120)
    out["proteome_entries"] = int(h["X-Total-Results"])
    out["uniprot_release"] = h.get("X-UniProt-Release")
    base = f"taxId:{taxid} AND isUniProtReferenceProteome:true"
    out["afdb_models"] = afdb_count(base)
    for lo in range(0, 100, PLDDT_STEP):
        hi = lo + PLDDT_STEP
        # Solr ranges: [lo TO hi} includes lo, excludes hi; the last band includes 100.
        rng = f"[{lo} TO {hi}]" if hi == 100 else f"[{lo} TO {hi}}}"
        q = f"{base} AND globalMetricValue:{rng}"
        out[f"plddt_{lo:02d}_{hi:03d}"] = afdb_count(q)
    return out


def step_alphafold(raw):
    pre = prebuilt_species(raw)
    done_f = raw / "alphafold_counts.jsonl"
    done = {}
    if done_f.exists():
        for line in done_f.read_text().splitlines():
            r = json.loads(line)
            if "error" not in r:
                done[r["upid"]] = r
    todo = [r for r in pre.iter_rows(named=True) if r["upid"] not in done]
    print(f"alphafold: {pre.height} pre-built species, {len(todo)} still to fetch")

    def safe(row):
        try:
            return af_counts(row)
        except Exception as e:  # one failed species is retried on the next run, not fatal
            return {"upid": row["upid"], "error": str(e)}
    with open(done_f, "a") as fh, cf.ThreadPoolExecutor(6) as ex:
        for i, r in enumerate(ex.map(safe, todo)):
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            if i % 100 == 0:
                print(f"alphafold: {i + 1} of {len(todo)}")
    rows = [json.loads(x) for x in done_f.read_text().splitlines()]
    errors = {r["upid"] for r in rows if "error" in r} - {r["upid"] for r in rows if "error" not in r}
    if errors:
        print(f"alphafold: {len(errors)} species failed; run the step again to retry them")
    t = pl.DataFrame([r for r in rows if "error" not in r]).unique("upid", keep="last").sort("upid")
    t.write_parquet(OUT / "alphafold_counts.parquet")
    record("alphafold_counts", done_f, f"{AFDB}?q=taxId:<taxon> AND isUniProtReferenceProteome:true "
           f"AND globalMetricValue:[lo TO hi}} (and UniProt proteome:<UP> AND database:alphafolddb)",
           "AlphaFold DB search API (models version 6) and UniProt " + str(t["uniprot_release"][0]),
           species=t.height, band_width=PLDDT_STEP)
    print(f"alphafold: {t.height} species")


def step_pairs(raw):
    """TimeTree pairwise age between each pre-built species and its nearest reviewed relative
    (chosen by scripts/organism/lookups.py, which writes data/organism/relatives.parquet)."""
    rel = pl.read_parquet(OUT / "relatives.parquet").drop_nulls("relative_taxid")
    pairs = sorted({(a, b) for a, b in rel.select("species_taxid", "relative_taxid").iter_rows()})
    done_f = raw / "timetree_pairs.jsonl"
    done = set()
    if done_f.exists():
        for line in done_f.read_text().splitlines():
            r = json.loads(line)
            if r["mya"] is not None:  # a pair with no age is asked again on the next run
                done.add((r["a"], r["b"]))
    todo = [p for p in pairs if p not in done]
    print(f"pairs: {len(pairs)} pairs, {len(todo)} to fetch")
    with open(done_f, "a") as fh, cf.ThreadPoolExecutor(6) as ex:
        for (a, b), age in zip(todo, ex.map(lambda ab: pairwise_age(*ab), todo)):
            fh.write(json.dumps({"a": a, "b": b, "mya": age}) + "\n")
    rows = [json.loads(x) for x in done_f.read_text().splitlines()]
    t = pl.DataFrame(rows, schema={"a": pl.Int64, "b": pl.Int64, "mya": pl.Float64}) \
        .unique(["a", "b"], keep="last").sort("a", "b")
    t.write_parquet(OUT / "timetree_pairs.parquet")
    hdr, body = get(f"{TIMETREE}/taxon/9606")
    record("timetree_pairs", done_f, f"{TIMETREE}/pairwise/<species>/<relative>",
           f'TimeTree 5, database version {json.loads(body).get("version")}', pairs=t.height,
           pairs_with_age=int(t["mya"].is_not_null().sum()))
    print(f"pairs: {int(t['mya'].is_not_null().sum())} of {t.height} have an age")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=pathlib.Path, required=True)
    ap.add_argument("--step", required=True,
                    choices=["taxonomy", "swissprot", "proteomes", "ebp", "timetree", "pairs", "alphafold"])
    a = ap.parse_args()
    a.raw.mkdir(parents=True, exist_ok=True)
    globals()[f"step_{a.step}"](a.raw)


if __name__ == "__main__":
    main()
