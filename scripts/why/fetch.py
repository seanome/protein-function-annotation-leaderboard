"""Download every input of the "Why this project" numbers and record where each came from.

Each download is written under --raw, and one row per file goes into
results/why/sources.json: URL, local name, bytes, sha256, the time it was fetched, and the
release the server reported (UniProt's X-UniProt-Release header, InterPro's API version,
Catalogue of Life's dataset version, the Last-Modified date of a flat file).

Run one step at a time, for example:
    python3 -I scripts/why/fetch.py --raw ~/data/why-numbers --step uniprot
Steps: uniprot, taxonomy, sifts, interpro, alphafold, col, worms, four.
"""

import argparse
import concurrent.futures as cf
import datetime as dt
import gzip
import hashlib
import http.client
import json
import pathlib
import re
import subprocess
import time
import urllib.parse
import urllib.request

UNIPROT_FTP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase"
PROTEOMES = {"human": ("UP000005640", 9606), "oyster": ("UP000005408", 29159)}
# The four proteins of step 6, resolved by UniProt search on 2026-10-15.
FOUR = {
    "Q25460": "mussel foot protein 1 (Mytilus edulis)",
    "P24856": "antifreeze glycoprotein polyprotein (Notothenia neglecta)",
    "P05484": "omega-conotoxin MVIIA (Conus magus)",
    "P42212": "green fluorescent protein (Aequorea victoria)",
}
COL_DATASET = "3LR"  # Catalogue of Life, latest release
COL_GROUPS = {  # COL taxon id -> name; Bacteria and Archaea are domains in COL
    "N": "Animalia",
    "P": "Plantae",
    "F": "Fungi",
    "C": "Chromista",
    "Z": "Protozoa",
    "CRRY6": "Bacteria",
    "CRLT8": "Archaea",
}
UA = {"User-Agent": "seanome-why-numbers/1.0 (olgabot; protein-function-annotation-leaderboard)"}


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


class Manifest:
    def __init__(self, path):
        self.path = pathlib.Path(path)

    def add(self, name, raw_path, url, release=None, **extra):
        # Re-read before each write: several steps run at once and share this file.
        self.rows = json.loads(self.path.read_text()) if self.path.exists() else {}
        self.rows[name] = {
            "file": raw_path.name,
            "url": url,
            "bytes": raw_path.stat().st_size,
            "sha256": sha256(raw_path),
            "fetched_utc": now(),
            "release": release,
            **extra,
        }
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.rows, indent=1, sort_keys=True) + "\n")


def get(url, tries=10, timeout=300):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.status, dict(resp.headers), resp.read()
        except urllib.error.HTTPError as err:
            if err.code == 404:
                return 404, dict(err.headers), b""
            if attempt == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException):
            if attempt == tries - 1:
                raise
        # InterPro answers 408 under load; wait up to a minute between tries.
        time.sleep(min(60, 2 ** attempt))


def download(url, dest):
    """Large files go through curl: urllib dropped UniProt's chunked stream part way."""
    head = dest.with_name(dest.name + ".headers")
    subprocess.run(
        ["curl", "-sS", "--fail", "--retry", "5", "--retry-all-errors", "-A", UA["User-Agent"],
         "-D", str(head), "-o", str(dest), url],
        check=True,
    )
    # Keep the headers of the last response (after any redirect).
    block = head.read_text().strip().split("\r\n\r\n")[-1]
    headers = {}
    for line in block.splitlines()[1:]:
        key, _, value = line.partition(":")
        headers[key.strip().lower()] = value.strip()
    head.unlink()
    return headers


def uniprot_paged(query, fields, dest):
    """UniProt REST search, 500 rows a page, following the Link header to the end.

    The /stream endpoint dropped the connection part way on every try (HTTP/2 and 1.1),
    so this asks for small pages and checks that every page reports the same release and
    that the row count matches X-Total-Results.
    """
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": query, "fields": fields, "format": "tsv", "size": 500}
    )
    first, releases, n, total = url, set(), 0, None
    with gzip.open(dest, "wt") as out:
        while url:
            status, headers, body = get(url, tries=8)
            headers = {k.lower(): v for k, v in headers.items()}
            releases.add((headers.get("x-uniprot-release"), headers.get("x-uniprot-release-date")))
            total = int(headers["x-total-results"])
            lines = body.decode().splitlines()
            if n == 0:
                out.write(lines[0] + "\n")
            for line in lines[1:]:
                out.write(line + "\n")
                n += 1
            m = re.search(r'<([^>]+)>; rel="next"', headers.get("link", ""))
            url = m.group(1) if m else None
    if len(releases) != 1 or n != total:
        raise RuntimeError(f"{query}: releases {releases}, rows {n}, X-Total-Results {total}")
    (rel, date), = releases
    return first, f"UniProt {rel} ({date})", n


def step_uniprot(raw, man):
    jobs = {
        # Step 2: every UniProtKB entry (Swiss-Prot and TrEMBL) that cites a PDB entry.
        "uniprotkb_with_pdb": (
            "database:pdb",
            "accession,reviewed,organism_id,organism_name,xref_pdb",
        ),
    }
    for tag, (proteome, taxid) in PROTEOMES.items():
        # Step 4: per-protein fields; the canonical set is taken from the FTP FASTA below.
        jobs[f"{tag}_proteome_fields"] = (
            f"proteome:{proteome}",
            "accession,reviewed,length,gene_primary,protein_name,xref_interpro,xref_alphafolddb,xref_pdb",
        )
    for name, (query, fields) in jobs.items():
        dest = raw / f"{name}.tsv.gz"
        url, release, n = uniprot_paged(query, fields, dest)
        man.add(name, dest, url, release, rows=n)
        print(name, release, n, flush=True)
    # The release notes and the one-protein-per-gene FASTA files of the same release.
    flat = {
        "reference_proteomes_readme": f"{UNIPROT_FTP}/reference_proteomes/README",
        "swissprot_fasta": f"{UNIPROT_FTP}/complete/uniprot_sprot.fasta.gz",
        # Steps 1 and 3: every reviewed entry with its OX taxon and DR PDB lines.
        "swissprot_dat": f"{UNIPROT_FTP}/complete/uniprot_sprot.dat.gz",
        "uniprot_relnotes": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt",
    }
    for tag, (proteome, taxid) in PROTEOMES.items():
        flat[f"{tag}_canonical_fasta"] = (
            f"{UNIPROT_FTP}/reference_proteomes/Eukaryota/{proteome}/{proteome}_{taxid}.fasta.gz"
        )
    for name, url in flat.items():
        dest = raw / url.rsplit("/", 1)[1]
        headers = download(url, dest)
        man.add(name, dest, url, f"Last-Modified {headers.get('last-modified')}")
        print(name, dest.stat().st_size)


# InterPro and its member databases, as UniProt cross-references them from InterPro's own
# match data. AntiFam is left out: it flags sequences that are not real proteins.
INTERPRO_XREFS = (
    "xref_interpro,xref_pfam,xref_panther,xref_cdd,xref_gene3d,xref_supfam,xref_smart,"
    "xref_prosite,xref_prints,xref_pirsf,xref_hamap,xref_sfld,xref_ncbifam"
)


def step_members(raw, man):
    """Step 4 InterPro column. The InterPro API's proteome listing timed out (HTTP 408)
    on the human proteome every time, so the matches are read from UniProt."""
    for tag, (proteome, _) in PROTEOMES.items():
        name = f"{tag}_interpro_xrefs"
        dest = raw / f"{name}.tsv.gz"
        url, release, n = uniprot_paged(f"proteome:{proteome}", "accession," + INTERPRO_XREFS, dest)
        man.add(name, dest, url, release, rows=n)
        print(name, release, n, flush=True)


def step_interpro_oyster(raw, man):
    """The InterPro API listing for the oyster only, as a check on the UniProt route."""
    PROTEOMES_SAVED = dict(PROTEOMES)
    PROTEOMES.clear()
    PROTEOMES["oyster"] = PROTEOMES_SAVED["oyster"]
    try:
        step_interpro(raw, man)
    finally:
        PROTEOMES.clear()
        PROTEOMES.update(PROTEOMES_SAVED)


def step_taxonomy(raw, man):
    url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    dest = raw / "taxdump.tar.gz"
    headers = download(url, dest)
    man.add("ncbi_taxdump", dest, url, f"Last-Modified {headers.get('last-modified')}")


def step_sifts(raw, man):
    url = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz"
    dest = raw / "sifts_uniprot_pdb.tsv.gz"
    headers = download(url, dest)
    man.add("sifts_uniprot_pdb", dest, url, f"Last-Modified {headers.get('last-modified')}")


def interpro_release():
    _, _, body = get("https://www.ebi.ac.uk/interpro/api/")
    db = json.loads(body)["databases"]["interpro"]
    return f"InterPro {db['version']} ({db['releaseDate'][:10]})"


def step_interpro(raw, man):
    release = interpro_release()
    for tag, (proteome, _) in PROTEOMES.items():
        # Every UniProtKB protein of the proteome with a match to any InterPro entry or
        # member-database signature (Pfam, PANTHER, CDD, ...), integrated or not.
        url = (
            "https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/all/proteome/uniprot/"
            f"{proteome}/?page_size=200"
        )
        first = url
        dest = raw / f"interpro_{tag}_proteins_with_match.jsonl.gz"
        n = 0
        with gzip.open(dest, "wt") as out:
            while url:
                _, _, body = get(url)
                page = json.loads(body)
                for r in page["results"]:
                    m = r["metadata"]
                    out.write(json.dumps({"accession": m["accession"], "source_database": m["source_database"]}) + "\n")
                    n += 1
                url = page["next"]
        expected = page["count"]
        if n != expected:
            raise RuntimeError(f"{tag}: {n} rows written, API count {expected}")
        man.add(f"interpro_{tag}", dest, first, release, rows=n)
        print(tag, n)
    # Step 6: every match, with coordinates, on the four proteins.
    for acc in FOUR:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}/?page_size=200"
        dest = raw / f"interpro_matches_{acc}.json"
        status, _, body = get(url)
        # 204 (no content) means no match at all; keep an empty result so it is recorded.
        dest.write_bytes(body if status == 200 and body else b'{"count": 0, "results": []}')
        man.add(f"interpro_matches_{acc}", dest, url, release, http_status=status)


def canonical_accessions(raw, tag):
    proteome, taxid = PROTEOMES[tag]
    with gzip.open(raw / f"{proteome}_{taxid}.fasta.gz", "rt") as fh:
        return [line.split("|")[1] for line in fh if line.startswith(">")]


def afdb_one(acc):
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"
    status, _, body = get(url, timeout=60)
    if status == 404:
        return {"accession": acc, "http_status": 404}
    rows = json.loads(body)
    keep = (
        "entryId", "latestVersion", "globalMetricValue", "sequenceChecksum",
        "uniprotStart", "uniprotEnd", "sequence", "uniprotSequence", "modelCreatedDate",
    )
    return {"accession": acc, "http_status": status, "models": [{k: r.get(k) for k in keep} for r in rows]}


def step_alphafold(raw, man):
    accs = {acc: None for tag in PROTEOMES for acc in canonical_accessions(raw, tag)}
    accs.update({acc: None for acc in FOUR})
    dest = raw / "afdb_predictions.jsonl.gz"
    started = now()
    with gzip.open(dest, "wt") as out, cf.ThreadPoolExecutor(12) as pool:
        for i, row in enumerate(pool.map(afdb_one, accs)):
            out.write(json.dumps(row) + "\n")
            if i % 5000 == 0:
                print(i, flush=True)
    man.add(
        "afdb_predictions", dest, "https://alphafold.ebi.ac.uk/api/prediction/{accession}",
        "AlphaFold DB API (version per model in latestVersion)", rows=len(accs), started_utc=started,
    )
    # Per-residue pLDDT of the four proteins, for the residue figure.
    for acc in FOUR:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-confidence_v6.json"
        dest = raw / f"afdb_confidence_{acc}.json"
        status, _, body = get(url)
        dest.write_bytes(body if status == 200 else b"null")
        man.add(f"afdb_confidence_{acc}", dest, url, "AlphaFold DB v6", http_status=status)


def col_count(dataset, params):
    url = f"https://api.checklistbank.org/dataset/{dataset}/nameusage/search?" + urllib.parse.urlencode(params)
    _, _, body = get(url)
    return url, json.loads(body)["total"]


def step_col(raw, man):
    _, _, body = get(f"https://api.checklistbank.org/dataset/{COL_DATASET}")
    meta = json.loads(body)
    key = meta["key"]
    release = f"Catalogue of Life {meta['alias']} (version {meta['version']}, doi {meta['doi']}, dataset key {key})"
    rows = []
    for taxon, name in COL_GROUPS.items():
        url, total = col_count(key, {"TAXON_ID": taxon, "rank": "species", "status": "accepted", "limit": 0})
        rows.append({"col_taxon_id": taxon, "col_group": name, "accepted_species": total, "url": url})
    url, total = col_count(key, {"rank": "species", "status": "accepted", "limit": 0})
    rows.append({"col_taxon_id": None, "col_group": "all", "accepted_species": total, "url": url})
    dest = raw / "col_accepted_species.json"
    dest.write_text(json.dumps({"dataset": meta, "counts": rows}, indent=1))
    man.add("col_accepted_species", dest, f"https://api.checklistbank.org/dataset/{key}", release)


def step_worms(raw, man, names_file):
    """Match metazoan species names to WoRMS and keep its habitat flags."""
    names = [n for n in pathlib.Path(names_file).read_text().splitlines() if n]
    dest = raw / "worms_matches.jsonl.gz"
    base = "https://www.marinespecies.org/rest/AphiaRecordsByMatchNames"
    started = now()

    def batch(chunk):
        q = urllib.parse.urlencode([("scientificnames[]", n) for n in chunk] + [("marine_only", "false")])
        status, _, body = get(f"{base}?{q}", timeout=120)
        hits = json.loads(body) if status == 200 and body else [[] for _ in chunk]
        return [{"query": n, "matches": h or []} for n, h in zip(chunk, hits)]

    chunks = [names[i : i + 50] for i in range(0, len(names), 50)]
    with gzip.open(dest, "wt") as out, cf.ThreadPoolExecutor(4) as pool:
        for rows in pool.map(batch, chunks):
            for row in rows:
                out.write(json.dumps(row) + "\n")
    man.add(
        "worms_matches", dest, base + "?scientificnames[]=...&marine_only=false",
        f"World Register of Marine Species REST API, queried {started}", rows=len(names),
    )


def step_four(raw, man):
    for acc in FOUR:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        dest = raw / f"uniprot_{acc}.json"
        headers = download(url, dest)
        man.add(f"uniprot_{acc}", dest, url, f"UniProt {headers.get('x-uniprot-release')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=pathlib.Path)
    ap.add_argument("--manifest", default="results/why/sources.json")
    ap.add_argument("--step", required=True)
    ap.add_argument("--worms-names", help="species names, one per line (step worms)")
    ap.add_argument("--file", type=pathlib.Path, help="step register: a file made here or on Sherlock")
    ap.add_argument("--name", help="step register: its key in sources.json")
    ap.add_argument("--origin", help="step register: how it was made (command, job id)")
    args = ap.parse_args()
    args.raw.mkdir(parents=True, exist_ok=True)
    man = Manifest(args.manifest)
    if args.step == "worms":
        step_worms(args.raw, man, args.worms_names)
    elif args.step == "register":
        man.add(args.name, args.file, args.origin, "computed", computed=True)
    else:
        globals()[f"step_{args.step}"](args.raw, man)


if __name__ == "__main__":
    main()
