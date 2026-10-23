"""UniProt and TimeTree lookups for the "solvers unreadable" test, read on the day they run.

Shared by model_organisms.py, atlas_proteomes.py and draw_controls.py. Standard library
only. Every TimeTree answer is cached in a JSON file next to the outputs, with the URL and
the UTC time it was read, so a second agent can redo the control draw from the same answers.
"""

import concurrent.futures
import csv
import gzip
import io
import json
import pathlib
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "seanome-protein-function-annotation-leaderboard (olga@seanome.org)"}
REST = "https://rest.uniprot.org"
FTP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase"
README_URL = FTP + "/reference_proteomes/README"
TIMETREE = "https://timetree.temple.edu/api/pairwise/9606/{taxid}"
HUMAN_TAXID = 9606

# BUSCO dataset by kingdom, PLAN.md "BUSCO lineage, chosen by kingdom before any run".
# The first lineage name found in UniProt's taxonomic lineage wins.
BUSCO_BY_KINGDOM = (("Metazoa", "metazoa_odb12"), ("Fungi", "fungi_odb12"), ("Viridiplantae", "viridiplantae_odb12"))
BUSCO_ALL = "eukaryota_odb12"


def now_utc():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def get(url, tries=6):
    """GET with retries; returns (headers as a lower-case dict, body text)."""
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                return {k.lower(): v for k, v in r.headers.items()}, r.read().decode()
        except Exception as e:  # network errors and HTTP 5xx: wait and retry
            if i == tries - 1:
                raise RuntimeError(f"{url}: {e}") from e
            time.sleep(5 * (i + 1))


def tsv_rows(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def busco_lineage(taxonomic_lineage):
    names = [s.strip() for s in taxonomic_lineage.split(",")]
    for kingdom, dataset in BUSCO_BY_KINGDOM:
        if kingdom in names:
            return dataset
    return BUSCO_ALL


def reference_proteomes(query):
    """UniProt proteomes search restricted to reference proteomes; returns (rows, release)."""
    url = f"{REST}/proteomes/search?" + urllib.parse.urlencode(
        {"query": f"reference:true AND ({query})", "fields": "upid,organism,organism_id,protein_count,lineage,busco", "format": "tsv", "size": 500}
    )
    headers, body = get(url)
    release = f"UniProt {headers.get('x-uniprot-release')} ({headers.get('x-uniprot-release-date')})"
    return tsv_rows(body), release, url


def species_taxid(binomial):
    """UniProt taxonomy id of a species named by binomial, and a note.

    First an exact scientific-name match at rank species; failing that, a species whose
    UniProt synonyms include the binomial exactly (e.g. Rana sylvatica -> Lithobates
    sylvaticus). Returns (None, note) when neither gives exactly one species."""
    url = f"{REST}/taxonomy/search?" + urllib.parse.urlencode(
        {"query": f'scientific:"{binomial}" AND rank:species', "fields": "id,scientific_name,rank", "format": "tsv"}
    )
    _, body = get(url)
    hits = [r for r in tsv_rows(body) if r["Scientific name"].lower() == binomial.lower()]
    if len(hits) == 1:
        return int(hits[0]["Taxon Id"]), ""
    url = f"{REST}/taxonomy/search?" + urllib.parse.urlencode(
        {"query": f'"{binomial}" AND rank:species', "fields": "id,scientific_name,rank,synonyms", "format": "tsv"}
    )
    _, body = get(url)
    syn = [r for r in tsv_rows(body) if binomial.lower() in [x.strip().lower() for x in (r.get("Synonyms") or "").split(",")]]
    if len(syn) == 1:
        return int(syn[0]["Taxon Id"]), f"matched as a UniProt synonym of {syn[0]['Scientific name']}"
    return None, f"binomial not found as one UniProt species ({len(hits)} exact, {len(syn)} synonym matches)"


def readme_table(dest):
    """The per-proteome table at the end of the reference proteomes README, saved to dest.

    Returns {upid: row} with the README's own column names, and the Last-Modified header.
    The table's release line ("Release 2026_03, 02-Sep-2026") is in the saved copy.
    """
    headers, body = get(README_URL)
    with gzip.open(dest, "wt") as fh:  # dest ends in .gz; 2.6 MB of text, kept for the redo
        fh.write(body)
    lines = body.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("Proteome_ID\t"))
    rows = tsv_rows("\n".join(lines[start:]))
    return {r["Proteome_ID"]: r for r in rows if r.get("Tax_ID")}, headers.get("last-modified")


def binomial_of(species_name):
    """'Magallana gigas (Pacific oyster) (Crassostrea gigas)' -> 'Magallana gigas'."""
    return " ".join(species_name.split("(")[0].split()[:2])


class TimeTree:
    """TimeTree pairwise times from human, cached in a JSON file (one entry per taxid)."""

    def __init__(self, cache_path):
        self.path = pathlib.Path(cache_path)
        self.cache = json.loads(self.path.read_text()) if self.path.exists() else {}

    @staticmethod
    def _read(taxid):
        url = TIMETREE.format(taxid=taxid)
        rec = {"url": url, "read_utc": now_utc()}
        if int(taxid) == HUMAN_TAXID:
            rec.update(adjusted_age=0.0, precomputed_age=0.0, name_b="Homo sapiens", note="human itself, set to 0")
            return rec
        _, body = get(url)
        rows = list(csv.DictReader(io.StringIO(body))) if body.strip() else []
        if rows and rows[0].get("adjusted_age"):
            r = rows[0]
            rec.update(
                adjusted_age=float(r["adjusted_age"]),
                precomputed_age=float(r["precomputed_age"]) if r.get("precomputed_age") else None,
                name_b=r["scientific_name_b"],
                studies=int(r["all_total"]) if r.get("all_total") else None,
            )
        else:
            rec.update(adjusted_age=None, precomputed_age=None, note="TimeTree returned no time")
        return rec

    def lookup(self, taxid):
        """Returns the cached record: adjusted_age and precomputed_age in My, or None if
        TimeTree gave no row. Human against itself is 0 My (TimeTree returns no row)."""
        key = str(taxid)
        if key not in self.cache:
            self.cache[key] = self._read(taxid)
            self.save()
        return self.cache[key]

    def prefetch(self, taxids, workers=4):
        """Read the uncached taxids, a few at a time in parallel. The answers do not depend
        on the order they are read in, so the draw is the same as one lookup at a time."""
        todo = list(dict.fromkeys(str(t) for t in taxids if str(t) not in self.cache))
        if not todo:
            return
        with concurrent.futures.ThreadPoolExecutor(workers) as ex:
            for key, rec in zip(todo, ex.map(self._read, todo)):
                self.cache[key] = rec
        self.save()

    def save(self):
        self.path.write_text(json.dumps(self.cache, indent=1, sort_keys=True) + "\n")
