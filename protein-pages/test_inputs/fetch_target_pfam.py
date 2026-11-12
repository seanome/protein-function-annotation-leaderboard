"""Pfam domains, with positions, of the Swiss-Prot proteins drawn on the test pages for which
neither hmmscan.domtblout nor interproscan.tsv gives any domain. They come from the InterPro
REST API, one request per protein:

    https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/<accession>/

(InterPro API docs: https://github.com/ProteinsWebTeam/interpro7-api/tree/master/docs). The
answer lists each Pfam entry matching the protein, with its start and end on it. A protein with
no Pfam match comes back as HTTP 204 (no content); it is written as one row with Pfam id "-",
so build.py can say "none" instead of "not known".

An isoform (P12345-2) is not a protein of that endpoint; its matches come from

    https://www.ebi.ac.uk/interpro/api/protein/uniprot/<canonical accession>/?isoforms=<isoform>

whose "features" hold every member-database match on the isoform; the Pfam ones are kept.

Usage (from protein-pages/):

    python test_inputs/fetch_target_pfam.py test_inputs/data

Writes test_inputs/data/target_pfam_interpro.tsv. Every answer is cached in
test_inputs/interpro_cache.jsonl (not committed), so a rerun only asks for what is missing.
"""
import concurrent.futures
import csv
import datetime
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import build  # noqa: E402

API = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/%s/"
API_ISOFORM = "https://www.ebi.ac.uk/interpro/api/protein/uniprot/%s/?isoforms=%s"
OUT_NAME = "target_pfam_interpro.tsv"
COLUMNS = ["accession", "pfam", "pfam_name", "start", "end", "source"]


def missing_targets(data):
    """Swiss-Prot accessions of drawn search hits (E-value at most 10) that have no domain in
    hmmscan.domtblout or interproscan.tsv, read with build.py's own readers."""
    data = Path(data)

    def inp(name):
        p = data / name
        return p if p.exists() else data / (name + ".gz")
    prots = build.read_fasta(inp("proteins.fasta"))
    prots = {h.split()[0]: s for h, s in prots.items()}
    targets, known = set(), set()
    for t in build.TRACKS:
        p = inp(t["file"])
        if not p.exists():
            continue
        if t["kind"] == "m8":
            got = build.parse_m8(p, t["key"])
        elif t["kind"] == "domtbl_seq":
            got = build.parse_domtbl_seq(p, t["key"])
        elif t["kind"] == "kmerseek":
            got, _ = build.parse_kmerseek(p)
        elif t["kind"] == "domtbl_pfam":
            _, tp, _ = build.parse_domtbl_pfam(p, prots)
            known.update(tp)
            continue
        elif t["kind"] == "interproscan":
            _, ta, _ = build.parse_interproscan(p, prots)
            known.update(ta)
            continue
        else:
            continue
        for q, hs in got.items():
            if q not in prots:
                continue
            for x in hs:
                h = build.materialize(x)
                if h["passes"] and h["target"]:
                    targets.add(h["target"])
    return sorted(targets - known)


def fetch(acc, tries=8):
    """(accession, HTTP status, parsed JSON or None). Retries on errors and 408/429/5xx; on 429
    (too many requests) it waits as long as the server's Retry-After asks, or 10 s."""
    url = API_ISOFORM % (acc.split("-")[0], acc) if "-" in acc else API % acc
    for i in range(tries):
        wait = 2 * 2 ** min(i, 5)
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                       "User-Agent": "protein-pages test-input builder"})
            with urllib.request.urlopen(req, timeout=30) as r:
                if r.status == 204:
                    return acc, 204, None
                body = r.read()
                return acc, r.status, json.loads(body) if body else None
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):
                return acc, e.code, None
            if e.code not in (408, 429, 500, 502, 503, 504):
                return acc, e.code, None
            if e.code == 429:
                ra = e.headers.get("Retry-After", "")
                wait = int(ra) if ra.isdigit() else 10
        except (urllib.error.URLError, TimeoutError, ConnectionError, json.JSONDecodeError):
            pass
        time.sleep(wait)
    return acc, None, None


def rows_from(acc, status, js):
    if status == 204:
        return [[acc, "-", "", "", "", "InterPro API: no Pfam match"]]
    if status == 404:
        return [[acc, "-", "", "", "", "InterPro API: protein not found"]]
    out = []
    for f in ((js or {}).get("features") or {}).values():  # an isoform's answer
        if f.get("source_database") == "pfam":
            for loc in f.get("locations") or []:
                for fr in loc.get("fragments") or []:
                    out.append([acc, f["accession"], f.get("name") or "", fr["start"], fr["end"], "InterPro API, isoform"])
    for res in (js or {}).get("results", []):
        md = res["metadata"]
        for p in res.get("proteins", []):
            for loc in p.get("entry_protein_locations") or []:
                for fr in loc.get("fragments") or []:
                    out.append([acc, md["accession"], md.get("name") or "", fr["start"], fr["end"], "InterPro API"])
    if not out:
        out = [[acc, "-", "", "", "", "InterPro API: no Pfam match"]]
    return sorted(out, key=lambda r: (r[3], r[4], r[1]))


def main():
    data = Path(sys.argv[1] if len(sys.argv) > 1 else HERE / "data")
    cache_path = HERE / "interpro_cache.jsonl"
    cache = {}
    if cache_path.exists():
        for line in cache_path.read_text().splitlines():
            d = json.loads(line)
            if d["status"] in (200, 204, 404):
                cache[d["acc"]] = d
    accs = missing_targets(data)
    todo = [a for a in accs if a not in cache]
    print("%d targets have no domain data; %d already cached, %d to ask" % (len(accs), len(accs) - len(todo), len(todo)), flush=True)
    with open(cache_path, "a") as fh, concurrent.futures.ThreadPoolExecutor(3) as pool:
        for n, (acc, status, js) in enumerate(pool.map(fetch, todo), 1):
            d = dict(acc=acc, status=status, json=js)
            fh.write(json.dumps(d, separators=(",", ":")) + "\n")
            fh.flush()
            if status in (200, 204, 404):
                cache[acc] = d
            if n % 250 == 0:
                print("  %d of %d" % (n, len(todo)), flush=True)
    failed = [a for a in accs if a not in cache]
    rows = []
    for a in accs:
        if a in cache:
            rows.extend(rows_from(a, cache[a]["status"], cache[a]["json"]))
    with open(data / OUT_NAME, "w", newline="") as fh:
        fh.write("# Pfam domains of Swiss-Prot targets from the InterPro REST API (%s; isoforms: %s), read %s "
                 "by fetch_target_pfam.py. pfam '-' means InterPro lists no Pfam match.\n"
                 % (API % "<accession>", API_ISOFORM % ("<accession>", "<isoform>"), datetime.date.today().isoformat()))
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(rows)
    n_dom = len({r[0] for r in rows if r[1] != "-"})
    print("wrote %d rows: %d proteins with a Pfam match, %d with none, %d not answered%s" % (
        len(rows), n_dom, len({r[0] for r in rows}) - n_dom, len(failed),
        (" (e.g. %s; rerun to retry)" % failed[0]) if failed else ""), flush=True)


if __name__ == "__main__":
    main()
