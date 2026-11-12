"""Download the inputs on Sherlock, all from one UniProt release.

Python 3.6 standard library (Sherlock's /usr/bin/python3). Two modes:

  --setup (setup.sbatch, once): relnotes.txt, which fixes the release for the whole run, and
      Swiss-Prot (uniprot_sprot.fasta.gz) of that release.
  --row N (fetch.sbatch, one array task per proteome; N is 0-based, header not counted):
      <upid>_<taxid>.fasta.gz, the one-protein-per-gene FASTA from the release's reference
      proteomes FTP folder (the proteins counted), and
      <upid>_uniprot.tsv.gz, UniProt fields for every entry of the proteome: accession,
      protein name, and the InterPro cross-references labels.py used (INTERPRO_XREFS below,
      copied from scripts/why/fetch.py). Paged /search, 500 rows a page, as fetch.py does,
      because /stream dropped the connection on large downloads.

Every UniProt page must report the release written by --setup, and the row count must equal
X-Total-Results, or the task stops. A file already in its manifest with a matching sha256 is
not downloaded again, so a rerun only fetches what is missing.

    python3 -I fetch_inputs.py --setup --out $RUN/inputs
    python3 -I fetch_inputs.py --row 3 --proteomes proteomes.tsv --out $RUN/inputs
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request

UA = "seanome-protein-function-annotation-leaderboard (olga@seanome.org)"
FTP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release"
# Copied from scripts/why/fetch.py, INTERPRO_XREFS.
INTERPRO_XREFS = (
    "xref_interpro,xref_pfam,xref_panther,xref_cdd,xref_gene3d,xref_supfam,xref_smart,"
    "xref_prosite,xref_prints,xref_pirsf,xref_hamap,xref_sfld,xref_ncbifam"
)
FIELDS = "accession,protein_name," + INTERPRO_XREFS


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def curl(url, dest):
    tmp = dest + ".part"
    subprocess.check_call(["curl", "-sS", "--fail", "-L", "--retry", "5", "--retry-delay", "10", "-A", UA, "-o", tmp, url])
    os.replace(tmp, dest)


def get(url, tries=8):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=300) as r:
                return {k.lower(): v for k, v in r.headers.items()}, r.read().decode()
        except Exception as e:
            if i == tries - 1:
                raise
            print("retry", i + 1, url, e, file=sys.stderr, flush=True)
            time.sleep(10 * (i + 1))


def uniprot_paged(query, dest, release):
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": query, "fields": FIELDS, "format": "tsv", "size": 500}
    )
    first, n, total, seen = url, 0, None, set()
    tmp = dest + ".part"
    with gzip.open(tmp, "wt") as out:
        while url:
            headers, body = get(url)
            seen.add(headers.get("x-uniprot-release"))
            total = int(headers["x-total-results"])
            lines = body.splitlines()
            if n == 0:
                out.write(lines[0] + "\n")
            for line in lines[1:]:
                out.write(line + "\n")
                n += 1
            m = re.search(r'<([^>]+)>; rel="next"', headers.get("link", ""))
            url = m.group(1) if m else None
    if seen != {release} or n != total:
        raise RuntimeError("{}: releases {}, expected {}; rows {}, X-Total-Results {}".format(query, seen, release, n, total))
    os.replace(tmp, dest)
    return first, n


def count_fasta(path):
    with gzip.open(path, "rt") as fh:
        return sum(1 for line in fh if line.startswith(">"))


class Manifest:
    def __init__(self, path):
        self.path = path
        self.d = json.load(open(path)) if os.path.exists(path) else {}

    def done(self, name, path):
        return name in self.d and os.path.exists(path) and self.d[name]["sha256"] == sha256(path)

    def record(self, name, path, url, **extra):
        self.d[name] = dict(file=os.path.basename(path), url=url, sha256=sha256(path), bytes=os.path.getsize(path), fetched_utc=now(), **extra)
        with open(self.path, "w") as fh:
            json.dump(self.d, fh, indent=1, sort_keys=True)


def setup(out):
    man = Manifest(os.path.join(out, "setup_manifest.json"))
    rel_path = os.path.join(out, "relnotes.txt")
    if not man.done("relnotes", rel_path):
        curl(FTP + "/relnotes.txt", rel_path)
        # First line: "UniProt Release 2026_03".
        release = re.search(r"Release (\d{4}_\d{2})", open(rel_path).read()).group(1)
        man.record("relnotes", rel_path, FTP + "/relnotes.txt", release=release)
    release = man.d["relnotes"]["release"]
    sprot = os.path.join(out, "uniprot_sprot.fasta.gz")
    if not man.done("swissprot_fasta", sprot):
        url = FTP + "/knowledgebase/complete/uniprot_sprot.fasta.gz"
        curl(url, sprot)
        # FTP files do not carry the release; check relnotes again after the download.
        _, notes = get(FTP + "/relnotes.txt")
        if "Release " + release not in notes:
            raise RuntimeError("UniProt release changed during setup; rerun setup")
        man.record("swissprot_fasta", sprot, url, release=release, proteins=count_fasta(sprot))
    print("UniProt release", release, "; Swiss-Prot", man.d["swissprot_fasta"]["proteins"], "proteins", flush=True)


def one(out, proteomes, row):
    release = json.load(open(os.path.join(out, "setup_manifest.json")))["relnotes"]["release"]
    with open(proteomes) as fh:
        r = list(csv.DictReader(fh, delimiter="\t"))[row]
    upid, taxid = r["upid"], r["proteome_taxid"]
    man = Manifest(os.path.join(out, upid + "_manifest.json"))
    fa = os.path.join(out, "{}_{}.fasta.gz".format(upid, taxid))
    if not man.done("fasta", fa):
        curl(r["fasta_url"], fa)
        _, notes = get(FTP + "/relnotes.txt")
        if "Release " + release not in notes:
            raise RuntimeError("UniProt release is no longer " + release)
        man.record("fasta", fa, r["fasta_url"], release=release, proteins=count_fasta(fa))
    tsv = os.path.join(out, upid + "_uniprot.tsv.gz")
    if not man.done("uniprot", tsv):
        url, n = uniprot_paged("proteome:" + upid, tsv, release)
        man.record("uniprot", tsv, url, release=release, rows=n)
    print(upid, r["species"], man.d["fasta"]["proteins"], "proteins;", man.d["uniprot"]["rows"], "UniProt rows", flush=True)


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--setup", action="store_true")
    g.add_argument("--row", type=int)
    ap.add_argument("--proteomes")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    if args.setup:
        setup(args.out)
    else:
        one(args.out, args.proteomes, args.row)


if __name__ == "__main__":
    main()
