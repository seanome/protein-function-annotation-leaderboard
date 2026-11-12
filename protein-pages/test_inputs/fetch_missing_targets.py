"""Add to swissprot.fasta the Swiss-Prot proteins a track names that the file does not hold yet,
fetched from UniProt REST (current release), so build.py can align them and give their length.

    python test_inputs/fetch_missing_targets.py test_inputs/data prostt5_foldseek.m8.gz

The track file is an m8 table in build.py's column order; its last column is the accession.
"""
import gzip
import sys
import time
import urllib.request
from pathlib import Path

UNIPROT = "https://rest.uniprot.org/uniprotkb/accessions?accessions={}&format=fasta&size=500"


def get(url, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "protein-pages/fetch_missing_targets (olgabot)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode()
        except Exception as e:  # network errors; UniProt sometimes answers 5xx
            if i == tries - 1:
                raise
            print("retry", i + 1, e, file=sys.stderr)
            time.sleep(5 * (i + 1))


def main(folder, track):
    folder = Path(folder)
    fa = folder / "swissprot.fasta"
    have = set()
    for line in open(fa):
        if line.startswith(">"):
            h = line[1:].split()[0]
            have.add(h.split("|")[1] if "|" in h else h)
    path = folder / track
    opener = gzip.open if track.endswith(".gz") else open
    want = sorted({line.rstrip("\n").split("\t")[-1] for line in opener(path, "rt")} - have)
    added = set()
    with open(fa, "a") as out:
        for i in range(0, len(want), 400):
            text = get(UNIPROT.format(",".join(want[i:i + 400])))
            keep = False
            for line in text.splitlines():
                if line.startswith(">"):
                    acc = line[1:].split("|")[1]
                    keep = acc in want and acc not in added
                    if keep:
                        added.add(acc)
                if keep:
                    out.write(line + "\n")
            print(len(added), "of", len(want), file=sys.stderr)
    missing = [a for a in want if a not in added]
    print("added %d; %d not found in UniProt: %s" % (len(added), len(missing), missing[:10]), file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:])
