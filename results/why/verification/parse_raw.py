"""Parse the big raw files once into small parquet caches in this folder.

Reads only from /Users/olga/data/why-numbers/2026_03/ and writes only here.
"""
import gzip
import io
import sys
import tarfile
from pathlib import Path

import polars as pl

R = Path("/Users/olga/data/why-numbers/2026_03")
OUT = Path(__file__).resolve().parent / "cache"
OUT.mkdir(exist_ok=True)


def parse_sprot():
    acc, ox, pdb, n_ox_lines, os_names = [], [], [], [], []
    cur_acc = None
    cur_ox = []
    cur_pdb = 0
    cur_os = []
    with gzip.open(R / "uniprot_sprot.dat.gz", "rt") as fh:
        for line in fh:
            tag = line[:2]
            if tag == "AC":
                if cur_acc is None:
                    cur_acc = line[5:].split(";")[0].strip()
            elif tag == "OX":
                # OX   NCBI_TaxID=9606;  (may carry {ECO...})
                v = line[5:].strip()
                v = v.split("NCBI_TaxID=")[1]
                v = v.split(";")[0].split("{")[0].strip()
                cur_ox.append(int(v))
            elif tag == "OS":
                cur_os.append(line[5:].strip())
            elif tag == "DR":
                if line[5:9] == "PDB;":
                    cur_pdb += 1
            elif tag == "//":
                acc.append(cur_acc)
                ox.append(cur_ox[0] if cur_ox else None)
                n_ox_lines.append(len(cur_ox))
                pdb.append(cur_pdb)
                os_names.append(" ".join(cur_os))
                cur_acc, cur_ox, cur_pdb, cur_os = None, [], 0, []
    df = pl.DataFrame(
        {"acc": acc, "ox": ox, "n_ox": n_ox_lines, "n_pdb": pdb, "os": os_names},
        schema={"acc": pl.String, "ox": pl.Int64, "n_ox": pl.Int32, "n_pdb": pl.Int32, "os": pl.String},
    )
    df.write_parquet(OUT / "sprot.parquet")
    print("sprot", df.shape)


def parse_taxdump():
    with tarfile.open(R / "taxdump.tar.gz") as tf:
        def read(name):
            return tf.extractfile(name).read().decode()

        nodes = read("nodes.dmp")
        names = read("names.dmp")
        merged = read("merged.dmp")
        deleted = read("delnodes.dmp")
    tid, par, rank = [], [], []
    for ln in nodes.splitlines():
        f = ln.split("\t|\t")
        tid.append(int(f[0]))
        par.append(int(f[1]))
        rank.append(f[2])
    pl.DataFrame({"taxid": tid, "parent": par, "rank": rank}).write_parquet(OUT / "nodes.parquet")
    nt, nn = [], []
    for ln in names.splitlines():
        f = ln.rstrip("\t|").split("\t|\t")
        if f[3].strip() == "scientific name":
            nt.append(int(f[0]))
            nn.append(f[1])
    pl.DataFrame({"taxid": nt, "name": nn}).write_parquet(OUT / "names.parquet")
    mo, mn = [], []
    for ln in merged.splitlines():
        f = ln.rstrip("\t|").split("\t|\t")
        mo.append(int(f[0]))
        mn.append(int(f[1].strip()))
    pl.DataFrame({"old": mo, "new": mn}).write_parquet(OUT / "merged.parquet")
    dl = [int(ln.split("\t")[0]) for ln in deleted.splitlines() if ln.strip()]
    pl.DataFrame({"taxid": dl}).write_parquet(OUT / "deleted.parquet")
    print("tax", len(tid), len(nt), len(mo), len(dl))


if __name__ == "__main__":
    which = sys.argv[1:] or ["tax", "sprot"]
    if "tax" in which:
        parse_taxdump()
    if "sprot" in which:
        parse_sprot()
