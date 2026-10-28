#!/usr/bin/env python3
"""Make test_inputs/data/ from the full test-set inputs, cut to what build.py reads.

    python cut_test_data.py <inputs folder from make_test_inputs.sbatch> data/

The inputs folder (PROVENANCE.md says where it is) is 478 MB, mostly because hmmscan.domtblout
holds the Pfam domains of all 25,729 Swiss-Prot proteins the searches matched, at every
E-value, and jackhmmer.out holds all three rounds. This keeps:
  hmmscan.domtblout   every row of your 51 proteins; for a Swiss-Prot protein only the domains
                      with i-Evalue at most 0.001, the only ones build.py uses (TARGET_PFAM_E)
  jackhmmer.out       the last round of each query, the only one build.py reads
  swissprot.fasta     the Swiss-Prot proteins named in a search table
and copies every other file as it is. Files over 20 MB after cutting are written gzipped;
build.py reads name.gz when name is missing.
"""
import gzip
import shutil
import sys
from pathlib import Path

TARGET_PFAM_E = 1e-3
src, dst = Path(sys.argv[1]), Path(sys.argv[2])
dst.mkdir(parents=True, exist_ok=True)
proteins = {l[1:].split()[0] for l in open(src / "proteins.fasta") if l.startswith(">")}


def write(name, lines):
    tmp = dst / name
    with open(tmp, "w") as fh:
        fh.writelines(lines)
    if tmp.stat().st_size > 20e6:
        with open(tmp, "rb") as a, gzip.open(str(tmp) + ".gz", "wb", compresslevel=9) as b:
            shutil.copyfileobj(a, b)
        tmp.unlink()
        name += ".gz"
    print("%-26s %6.1f MB" % (name, (dst / name).stat().st_size / 1e6))


def hmmscan_rows():
    for line in open(src / "hmmscan.domtblout"):
        if line.startswith("#"):
            yield line
            continue
        f = line.split(None, 13)
        if f[3] in proteins or float(f[12]) <= TARGET_PFAM_E:
            yield line


def jackhmmer_last_round():
    block = []
    for line in open(src / "jackhmmer.out"):
        if line.startswith("Query:"):
            block = [line]
        elif line.startswith("@@ Round:"):
            block = block[:1] + [line]
        elif line.startswith("//"):
            yield from block + [line]
            block = []
        elif block:
            block.append(line)


named = set()
for name, col in (("mmseqs2_3rounds.m8", 1), ("phmmer.domtblout", 0), ("jackhmmer.domtblout", 0), ("kmerseek.tsv", None)):
    for line in open(src / name):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t") if col is None or name.endswith(".m8") else line.split()
        if col is None:
            if f[0] in ("query", "query_name"):
                head = f
                continue
            t = dict(zip(head, f)).get("target") or dict(zip(head, f)).get("target_name")
        else:
            t = f[col]
        named.add(t.split("|")[1] if t.startswith(("sp|", "tr|")) else t)


def swissprot():
    keep = False
    for line in open(src / "swissprot.fasta"):
        if line.startswith(">"):
            h = line[1:].split()[0]
            keep = (h.split("|")[1] if h.startswith(("sp|", "tr|")) else h) in named
        if keep:
            yield line


write("hmmscan.domtblout", hmmscan_rows())
write("jackhmmer.out", jackhmmer_last_round())
write("swissprot.fasta", swissprot())
for name in ("proteins.fasta", "mmseqs2_3rounds.m8", "phmmer.domtblout", "jackhmmer.domtblout", "kmerseek.tsv",
             "phmmer.out", "hmmscan.out", "protein_names.tsv", "run_settings.tsv", "kmerseek_alphabets.rs"):
    write(name, open(src / name))
