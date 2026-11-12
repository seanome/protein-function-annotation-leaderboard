#!/usr/bin/env python3
"""Cut a 50-protein test set for build.py out of real Botryllus schlosseri tool outputs on Sherlock.

No tool is run here except at the end: hmmscan against Pfam, on the 50 proteins and on the
Swiss-Prot proteins they match (see make_test_inputs.sbatch). Every other file is a filtered
copy of an output that already existed. Where a file needed a change to fit build.py's input
format, the change is listed in PROVENANCE.md next to this script.

Proteins and their ids come from the new Botryllus proteome (FUN ids, FUN000003_FUN000003).
The 50 are drawn from the 7,216 FUN proteins that were searched with phmmer, jackhmmer,
MMseqs2 three rounds and kmerseek against Swiss-Prot minus Ascidiacea (2024-kmerseek-analysis
notebook 281 and the dark-set run). No one-round MMseqs2, Foldseek, ProstT5 + Foldseek,
hhsearch or InterProScan output exists for these proteins, so those files are not written and
build.py draws their tracks as "not run".

Python 3.6, standard library only (Sherlock's /usr/bin/python3).
"""
import glob
import gzip
import os
import random
import sys

S = os.environ["SCRATCH"]
FUN_CHUNKS = sorted(glob.glob(S + "/botryllus-mhc/qpos/chunks/chunk_*.fasta"))
QPOS = S + "/botryllus-mhc/qpos/out"
KMERSEEK = S + "/invertebrate-dark-set/results-0.4/dark-regions/botryllus_dark_regions.csv.gz"
KMERSEEK_MOLTYPE, KMERSEEK_K = "hp_lehninger2", "17"
SPROT = S + "/botryllus-mhc-swissprot/db/swissprot/uniprot_sprot.fasta.gz"
REFERENCE = S + "/invertebrate-dark-set/reference/minus_Ascidiacea/reference.fasta"
E = 1e-3  # used only to draw the 50 (none / one / two or more programs with a hit)
E_TARGETS = 10.0  # targets of hits at or below this get hmmscan, so build.py can name their Pfam family
N_NONE, N_ONE, N_TOTAL = 16, 16, 50
# Added to the 50 drawn at random, after the draw, so the draw itself does not change.
# FUN008084_FUN008084 is the Botryllus histocompatibility factor (BHF): its 252 residues equal
# UniProt S5RI84, "Histocompatibility factor" (Voskoboynik et al. 2013, PubMed 23888037).
ALWAYS = ["FUN008084_FUN008084"]
OUT = sys.argv[1]


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def lines(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            yield line


fun = {}
for p in FUN_CHUNKS:
    for h, s in read_fasta(p).items():
        fun[h.split()[0]] = s.rstrip("*")
both = set(fun)
print("FUN proteins", len(fun), flush=True)

# Per tool, rows for the shared proteins, and programs with a hit at E <= 1e-3.
fam_hits = {f: set() for f in both}
rows = {"mmseqs2_3rounds": [], "phmmer": [], "jackhmmer": [], "kmerseek": []}
header = {}

import csv
kh = next(lines(KMERSEEK)).rstrip("\n").split(",")
header["kmerseek"] = kh
iq, ik, im, ie = kh.index("query_name"), kh.index("ksize"), kh.index("moltype"), kh.index("region_evalue")


def collect(keep):
    """With keep=None, record which programs hit each protein at E <= 1e-3 and store no rows
    (the jackhmmer tables alone are several GB). Otherwise store the rows of proteins in keep."""
    for p in sorted(glob.glob(QPOS + "/mmseqs2.chunk_*.tsv.gz")):
        for line in lines(p):
            f = line.split("\t")
            if keep is None:
                if f[0] in fam_hits and float(f[7]) <= E:
                    fam_hits[f[0]].add("mmseqs2")
            elif f[0] in keep:
                rows["mmseqs2_3rounds"].append(f)
    for tool in ("phmmer", "jackhmmer"):
        for p in sorted(glob.glob(QPOS + "/" + tool + ".chunk_*.domtblout.gz")):
            for line in lines(p):
                if line.startswith("#"):
                    continue
                f = line.split()
                if keep is None:
                    if f[3] in fam_hits and float(f[12]) <= E:
                        fam_hits[f[3]].add(tool)
                elif f[3] in keep:
                    rows[tool].append(line)
    it = lines(KMERSEEK)
    next(it)
    for f in csv.reader(it):
        if f[im] != KMERSEEK_MOLTYPE or f[ik] != KMERSEEK_K:
            continue
        if keep is None:
            if f[iq] in fam_hits and float(f[ie]) <= E:
                fam_hits[f[iq]].add("kmerseek")
        elif f[iq] in keep:
            rows["kmerseek"].append(f)


collect(None)

rng = random.Random(20261027)
by_n = {}
for p in sorted(both):
    by_n.setdefault(len(fam_hits[p]), []).append(p)
print("proteins by number of programs with a hit:", {k: len(v) for k, v in sorted(by_n.items())}, flush=True)
none = rng.sample(by_n.get(0, []), min(N_NONE, len(by_n.get(0, []))))
# exactly one program: spread over which program it is
one_pool = by_n.get(1, [])
one_by_fam = {}
for p in one_pool:
    one_by_fam.setdefault(next(iter(fam_hits[p])), []).append(p)
one = []
while len(one) < N_ONE and any(one_by_fam.values()):
    for fam in sorted(one_by_fam):
        if one_by_fam[fam] and len(one) < N_ONE:
            one.append(one_by_fam[fam].pop(rng.randrange(len(one_by_fam[fam]))))
many_pool = [p for n, ps in by_n.items() if n >= 2 for p in ps]
many = rng.sample(many_pool, N_TOTAL - len(none) - len(one))
chosen = sorted(set(none + one + many + [p for p in ALWAYS if p in fun]))
missing = [p for p in ALWAYS if p not in fun]
if missing:
    sys.exit("not in the searched proteins: " + ", ".join(missing))
cs = set(chosen)
print("chosen", len(chosen), "none", len(none), "one", len(one), "two or more", len(many),
      "always", ALWAYS, flush=True)
collect(cs)
print({k: len(v) for k, v in rows.items()}, flush=True)

os.makedirs(OUT, exist_ok=True)
with open(OUT + "/proteins.fasta", "w") as fh:
    for p in chosen:
        s = fun[p]
        fh.write(">" + p + "\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")

SPEC = "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,theader".split(",")

h = "query,target,qstart,qend,tstart,tend,bits,evalue,qlen,tlen".split(",")
with open(OUT + "/mmseqs2_3rounds.m8", "w") as fh:
    for f in rows["mmseqs2_3rounds"]:
        if f[0] not in cs:
            continue
        d = dict(zip(h, [x.strip() for x in f]))
        fh.write("\t".join(d.get(c, "") for c in SPEC) + "\n")
DOM_HEAD = (
    "#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord\n"
    "# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target\n"
)
for tool in ("phmmer", "jackhmmer"):
    with open(OUT + "/" + tool + ".domtblout", "w") as fh:
        fh.write(DOM_HEAD)
        for line in rows[tool]:
            if line.split()[3] in cs:
                fh.write(line)
with open(OUT + "/kmerseek.tsv", "w") as fh:
    fh.write("\t".join(header["kmerseek"]) + "\n")
    for f in rows["kmerseek"]:
        if f[iq] in cs:
            fh.write("\t".join(x.replace("\t", " ") for x in f) + "\n")

# Swiss-Prot sequences for every target any file names; and the targets of hits at E <= 10,
# which the sbatch runs hmmscan on so the agreement rule can name a target's Pfam family.
targets, strong = set(), set()


def acc(t):
    if t.startswith("AF-"):
        return t.split("-")[1]
    if t.startswith("sp|") or t.startswith("tr|"):
        return t.split("|")[1]
    return t


for p in ("mmseqs2_3rounds.m8",):
    for line in open(OUT + "/" + p):
        f = line.rstrip("\n").split("\t")
        targets.add(acc(f[1]))
        if float(f[8]) <= E_TARGETS:
            strong.add(acc(f[1]))
for p in ("phmmer.domtblout", "jackhmmer.domtblout"):
    for line in open(OUT + "/" + p):
        if line.startswith("#"):
            continue
        f = line.split()
        targets.add(acc(f[0]))
        if float(f[12]) <= E_TARGETS:
            strong.add(acc(f[0]))
it = csv.reader(open(OUT + "/kmerseek.tsv"), delimiter="\t")
next(it)
it_t = kh.index("target_name")
for f in it:
    targets.add(acc(f[it_t]))
    if float(f[ie]) <= E_TARGETS:
        strong.add(acc(f[it_t]))

found = {}
for hdr, s in read_fasta(SPROT).items():
    a = hdr.split("|")[1]
    if a in targets:
        found[a] = (hdr, s)  # whole header: the description, OS= and GN= are shown in the pages
missing = targets - set(found)
ref = read_fasta(REFERENCE)
from_ref = 0
for hdr, s in ref.items():
    a = hdr.split()[0]
    if a in missing:
        found[a] = (" ".join(["sp|" + a + "|"] + hdr.split()[1:]), s)
        from_ref += 1
print("targets", len(targets), "in Swiss-Prot 2026_03", len(targets) - len(missing),
      "taken from the dark-set reference", from_ref, "not found", len(targets - set(found)), flush=True)
with open(OUT + "/swissprot.fasta", "w") as fh:
    for a in sorted(found):
        hdr, s = found[a]
        fh.write(">" + hdr + "\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
with open(OUT + "/../hmmscan_targets.fasta", "w") as fh:
    for a in sorted(strong & set(found)):
        hdr, s = found[a]
        fh.write(">" + a + "\n" + s + "\n")
print("hmmscan targets", len(strong & set(found)), flush=True)

with open(OUT + "/../selection.tsv", "w") as fh:
    fh.write("protein\tgroup_at_selection\tfamilies_with_hit_at_selection\n")
    for p in chosen:
        n = len(fam_hits.get(p, ()))
        grp = "none" if n == 0 else ("one" if n == 1 else "two_or_more")
        if p in ALWAYS:
            grp = "always_" + grp
        fh.write("%s\t%s\t%s\n" % (p, grp, ",".join(sorted(fam_hits.get(p, ())))))
