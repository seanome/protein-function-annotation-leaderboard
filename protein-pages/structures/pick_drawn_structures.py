"""Copy the predicted structures of the drawn test proteins out of Biomni's zip of structures.

    python pick_drawn_structures.py botryllus_structures_le200aa.zip ../full_inputs/proteome_lengths.tsv.gz \
        ../test_inputs/data/proteins.fasta ../test_inputs/data/structures

The zip (from Biomni, given by Olga on 2026-11-03) holds one PDB file per protein of 200 aa or
less, 17,990 in all, named by the Biomni protein ID (FUN015821.2.pdb). The drawn proteins are
named differently, so they are matched by the MD5 of their sequence and written as
<drawn ID>.pdb. A drawn protein that already has a structure (BHF's, from ESMFold2) keeps it.
Each structure's residues are checked against the drawn sequence before it is written, and a
first REMARK line says where it came from (build.py shows it on the page).
"""
import csv
import gzip
import hashlib
import sys
import zipfile
from pathlib import Path

# The first line, which the pages show as how the structure was made. Biomni's files do not name
# the program; Biomni said on 2026-11-03 that it was ESMFold2's fast model (esmcfold2-fast, the
# folding head on the ESMC 6B language model, `esmc fold --model fast`, default settings), one
# structure per protein, no templates or MSA. Its pLDDT is on a 0 to 1 scale; build.py rescales it.
REMARK = ("REMARK   1 ESMFold2, fast model (esmcfold2-fast: folding head on ESMC 6B, esmc fold --model fast, "
          "default settings, no templates or MSA), run by Biomni for proteins of 200 aa or less "
          "(botryllus_structures_le200aa.zip, 2026-11-03)")
THREE = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I",
             LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")


def fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        else:
            out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def main(zip_path, lengths, drawn_fasta, out_dir):
    out_dir = Path(out_dir)
    by_md5 = {}
    for r in sorted(csv.DictReader(gzip.open(lengths, "rt"), delimiter="\t"), key=lambda r: r["protein_id"]):
        by_md5.setdefault(r["md5"], r["protein_id"])
    z = zipfile.ZipFile(zip_path)
    have = {Path(n).stem: n for n in z.namelist() if n.endswith(".pdb")}
    n = 0
    for did, seq in fasta(drawn_fasta).items():
        bid = by_md5.get(hashlib.md5(seq.encode()).hexdigest())
        if bid not in have:
            continue
        if any((out_dir / (x + ".pdb")).exists() for x in (did, did.split("_")[0])):
            print("%s already has a structure; kept" % did, file=sys.stderr)
            continue
        text = z.read(have[bid]).decode()
        res = "".join(THREE.get(l[17:20], "X") for l in text.splitlines() if l.startswith("ATOM") and l[12:16].strip() == "CA")
        if res != "".join(c if c in THREE.values() else "X" for c in seq):
            sys.exit("%s.pdb does not have the residues of %s" % (bid, did))
        (out_dir / (did + ".pdb")).write_text(REMARK + "\n" + text)
        n += 1
    print("%d structures written" % n, file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:])
