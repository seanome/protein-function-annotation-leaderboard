"""How often FANTASIA's nearest protein has the same gene as DIAMOND's best Swiss-Prot hit, by cosine distance.

    python fantasia_vs_diamond.py . > fantasia_vs_diamond.tsv

FANTASIA was run on the 26,676 proteins with no Swiss-Prot match or in the immune table, so the
proteins with both a FANTASIA answer and a DIAMOND hit are mostly immune-table proteins. Of those,
only pairs where UniProt gives a gene name for both accessions are counted. "Same gene" is the same
gene name, ignoring case; "same gene family" also drops a trailing number and letter (NLRC5 and
NLRC3 are both NLRC). FANTASIA's reliability index is 1 minus the cosine distance in every row, so
it adds nothing and is not used.
"""
import csv
import gzip
import re
import sys
from pathlib import Path

BINS = [(0, 0.1), (0.1, 0.15), (0.15, 0.2), (0.2, 0.25), (0.25, 0.3), (0.3, 1.0)]

inp = Path(sys.argv[1])
gene = {}
for line in gzip.open(inp / "uniprot_names.tsv.gz", "rt"):
    acc, _name, g, _org = line.rstrip("\n").split("\t")
    if g:
        gene[acc] = g.upper()
family = lambda g: re.sub(r"\d+[A-Z]?$", "", g)

pairs = []
for r in csv.DictReader(gzip.open(inp / "botryllus_protein_annotations_v2.tsv.gz", "rt"), delimiter="\t"):
    if r["fantasia_distance"]:
        assert abs(float(r["fantasia_ri"]) - (1 - float(r["fantasia_distance"]))) < 1e-6
        if r["sp_acc"] in gene and r["fantasia_uniprot"] in gene:
            pairs.append((float(r["fantasia_distance"]), gene[r["sp_acc"]], gene[r["fantasia_uniprot"]]))

print("distance_from\tdistance_to\tproteins\tsame_gene\tsame_gene_family\tsame_gene_pct\tsame_family_pct")
for lo, hi in BINS:
    s = [p for p in pairs if lo <= p[0] < hi]
    g = sum(a == b for _, a, b in s)
    f = sum(family(a) == family(b) for _, a, b in s)
    print(f"{lo}\t{hi}\t{len(s)}\t{g}\t{f}\t{100 * g / len(s):.0f}\t{100 * f / len(s):.0f}")
