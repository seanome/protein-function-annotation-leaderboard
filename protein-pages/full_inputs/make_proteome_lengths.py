"""Write proteome_lengths.tsv: for each protein of botryllus_new_genome_protein.fasta, the ID the
Biomni tables use, its length and the MD5 of its sequence. The sequences themselves are not
committed (the proteome is unpublished); the MD5 is enough to match a protein in another FASTA.

    python make_proteome_lengths.py botryllus_new_genome_protein.fasta > proteome_lengths.tsv

The FASTA repeats a header for each isoform of a gene (three ">FUN015821_FUN015821"). The Biomni
tables name the first copy by the part before the underscore (FUN015821) and later copies
FUN015821.2, FUN015821.3, in file order. tests/test_build_full.py checks this naming against the
protein lengths DIAMOND printed in botryllus_allorecognition_hits.tsv.
"""
import hashlib
import sys
from collections import Counter

seen = Counter()
recs = []
for line in open(sys.argv[1]):
    line = line.strip()
    if line.startswith(">"):
        gene = line[1:].split()[0].split("_")[0]
        seen[gene] += 1
        recs.append([gene if seen[gene] == 1 else "%s.%d" % (gene, seen[gene]), []])
    elif line:
        recs[-1][1].append(line)
print("protein_id\tlength\tmd5")
for pid, parts in recs:
    s = "".join(parts).rstrip("*")
    print("%s\t%d\t%s" % (pid, len(s), hashlib.md5(s.encode()).hexdigest()))
