"""Write unifire_predictions.tsv: UniFIRE's predictions for the Botryllus proteome, one row per
prediction, under the protein IDs the Biomni tables use.

    python make_unifire_table.py proteome_lengths.tsv.gz <unifire_whole_proteome folder> | gzip > unifire_predictions.tsv.gz

UniFIRE 2026.3 (ARBA, UniRule and PIRSR rules) was run on InterProScan 5.59 output by
olgabot/botryllus-mhc PR 13 (branch olgabot/annotation-timing-pilot, commit 257832d, folder
nextflow_runs/timing-pilot/unifire_whole_proteome). Proteins over 5,000 amino acids were not run.

The two sets of IDs number repeated FASTA headers differently. UniFIRE keeps the whole header,
with "/" changed to "_", and adds .1, .2, ... to every copy of a header that appears more than
once (FUN000051_GALNT6.1, FUN000051_GALNT6.2). The Biomni tables keep the part before the first
underscore and add .2, .3 to the later copies only (FUN000051, FUN000051.2); see
make_proteome_lengths.py. Each such part belongs to one header in the FASTA (checked on
2026-10-31), so the copy number alone links the two.
"""
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path

SOURCES = [("arba", "ARBA"), ("unirule", "UniRule"), ("unirule-pirsr", "PIRSR")]


def id_mapper(lengths):
    """A function from a whole-header protein name (FUN000051_GALNT6.2) to the Biomni ID
    (FUN000051.2), checked against proteome_lengths.tsv."""
    ids = [r["protein_id"] for r in csv.DictReader(gzip.open(lengths, "rt"), delimiter="\t")]
    copies = Counter(i.split(".")[0] for i in ids)
    known = set(ids)

    def biomni_id(uid):
        gene = uid.split("_")[0]
        n = int(re.search(r"\.(\d+)$", uid).group(1)) if copies[gene] > 1 else 1
        pid = gene if n == 1 else "%s.%d" % (gene, n)
        if pid not in known:
            sys.exit("protein %s has no Biomni ID (%s)" % (uid, pid))
        return pid
    return biomni_id


def main(lengths, folder):
    biomni_id = id_mapper(lengths)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["protein_id", "source", "rule", "type", "value", "start", "end"])
    for f, name in SOURCES:
        for r in csv.DictReader(gzip.open(Path(folder) / ("predictions_%s.out.gz" % f), "rt"), delimiter="\t"):
            w.writerow([biomni_id(r["ProteinId"]), name, r["Evidence"], r["AnnotationType"], r["Value"], r["Start"], r["End"]])


if __name__ == "__main__":
    main(*sys.argv[1:])
