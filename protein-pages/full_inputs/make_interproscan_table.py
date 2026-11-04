"""Write interproscan.tsv: the whole-proteome InterProScan matches under the protein IDs the Biomni
tables use, in InterProScan's TSV layout (protein, MD5, length, member database, signature,
description, start, end, score, status, date, InterPro entry, InterPro description).

    python make_interproscan_table.py proteome_lengths.tsv.gz interproscan_whole_proteome.tsv.gz | gzip > interproscan.tsv.gz

The input is ../test_inputs/interproscan_xml_to_tsv.py run on the XML of the InterProScan 5.59
run (nf-core/proteinannotator, olgabot/botryllus-mhc PR 13; on Sherlock in
$SCRATCH/botryllus-annotation-full/results/proteinannotator/interproscan). UniFIRE's predictions
come from the same run. Proteins over 5,000 amino acids were not run. Each row's MD5 is checked
against the protein's MD5 in proteome_lengths.tsv.
"""
import csv
import gzip
import sys

from make_unifire_table import id_mapper


def main(lengths, tsv):
    biomni_id = id_mapper(lengths)
    md5 = {r["protein_id"]: r["md5"] for r in csv.DictReader(gzip.open(lengths, "rt"), delimiter="\t")}
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    for f in csv.reader(gzip.open(tsv, "rt"), delimiter="\t"):
        pid = biomni_id(f[0])
        if md5[pid] != f[1]:
            sys.exit("%s (%s): sequence MD5 differs from proteome_lengths.tsv" % (f[0], pid))
        w.writerow([pid] + f[1:])


if __name__ == "__main__":
    main(*sys.argv[1:])
