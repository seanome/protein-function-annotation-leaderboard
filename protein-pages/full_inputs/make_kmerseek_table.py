"""Write kmerseek_best_per_protein.tsv: for each Botryllus protein, its best kmerseek hit against
Swiss-Prot, under the protein IDs the Biomni tables use.

    python make_kmerseek_table.py proteome_lengths.tsv.gz kmerseek_botryllus_swissprot.per_protein.tsv.gz | gzip > kmerseek_best_per_protein.tsv.gz

The input is results/summary/kmerseek_botryllus_swissprot.per_protein.tsv.gz of the kmerseek run
on Sherlock in $SCRATCH/botryllus-kmerseek-swissprot: every protein of the proteome searched
against all of UniProt Swiss-Prot 2026_03 at 38 alphabet and k-mer size settings. The best hit is
the lowest region E-value over the 35 settings whose E-value fit worked. That run named proteins
by their whole FASTA header, with ".1", ".2" on repeated headers and commas changed to "-"
(query/query_renaming.tsv); make_unifire_table.id_mapper turns those names into Biomni IDs.
"""
import csv
import gzip
import sys

from make_unifire_table import id_mapper


def main(lengths, per_protein):
    biomni_id = id_mapper(lengths)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["protein_id", "n_settings", "n_settings_strong", "alphabet", "ksize", "evalue", "target",
                "target_entry", "target_name", "target_organism"])
    seen = set()
    for r in csv.DictReader(gzip.open(per_protein, "rt"), delimiter="\t"):
        pid = biomni_id(r["query_name"])
        if pid in seen:
            sys.exit("%s twice" % pid)
        seen.add(pid)
        w.writerow([pid, r["n_fitted_pairs"], r["n_pairs_with_strong_match"], r["best_alphabet"], r["best_ksize"],
                    r["best_evalue"], r["best_target"], r["best_target_entry"], r["best_target_protein_name"],
                    r["best_target_organism"]])


if __name__ == "__main__":
    main(*sys.argv[1:])
