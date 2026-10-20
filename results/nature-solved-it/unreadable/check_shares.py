#!/usr/bin/env python3
"""Second-agent check: recount every proteome's unlabelled share from the raw tables.

Written independently of scripts/nature_solved_it/unreadable/count_unreadable.py (not read).
Definition, from PLAN.md: a protein is unlabelled if it has neither an InterPro match (any
non-empty InterPro or member-database column in the UniProt table) nor an MMseqs2 hit to
Swiss-Prot with fident > 0.30 (a hit to its own entry counts). The proteins are the
accessions in the one-protein-per-gene FASTA.

Usage (Python 3.6, standard library only):
    python3 check_shares.py <unreadable dir> > check_raw.tsv
"""
import glob
import gzip
import json
import os
import sys

MEMBER_COLUMNS = ["InterPro", "Pfam", "PANTHER", "CDD", "Gene3D", "SUPFAM", "SMART",
                  "PROSITE", "PRINTS", "PIRSF", "HAMAP", "SFLD", "NCBIfam"]


def fasta_accessions(path):
    accs = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                # >sp|P12345|NAME_SPECIES ... or >tr|...
                accs.append(line[1:].split("|")[1])
    return accs


def main(base):
    out_dirs = sorted(glob.glob(os.path.join(base, "out", "UP*")))
    cols = ["upid", "n_fasta_records", "n_proteins", "n_interpro", "n_sprot_hit_gt30",
            "n_unlabelled", "unlabelled_pct", "n_uncharacterized_name",
            "n_query_ids_not_in_fasta", "n_fasta_missing_from_uniprot_table",
            "summary_n_proteins", "summary_n_unlabelled", "summary_unlabelled_pct"]
    print("\t".join(cols))
    for d in out_dirs:
        upid = os.path.basename(d)
        fa = glob.glob(os.path.join(base, "inputs", upid + "_*.fasta.gz"))
        assert len(fa) == 1, (upid, fa)
        accs = fasta_accessions(fa[0])
        acc_set = set(accs)

        interpro = set()
        uncharacterized = set()
        seen_in_table = set()
        with gzip.open(os.path.join(base, "inputs", upid + "_uniprot.tsv.gz"), "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = [header.index(c) for c in MEMBER_COLUMNS]
            i_entry = header.index("Entry")
            i_name = header.index("Protein names")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                acc = f[i_entry]
                if acc not in acc_set:
                    continue
                seen_in_table.add(acc)
                if any(f[i].strip().strip(";") for i in idx):
                    interpro.add(acc)
                if f[i_name] == "Uncharacterized protein":
                    uncharacterized.add(acc)

        hit = set()
        stray = set()
        with gzip.open(os.path.join(d, "mmseqs_vs_sprot.tsv.gz"), "rt") as fh:
            for line in fh:
                f = line.split("\t", 3)
                if float(f[2]) > 0.30:
                    q = f[0]
                    if q in acc_set:
                        hit.add(q)
                    else:
                        stray.add(q)

        unlabelled = acc_set - interpro - hit
        n = len(acc_set)
        s = {}
        sp = os.path.join(d, "summary.json")
        if os.path.exists(sp):
            with open(sp) as fh:
                s = json.load(fh)
        row = [upid, len(accs), n, len(interpro), len(hit), len(unlabelled),
               "%.4f" % (100.0 * len(unlabelled) / n), len(uncharacterized), len(stray),
               len(acc_set - seen_in_table), s.get("n_proteins", ""),
               s.get("n_unlabelled", ""), s.get("unlabelled_pct", "")]
        print("\t".join(str(x) for x in row))
        sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1])
