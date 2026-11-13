"""One proteome's numbers: the unlabelled share, the "Uncharacterized protein" share and BUSCO.

Runs at the end of each array task of unreadable.sbatch. Python 3.6 standard library.

The unlabelled share is the leaderboard's (scripts/why/labels.py, `{tag}_no_interpro_no_sprot_pct`):
"share of proteins with neither an InterPro match nor a Swiss-Prot hit above 30% identity
(a hit to its own entry counts)". Same pieces as labels.py:
  - proteins: accessions in the one-protein-per-gene FASTA;
  - InterPro match: any non-empty InterPro or member-database column in UniProt's
    cross-references (labels.py interpro_xref_set);
  - Swiss-Prot hit: an MMseqs2 row with fident > 0.30 (labels.py MIN_FIDENT), self-hits kept;
  - accessions from "sp|P12345|NAME" headers reduced to "P12345" (labels.py acc_of);
  - share rounded to 2 decimals, as labels.py does.
Like labels.py, it stops if a FASTA accession is missing from the UniProt table.

    python3 -I count_unreadable.py --fasta X.fasta.gz --uniprot X_uniprot.tsv.gz \
        --mmseqs mmseqs_vs_sprot.tsv.gz --busco busco_metazoa_odb12 [--busco busco_eukaryota_odb12] \
        --out summary.json
"""

import argparse
import csv
import glob
import gzip
import json
import os
import re
import sys

MIN_FIDENT = 0.30
NOT_XREF = ("Entry", "Protein names")
UNCHARACTERIZED = "Uncharacterized protein"


def acc_of(s):
    parts = s.split("|")
    return parts[1] if len(parts) > 1 else s


def fasta_accessions(path):
    accs = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                accs.append(acc_of(line[1:].split()[0]))
    return accs


def busco_numbers(run_dir):
    """Completeness from full_table.tsv (one status per BUSCO id), checked against the
    C: percentage BUSCO prints in its short summary."""
    tables = glob.glob(os.path.join(run_dir, "run_*", "full_table.tsv"))
    if len(tables) != 1:
        raise RuntimeError("{}: expected one run_*/full_table.tsv, found {}".format(run_dir, tables))
    status = {}
    with open(tables[0]) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            status.setdefault(f[0], f[1])
    counts = {}
    for s in status.values():
        counts[s] = counts.get(s, 0) + 1
    n = len(status)
    complete = counts.get("Complete", 0) + counts.get("Duplicated", 0)
    pct = 100.0 * complete / n
    lineage = os.path.basename(os.path.dirname(tables[0]))[len("run_"):]
    printed = None
    for txt in glob.glob(os.path.join(run_dir, "short_summary*.txt")) + glob.glob(os.path.join(run_dir, "run_*", "short_summary*.txt")):
        m = re.search(r"C:([\d.]+)%", open(txt).read())
        if m:
            printed = float(m.group(1))
    if printed is not None and abs(printed - pct) > 0.1:
        raise RuntimeError("{}: full_table gives {:.2f}% complete, short summary {}%".format(run_dir, pct, printed))
    return {
        "lineage": lineage,
        "buscos": n,
        "complete_single": counts.get("Complete", 0),
        "complete_duplicated": counts.get("Duplicated", 0),
        "fragmented": counts.get("Fragmented", 0),
        "missing": counts.get("Missing", 0),
        "completeness_pct": round(pct, 2),
        "short_summary_c_pct": printed,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--uniprot", required=True)
    ap.add_argument("--mmseqs", required=True)
    ap.add_argument("--busco", action="append", default=[])
    ap.add_argument("--out", required=True)
    ap.add_argument("--meta", default="{}", help="JSON object copied into the summary (row of proteomes.tsv)")
    args = ap.parse_args()

    accs = fasta_accessions(args.fasta)
    if len(set(accs)) != len(accs):
        raise RuntimeError("duplicate accessions in " + args.fasta)
    proteins = set(accs)

    interpro, uncharacterized, in_table = set(), set(), set()
    with gzip.open(args.uniprot, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        xref_cols = [c for c in reader.fieldnames if c not in NOT_XREF]
        for r in reader:
            a = r["Entry"]
            in_table.add(a)
            if any((r[c] or "").strip() for c in xref_cols):
                interpro.add(a)
            if (r["Protein names"] or "").strip() == UNCHARACTERIZED:
                uncharacterized.add(a)
    missing = proteins - in_table
    if missing:
        raise RuntimeError("{} FASTA accessions missing from {}, e.g. {}".format(len(missing), args.uniprot, sorted(missing)[:5]))

    sprot = set()
    with gzip.open(args.mmseqs, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if float(f[2]) > MIN_FIDENT:
                sprot.add(acc_of(f[0]))
    stray = sprot - proteins
    if stray:
        raise RuntimeError("{} MMseqs2 query ids are not FASTA accessions, e.g. {}".format(len(stray), sorted(stray)[:5]))

    n = len(proteins)
    unlabelled = proteins - interpro - sprot
    summary = json.loads(args.meta)
    summary.update(
        n_proteins=n,
        n_interpro_match=len(proteins & interpro),
        n_sprot_hit=len(sprot),
        n_unlabelled=len(unlabelled),
        unlabelled_pct=round(100 * len(unlabelled) / n, 2),
        n_uncharacterized_name=len(proteins & uncharacterized),
        uncharacterized_name_pct=round(100 * len(proteins & uncharacterized) / n, 2),
        definition="share of proteins with neither an InterPro match nor a Swiss-Prot hit above 30% identity (a hit to its own entry counts)",
        busco=[busco_numbers(b) for b in args.busco],
    )
    with open(args.out, "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
    print(json.dumps({k: summary[k] for k in ("n_proteins", "n_unlabelled", "unlabelled_pct", "uncharacterized_name_pct")}))


if __name__ == "__main__":
    sys.exit(main())
