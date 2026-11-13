"""Write s4pred_states.tsv: s4pred's predicted secondary structure for every Botryllus protein,
under the protein IDs the Biomni tables use, with no residues.

    python make_s4pred_table.py proteome_lengths.tsv.gz s4pred_ss.fasta.gz | gzip -n > s4pred_states.tsv.gz

s4pred 1.2.1 (single-sequence prediction of helix, strand and coil) was run on the whole proteome
by olgabot/botryllus-mhc PR 13 (branch olgabot/annotation-timing-pilot, commit edcbf8b, file
nextflow_runs/timing-pilot/s4pred_whole_proteome/s4pred_ss.fasta.gz; head job 50520678 on
Sherlock, 2026-10-31). That file has three lines per protein: the name, the residues, and one
letter per residue (H helix, E strand, C coil). It has no length limit, so the 50 proteins over
5,000 amino acids are in it too.

Only the third line is kept, since the proteome is unpublished. It is written as runs: C16E5C1E2
is 16 coil, 5 strand, 1 coil, 2 strand. A protein is matched by name, with the rule
make_unifire_table.py uses (s4pred kept the whole FASTA header and numbered every copy of a
repeated header .1, .2, ...), and then checked: the residues s4pred printed must have the length
and the MD5 that proteome_lengths.tsv gives for that ID, and every protein of that table must
appear once. The script stops at the first protein that fails.
"""
import csv
import gzip
import hashlib
import itertools
import re
import sys
from collections import Counter


def runs(states):
    """'CCCEE' -> 'C3E2'."""
    return "".join("%s%d" % (c, len(list(g))) for c, g in itertools.groupby(states))


def main(lengths, ss):
    want = {r["protein_id"]: (int(r["length"]), r["md5"]) for r in csv.DictReader(gzip.open(lengths, "rt"), delimiter="\t")}
    copies = Counter(i.split(".")[0] for i in want)
    lines = [l.rstrip("\n") for l in gzip.open(ss, "rt")]
    if len(lines) % 3:
        sys.exit("%s: %d lines, not three per protein" % (ss, len(lines)))
    out = {}
    for name, seq, states in zip(lines[0::3], lines[1::3], lines[2::3]):
        if not name.startswith(">"):
            sys.exit("expected a name line, got %r" % name[:40])
        uid = name[1:]
        gene = uid.split("_")[0]
        n = 1
        if copies[gene] > 1:
            m = re.search(r"\.(\d+)$", uid)
            if not m:
                sys.exit("%s: its gene has %d proteins but the name has no copy number" % (uid, copies[gene]))
            n = int(m.group(1))
        pid = gene if n == 1 else "%s.%d" % (gene, n)
        seq = seq.rstrip("*")
        if pid not in want:
            sys.exit("%s: no protein %s in %s" % (uid, pid, lengths))
        if pid in out:
            sys.exit("%s: protein %s appears twice" % (uid, pid))
        if (len(seq), hashlib.md5(seq.encode()).hexdigest()) != want[pid]:
            sys.exit("%s: its residues are not those of %s (length or MD5 differs)" % (uid, pid))
        if len(states) != len(seq) or set(states) - set("HEC"):
            sys.exit("%s: %d states for %d residues, letters %s" % (uid, len(states), len(seq), sorted(set(states))))
        out[pid] = runs(states)
    missing = [p for p in want if p not in out]
    if missing:
        sys.exit("%d proteins have no prediction, first %s" % (len(missing), missing[0]))
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["protein_id", "states"])
    for pid in want:
        w.writerow([pid, out[pid]])


if __name__ == "__main__":
    main(*sys.argv[1:])
