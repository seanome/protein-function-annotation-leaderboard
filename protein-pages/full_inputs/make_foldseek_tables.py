"""From Biomni's ProstT5 + Foldseek search of the whole proteome, write the best hit per protein
(for build_full.py) and every hit of the drawn test proteins (for build.py).

    python make_foldseek_tables.py botryllus_prostt5_foldseek.m8 proteome_lengths.tsv.gz \
        ../test_inputs/data/proteins.fasta foldseek_best_per_protein.tsv.gz ../test_inputs/data/prostt5_foldseek.m8.gz

The input is botryllus_prostt5_foldseek.m8 from Biomni (report_botryllus_immune_annotation_v3.md):
3Di letters predicted from each sequence by ProstT5, searched with Foldseek 10 against the
AlphaFold models of Swiss-Prot, E-value cutoff 10. No header; the columns are query, target,
fident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits, and the queries are
the Biomni protein IDs. The best hit is the lowest E-value, then the highest bit score.

The drawn proteins are named differently (FUN019272_TTN.1 is the tables' FUN019272.2), so they
are matched by the MD5 of their sequence. Their rows are written in build.py's column order
(M8_COLS, then alntmscore, which ProstT5 mode does not give, then the target).
"""
import csv
import gzip
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build import af_accession  # noqa: E402  keeps an isoform's -2


def fasta_md5(path):
    out, name, seq = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out[name] = hashlib.md5("".join(seq).encode()).hexdigest()
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        out[name] = hashlib.md5("".join(seq).encode()).hexdigest()
    return out


def main(m8, lengths, drawn_fasta, best_out, drawn_out):
    md5 = {r["protein_id"]: r["md5"] for r in csv.DictReader(gzip.open(lengths, "rt"), delimiter="\t")}
    # One Biomni ID per drawn protein: identical sequences have identical hits, so the first is enough
    by_md5 = {}
    for pid in sorted(md5):
        by_md5.setdefault(md5[pid], pid)
    drawn = {by_md5[m]: did for did, m in fasta_md5(drawn_fasta).items() if m in by_md5}
    best, n = {}, 0
    with gzip.open(drawn_out, "wt") as dw:
        for line in open(m8):
            f = line.rstrip("\n").split("\t")
            q, e, bits = f[0], float(f[10]), float(f[11])
            if q not in md5:
                sys.exit("%s is in %s but not in the proteome" % (q, m8))
            n += 1
            b = best.get(q)
            if b is None or (e, -bits) < (float(b[10]), -float(b[11])):
                best[q] = f
            if q in drawn:
                acc = af_accession(f[1]) if f[1].startswith("AF-") else f[1]
                dw.write("\t".join([drawn[q], f[1], f[2], f[3], f[6], f[7], f[8], f[9], f[10], f[11], "", acc]) + "\n")
    with gzip.open(best_out, "wt") as bw:
        w = csv.writer(bw, delimiter="\t", lineterminator="\n")
        w.writerow(["protein_id", "target", "fident", "qstart", "qend", "tstart", "tend", "evalue", "bits"])
        for q in sorted(best):
            f = best[q]
            w.writerow([q, f[1], f[2], f[6], f[7], f[8], f[9], f[10], f[11]])
    print("%d rows; best hit for %d proteins; %d of %d drawn proteins matched" % (n, len(best), len(drawn), len(fasta_md5(drawn_fasta))),
          file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:])
