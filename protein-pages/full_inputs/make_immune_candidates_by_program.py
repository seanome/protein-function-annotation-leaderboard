"""The table behind the index page's "Possible immune genes, by program" panel.

Reads immune_gene_candidates.tsv from olgabot/botryllus-mhc PR 16
(immune_gene_candidates/results/), which lists every Botryllus protein that at least one
program calls a possible immune gene, with one evidence column per program. It writes
immune_candidates_by_program.tsv.gz with these columns:

  protein_id   the protein's ID on these pages (FUN000002, FUN002196.2, ...), matched by
               identical sequence (the md5 in proteome_lengths.tsv.gz); for a protein of the
               old proteome with no identical FUN protein, its old ID (g21357.t1), which has
               no page here
  proteome     new or old
  name         the gene name in the proteome FASTA header, empty when the header has none
  <program>    that program's first piece of evidence, cut to 110 characters, and
               "(+N more)" when it gave more; empty when the program does not call it immune

  python make_immune_candidates_by_program.py immune_gene_candidates.tsv \
      botryllus_new_genome_protein_unique_names.fasta proteome_lengths.tsv.gz \
      immune_candidates_by_program.tsv.gz
"""

import gzip
import hashlib
import re
import sys

import polars as pl

PROGRAMS = ["funannotate_name", "interproscan", "unifire", "mmseqs2_best_hit",
            "foldseek_best_hit", "mhc_search", "kmerseek_human", "kmerseek_swissprot",
            "known_botryllus"]


def fasta_md5(path):
    out, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out[name] = hashlib.md5("".join(buf).encode()).hexdigest()
            name, buf = line[1:].split()[0].split("|")[0], []
        else:
            buf.append(line.strip().upper().replace("*", ""))
    if name:
        out[name] = hashlib.md5("".join(buf).encode()).hexdigest()
    return out


def round_e(note):
    """E=0.000073952399665297429 -> E=7.4e-05, as the page shows every E-value"""
    return re.sub(r"E=([0-9.]+(?:[eE][-+]?[0-9]+)?)", lambda m: "E=%.2g" % float(m.group(1)), note)


def short(cell):
    if not cell:
        return ""
    notes = [round_e(n) for n in cell.split(" | ")]
    first = notes[0]
    if len(first) > 110:  # shorten the protein name, keep the species, E-value and setting after it
        cut = first.find(" [") if " [" in first else first.find(" E=")
        tail = first[cut:] if cut > 0 else ""
        first = first[:max(20, 109 - len(tail))].rstrip() + "…" + tail if tail else first[:109] + "…"
    return first + (" (+%d more)" % (len(notes) - 1) if len(notes) > 1 else "")


def main(cand, fasta, lengths, out):
    md5_of = fasta_md5(fasta)
    page_of = {m: p for p, m in pl.read_csv(lengths, separator="\t").select("protein_id", "md5").iter_rows()}
    df = pl.read_csv(cand, separator="\t", infer_schema_length=0, quote_char=None).fill_null("")
    rows, missing = [], 0
    for r in df.iter_rows(named=True):
        r = {k: v.strip('"') for k, v in r.items()}
        if r["proteome"] == "old only":
            pid, prot = r["protein_id"], "old"
        else:
            pid, prot = page_of.get(md5_of.get(r["protein_id"], "")), "new"
            if pid is None:
                missing += 1
                continue
        name = "" if prot == "old" else r["protein_id"].split("_", 1)[1].rsplit(".", 1)[0]
        name = "" if name == r["protein_id"].split("_", 1)[0] else name  # FUN000001_FUN000001: no name
        rows.append({"protein_id": pid, "proteome": prot, "name": name, **{p: short(r[p]) for p in PROGRAMS}})
    if missing:
        sys.exit("%d new-proteome proteins have no page ID" % missing)
    with gzip.open(out, "wt") as fh:
        fh.write(pl.DataFrame(rows).write_csv(separator="\t", quote_style="never"))
    print(len(rows), "proteins")


if __name__ == "__main__":
    main(*sys.argv[1:5])
