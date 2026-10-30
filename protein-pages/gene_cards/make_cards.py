"""Write one plain-text evidence sheet per protein: everything the Biomni tables say about it,
with names in place of bare accessions and GO IDs. The language model reads these sheets and
nothing else, so a sheet must hold all the evidence and no protein sequence.

  python gene_cards/make_cards.py full_inputs gene_cards/out --pilot 100

writes out/evidence.jsonl (protein_id, kind, evidence text) for the chosen proteins, and with
--pilot N, out/pilot_ids.tsv: N proteins drawn so that every kind of case is in the sample.
"""
import argparse
import collections
import csv
import gzip
import json
import random
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from build_full import inp, opentext, split_list  # noqa: E402

# Proteins named in the PR discussions, always in the pilot so the answers can be checked by eye
FIXED = ["FUN008084", "FUN009411", "FUN000002", "FUN035341", "FUN006619", "FUN040347"]
# How many of each kind go into the pilot, before the fixed ones
PILOT_MIX = {"conflict_different_gene": 30, "conflict_like_only": 8, "agree": 12, "unconfirmed": 10,
             "new": 10, "immune": 10, "fantasia_only": 14}


def tsv(folder, name):
    return list(csv.DictReader(opentext(inp(folder, name)), delimiter="\t"))


def norm(name):
    """A gene name without case, hyphens or a trailing "-like": Plin2-like and PLIN2 are one name."""
    return re.sub(r"[-_ ]?like$", "", name.strip(), flags=re.I).lower().replace("-", "").replace(" ", "")


def kind(r, immune):
    if r["name_status"] == "conflict":
        return "conflict_like_only" if norm(r["existing_name"]) == norm(r["proposed_name"]) else "conflict_different_gene"
    if r["name_status"] in ("agree", "unconfirmed", "new"):
        return r["name_status"]
    if r["protein_id"] in immune:
        return "immune"
    if not r["sp_acc"] and not r["pfam_domains"]:
        return "fantasia_only"
    return "other"


def load(folder):
    ann = tsv(folder, "botryllus_protein_annotations_v2.tsv")
    immune = {r["protein_id"]: r for r in tsv(folder, "botryllus_immune_genes_v2.tsv")}
    uni = {r["accession"]: r for r in tsv(folder, "uniprot_names.tsv")}
    go = {r["go_id"]: r for r in tsv(folder, "go_names.tsv")}
    lengths = {r["protein_id"]: int(r["length"]) for r in tsv(folder, "proteome_lengths.tsv")}
    refs = {r["accession"]: r["title"] for r in tsv(folder, "allorecognition_references.tsv")}
    jh = collections.defaultdict(list)
    for r in tsv(folder, "botryllus_jackhmmer_immune.tsv"):
        jh[r["target"]].append(r)
    allo = collections.defaultdict(list)
    for f in csv.reader(open(inp(folder, "botryllus_allorecognition_hits.tsv")), delimiter="\t"):
        allo[f[0]].append(f)
    members = collections.defaultdict(list)
    rep_of = {}
    for f in csv.reader(open(inp(folder, "botryllus_mmseqs_clusters.tsv")), delimiter="\t"):
        members[f[0]].append(f[1])
        rep_of[f[1]] = f[0]
    return dict(ann=ann, by_id={r["protein_id"]: r for r in ann}, immune=immune, uni=uni, go=go,
                lengths=lengths, refs=refs, jh=jh, allo=allo, members=members, rep_of=rep_of)


def uname(d, acc):
    u = d["uni"].get(acc)
    if not u:
        return acc + " (no longer in UniProt)"
    gene = (" [gene %s]" % u["gene"]) if u["gene"] else ""
    return "%s %s%s, %s" % (acc, u["protein_name"], gene, u["organism"])


def evidence(d, pid):
    r, L = d["by_id"][pid], d["lengths"].get(pid)
    out = ["Protein %s, %s amino acids, from the Botryllus schlosseri (colonial tunicate) genome." % (pid, L)]
    out.append("Name in the genome's FASTA header (from an earlier, unrecorded BLAST-based annotation): %s"
               % (r["existing_name"] or "none"))
    out.append("Name proposed by the new annotation from the best Swiss-Prot match or Pfam domains: %s"
               % (r["proposed_name"] or "none"))
    if r["sp_acc"]:
        out.append("Best Swiss-Prot match (DIAMOND blastp --sensitive, E-value cutoff 0.001): %s; %s%% identity, "
                   "the match covers %s%% of the Botryllus protein, E-value %s." % (uname(d, r["sp_acc"]), r["sp_pident"], r["sp_qcov"], r["sp_evalue"]))
        names = [d["go"][g]["name"] for g in split_list(r["go_terms"]) if g in d["go"]]
        if names:
            out.append("GO terms of that Swiss-Prot entry (annotated on the matched protein, not measured here): " + "; ".join(names[:25])
                       + (" (and %d more)" % (len(names) - 25) if len(names) > 25 else ""))
        kw = [k.strip() for k in re.split(r";+", r["sp_keywords"]) if k.strip() and k.strip() not in ("Reference proteome", "3D-structure", "Alternative splicing")]
        if kw:
            out.append("UniProt keywords of that entry: " + ", ".join(kw))
    else:
        out.append("No Swiss-Prot match at E-value 0.001 or below.")
    if r["pfam_domains"]:
        c = collections.Counter(p.split("|")[-1] for p in split_list(r["pfam_domains"]))
        out.append("Pfam domains (hmmsearch, gathering thresholds; order along the protein not recorded): "
                   + ", ".join("%s x%d" % (n, k) if k > 1 else n for n, k in c.items()))
    else:
        out.append("No Pfam domain at the gathering threshold.")
    im = d["immune"].get(pid)
    if im:
        out.append("Immune-gene candidate, categories: %s; evidence: %s" % (im["categories"], im["evidence"]))
    if d["jh"].get(pid):
        best = {}
        for j in d["jh"][pid]:
            if j["seed_family"] not in best or float(j["e_value"]) < float(best[j["seed_family"]]["e_value"]):
                best[j["seed_family"]] = j
        out.append("jackhmmer hits from immune seed proteins (best per seed family): " + "; ".join(
            "%s (seed %s, round %s, E-value %s)" % (j["seed_family"], j["seed"], j["iteration"], j["e_value"]) for j in best.values()))
    if d["allo"].get(pid):
        top = sorted(d["allo"][pid], key=lambda f: -float(f[5]))[:3]
        out.append("DIAMOND matches to Botryllus allorecognition reference proteins (Fester family and BHF): " + "; ".join(
            "%s %s, %s%% identity over %s columns, E-value %s, reference %s aa" % (f[1], d["refs"].get(f[1], "?"), f[2], f[3], f[4], f[7]) for f in top))
    if r["fantasia_uniprot"]:
        gos = [g.split("|", 1)[1] for g in split_list(r["fantasia_go_terms"]) if "|" in g]
        out.append("FANTASIA (protein-language-model embedding, nearest experimentally annotated UniProt protein, no sequence alignment): "
                   "nearest protein %s, cosine distance %.2f, reliability index %.2f (1 is best). GO terms copied from it: %s"
                   % (uname(d, r["fantasia_uniprot"]), float(r["fantasia_distance"]), float(r["fantasia_ri"]), "; ".join(gos) or "none"))
    rep = d["rep_of"].get(pid)
    mem = d["members"].get(rep, [])
    if len(mem) > 1:
        named = [m for m in mem if m != pid and (d["by_id"][m]["existing_name"] or d["by_id"][m]["proposed_name"])]
        out.append("MMseqs2 cluster of %d Botryllus proteins (representative %s). Names of other members: %s"
                   % (len(mem), rep, "; ".join("%s: header %s, proposed %s" % (m, d["by_id"][m]["existing_name"] or "none", d["by_id"][m]["proposed_name"] or "none")
                                              for m in named[:6]) or "none named"))
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--pilot", type=int, default=0, help="draw this many proteins covering every kind of case")
    ap.add_argument("--seed", type=int, default=20261030)
    a = ap.parse_args()
    d = load(a.inputs)
    kinds = {r["protein_id"]: kind(r, d["immune"]) for r in d["ann"]}
    a.out.mkdir(parents=True, exist_ok=True)
    if a.pilot:
        rng = random.Random(a.seed)
        ids = list(FIXED)
        for k, n in PILOT_MIX.items():
            pool = sorted(p for p, kk in kinds.items() if kk == k and p not in ids)
            ids += rng.sample(pool, n)
        ids = ids[:a.pilot]
        with open(a.out / "pilot_ids.tsv", "w") as fh:
            fh.write("protein_id\tkind\n" + "".join("%s\t%s\n" % (p, kinds[p]) for p in ids))
    else:
        ids = [r["protein_id"] for r in d["ann"]]
    with open(a.out / "evidence.jsonl", "w") as fh:
        for p in ids:
            fh.write(json.dumps({"protein_id": p, "kind": kinds[p], "evidence": evidence(d, p)}) + "\n")
    print(len(ids), "evidence sheets;", collections.Counter(kinds[p] for p in ids), file=sys.stderr)
    print("all proteins by kind:", collections.Counter(kinds.values()), file=sys.stderr)


if __name__ == "__main__":
    main()
