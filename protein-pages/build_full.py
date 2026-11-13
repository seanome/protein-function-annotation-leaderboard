#!/usr/bin/env python3
"""Pages for every protein of the Botryllus proteome (45,339), from the Biomni annotation tables.

    python build_full.py full_inputs --drawn site --drawn-fasta test_inputs/data/proteins.fasta --out site_full

Reads, from the input folder (see full_inputs/README.md for where each file came from):
    botryllus_protein_annotations.tsv[.gz]   one row per protein: names, best Swiss-Prot hit, GO, Pfam
    botryllus_immune_genes.tsv                immune-gene candidates, their categories and evidence
    botryllus_jackhmmer_immune.tsv            jackhmmer hits of 20 immune seed proteins
    botryllus_mmseqs_clusters.tsv             MMseqs2 clusters (representative, member)
    proteome_lengths.tsv[.gz]                 each protein's length and the MD5 of its sequence

None of these tables gives where on the protein a hit lies, so these pages list hits in tables
and draw nothing on a ruler. Proteins also in a site made by build.py (--drawn) open on that
page, where every program's hits are drawn, with these tables added under the drawing;
they are matched by sequence (--drawn-fasta), since the two FASTA files number isoforms differently.

Writes into --out:
    index.html          every protein: search, filters, counts
    protein.html        one page that shows any protein, named after "#" in its address
    full/rows.js        what the index shows for each protein, packed (assets/gz.js reads it)
    full/meta.js        the names protein pages need: immune categories, name checks, verdicts
    full/NNN.js         the rest of each protein's data, in 96 files; a protein's file number is a
                        hash of its ID (shard_of), so protein.html loads one file
    full_summary.json   the counts on the index
    drawn.html, proteins/, targets/, assets/, search.js  copied from --drawn; drawn.html is
                        build.py's index of the proteins drawn there
"""
import argparse
import base64
import csv
import gzip
import hashlib
import json
import re
import shutil
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import build  # esc, CSS_LINKS, the asset version

N_SHARDS = 96
MEMBERS_MAX = 200  # cluster members listed on a protein page; the rest are counted

# The Biomni table's name_status, checked against existing_name and proposed_name on every row
# (tests/test_build_full.py): same name ignoring case, different name, an existing name and no
# proposal, a proposal and no existing name.
def names_conflict(existing, proposed):
    """True when two names differ by more than capitals, punctuation and a trailing -like."""
    key = lambda s: re.sub(r"[^a-z0-9]", "", re.sub(r"[- ]?like$", "", s.strip().lower()))
    return bool(existing and proposed) and key(existing) != key(proposed)


ONLY_KMERSEEK_E = 1e-3  # the "only kmerseek" filters keep this E-value or below; rows.js holds the hit for those


# Gene groups for the index page's "gene group" filter. A protein is in a group when the best hit of
# DIAMOND, MMseqs2, either Foldseek search or kmerseek is a Swiss-Prot protein whose name matches the
# group's pattern, at any E-value the pages show: the filter is for looking, so weak hits are kept and
# each row says which program and E-value put it there. Only the protein's main name is matched: "Isoform 2 of "
# and the other names in brackets are dropped first, or Calnexin would be an MHC protein for its other name.
# MHC class I and class II are also listed on their own; "HLA and other MHC" holds both.
# The aging group is read another way: the best hit's gene symbol is one of the 307 human genes of
# GenAge build 21 (full_inputs/genage_human.tsv), capitals ignored, in whichever species the hit is.
# The allorecognition group is the one group that is not a match to another species: Biomni's DIAMOND
# search of the proteome against 96 Botryllus schlosseri proteins from NCBI (Fester family, Fester
# co-receptor and one partial BHF), so a hit can be the protein itself.
GENAGE, ALLO = "genage", "allorecognition"  # groups read from a table, not from a name pattern
GENE_GROUPS = [
    ("Period", r"^Period (circadian|clock) protein"),
    ("Cry and photolyases", r"^Cryptochrome|photo-?lyase"),
    ("other clock genes", r"^(Circadian locomoter output cycles protein kaput|Basic helix-loop-helix ARNT-like protein|"
     r"Neuronal PAS domain-containing protein 2|Nuclear receptor subfamily 1 group D|Nuclear receptor ROR-|Protein timeless|"
     r"TIMELESS-interacting protein|Casein kinase I isoform (delta|epsilon)|Nuclear factor interleukin-3-regulated protein|"
     r"Thyrotroph embryonic factor|Transcription factor VBP|Hepatic leukemia factor|D site-binding protein|"
     r"Circadian-associated transcriptional repressor|Protein cycle)"),
    ("HLA and other MHC", r"^(?!Minor ).*histocompatibility|^HLA class|^Beta-2-microglobulin|^MHC class I polypeptide-related"),
    # beta-2-microglobulin is the second chain of a class I molecule; MIC-A and MIC-B are class I relatives
    ("MHC class I", r"^(?!Minor ).*[Cc]lass I histocompatibility|^Beta-2-microglobulin|^MHC class I polypeptide-related"),
    # the gamma chain is the invariant chain (CD74), which escorts class II molecules
    ("MHC class II", r"^(?!Minor ).*[Cc]lass II histocompatibility"),
    ("allorecognition (Fester, Fester co-receptor, BHF)", ALLO),
    ("aging (GenAge human genes)", GENAGE),
]


def read_allorecognition(folder):
    """protein -> [E-value, reference accession, its NCBI title] of its best hit (lowest E-value, then
    highest bit score) in botryllus_allorecognition_hits.tsv; empty when the table is not there."""
    hits, titles = inp(folder, "botryllus_allorecognition_hits.tsv"), inp(folder, "allorecognition_references.tsv")
    if not hits.exists():
        return {}
    title = {r["accession"]: r["title"].replace(" [Botryllus schlosseri]", "") for r in csv.DictReader(opentext(titles), delimiter="\t")} if titles.exists() else {}
    best = {}
    for f in csv.reader(opentext(hits), delimiter="\t"):  # no header: query, reference, % identity, columns, E-value, bit score, ...
        key = (float(f[4]), -float(f[5]))
        if f[0] not in best or key < best[f[0]][0]:
            best[f[0]] = (key, [float(f[4]), f[1], title.get(f[1], "")])
    return {pid: v for pid, (_, v) in best.items()}


def read_genage(folder):
    """The GenAge human gene symbols, lower case; empty when the table is not there."""
    p = inp(folder, "genage_human.tsv")
    return {r["symbol"].lower() for r in csv.DictReader(opentext(p), delimiter="\t")} if p.exists() else set()


def gene_groups(r, pid, uni, aging, mm, fs, fst, km, allo=None):
    """[[group index, program, E-value, accession, name], ...] for one protein, best E-value first.
    uni: accession -> names (uniprot_names.tsv); aging: read_genage(); allo: read_allorecognition()."""
    hits = [("DIAMOND", r["sp_evalue"], r["sp_acc"], r["sp_name"])] if r["sp_acc"] else []
    for prog, d in (("MMseqs2", mm), ("ProstT5 + Foldseek", fs), ("Foldseek on structure", fst)):
        if d and d.get(pid):
            hits.append((prog, d[pid][4], d[pid][0], d[pid][6]))
    if km and km.get(pid):
        hits.append(("kmerseek", km[pid][1], km[pid][2], km[pid][3]))
    out = []
    for i, (_, pat) in enumerate(GENE_GROUPS):
        if pat == ALLO:
            if allo and allo.get(pid):
                out.append([i, "DIAMOND against Botryllus proteins", allo[pid][0], allo[pid][1], allo[pid][2] or "no title"])
            continue
        if pat == GENAGE:
            for prog, e, acc, name in hits:
                gene = uni.get(acc, {}).get("gene", "")
                if gene.lower() in aging:
                    out.append([i, prog, float(e), acc, "%s, gene %s" % (name or "no name", gene)])
            continue
        for prog, e, acc, name in hits:
            if name and re.search(pat, re.split(r" [(\[]", re.sub(r"^Isoform \S+ of ", "", name))[0]):
                out.append([i, prog, float(e), acc, name])
    return sorted(out, key=lambda x: (x[2], x[0]))


def other_program_found(r, pid, mm, fs, fst, jh, ips, uf):
    """True when any program but kmerseek found something on this protein, each at the cutoff the pages
    use. InterProScan's predictions (disorder, coiled coil, membrane, signal) are not findings."""
    sp = r["sp_acc"] and r["sp_evalue"] and float(r["sp_evalue"]) <= 1e-3
    fa = r.get("fantasia_distance") and float(r["fantasia_distance"]) < build.FANTASIA_MAX_DIST
    return bool(sp or r["pfam_domains"] or fa or any(d and d.get(pid) for d in (mm, fs, fst, jh, uf)) or
                any(x[0] not in build.IPS_DRAWN for x in (ips or {}).get(pid, [])))


def read_rules(folder):
    """{rule: ({signature: [(set, condition), ...]}, [conditions in each set])} from
    unifire_rule_conditions.tsv. A set holds when each of its conditions has one of its signatures."""
    sets = defaultdict(dict)
    if inp(folder, "unifire_rule_conditions.tsv").exists():
        with opentext(inp(folder, "unifire_rule_conditions.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sets[r["rule"]].setdefault(r["set"], []).append(r["signatures"].split(","))
    rules = {}
    for rule, by_set in sets.items():
        where, sizes = defaultdict(list), []
        for s, conds in enumerate(by_set.values()):
            sizes.append(len(conds))
            for c, sigs in enumerate(conds):
                for g in sigs:
                    where[g].append((s, c))
        rules[rule] = (where, sizes)
    return rules


def rule_places(preds, matches, rules):
    """UniFIRE's rows with a seventh field: which of this protein's InterProScan matches (their index in
    the sorted list) are signatures of a condition set that holds here. Rows with their own position keep it."""
    sig = defaultdict(set)  # each name a signature goes by -> its matches; G3DSA:3.40.50.300, PTHR1:SF5
    for i, m in enumerate(matches):
        for g in (m[1], m[1].split(":")[0], m[1].replace("G3DSA:", ""), m[6]):
            if g:
                sig[g].add(i)
    out, memo = [], {}
    for u in preds:
        if u[4] is not None:
            out.append(list(u[:6]))
            continue
        if u[1] not in memo:
            where, sizes = rules.get(u[1], ({}, []))
            seen = defaultdict(lambda: defaultdict(set))  # set -> condition -> matches
            for g, idx in sig.items():
                for s, c in where.get(g, ()):
                    seen[s][c] |= idx
            full = [s for s, conds in seen.items() if len(conds) == sizes[s]]
            # UniFIRE ran an earlier release of the rules than UniProt's API serves, so a set may no longer hold in
            # full: then the set that comes closest says where the rule points
            near = full or sorted(seen, key=lambda s: (-len(seen[s]) / sizes[s], s))[:1]
            memo[u[1]] = sorted({i for s in near for v in seen[s].values() for i in v})
        out.append(list(u[:6]) + [memo[u[1]]])
    return out


def unifire_hits(preds, matches, length):
    """One hit for each rule and place, for the grid: a site's own position, else where the rule's
    signatures lie (assets/protein-dense.js does the same on protein.html)."""
    hits = {}
    for u in preds:
        at = [(u[4], u[5] or u[4], "the rule's own position")] if u[4] else [
            (matches[i][3], matches[i][4], matches[i][1] + (" / " + matches[i][6] if matches[i][6] else "")) for i in (u[6] if len(u) > 6 else [])]
        for a, b, why in at or [(1, length, "whole protein: UniProt's current version of this rule names no signature that InterProScan found here")]:
            h = hits.setdefault((u[1], a, b), dict(k="unifire", q=[a, b], s=[u[0], u[1], u[2], why], e=None, f=[], v=[]))
            if u[3] not in h["v"]:
                h["v"].append(u[3])
    for h in hits.values():
        v = h.pop("v")
        h.update(l=v[0] + (" +%d" % (len(v) - 1) if len(v) > 1 else ""), n="; ".join(v))
    return list(hits.values())


def read_s4pred(folder, lengths):
    """protein -> s4pred's predicted secondary structure as runs ('C16E5C1E2': 16 coil, 5 strand, ...), from
    s4pred_states.tsv (full_inputs/make_s4pred_table.py); None when the table is not there. Each protein's runs
    must add up to its length, and the summary's three shares are of all residues."""
    if not inp(folder, "s4pred_states.tsv").exists():
        return None, None
    csv.field_size_limit(10 ** 9)  # the longest protein's runs are longer than csv's default limit
    ss, tot = {}, Counter()
    for r in csv.DictReader(opentext(inp(folder, "s4pred_states.tsv")), delimiter="\t"):
        pid, runs = r["protein_id"], re.findall(r"([HEC])(\d+)", r["states"])
        if pid not in lengths:
            sys.exit("%s is in s4pred_states.tsv but not in the proteome" % pid)
        if "".join(c + n for c, n in runs) != r["states"] or sum(int(n) for _, n in runs) != lengths[pid][0]:
            sys.exit("s4pred_states.tsv: the runs of %s do not add up to its %d residues" % (pid, lengths[pid][0]))
        ss[pid] = r["states"]
        for c, n in runs:
            tot[c] += int(n)
    n = sum(tot.values())
    return ss, dict(proteins=len(ss), residues=n, helix_pct=round(100 * tot["H"] / n, 1), strand_pct=round(100 * tot["E"] / n, 1),
                    coil_pct=round(100 * tot["C"] / n, 1))


def drawn_additions(r, d, length):
    """The script that adds DIAMOND's realigned best hit and UniFIRE's predictions to a drawn page's
    data, before assets/dense-bhf.js reads it. A DIAMOND hit too long to realign keeps its row with no position."""
    hits = []
    move = []
    if d.get("dm") or not r[6]:
        move.append("w_diamond")
    if d.get("dm"):
        a = d["dm"]
        hits.append(dict(k="diamond", l=r[6], t=r[6], q=a[:2], r=a[2:4], id=a[4:6], e=float(d["sp_e"]), f=None, n=r[5], o=d["sp_org"],
                         s=[d["sp_e"], "" if r[7] is None else str(r[7]), "" if r[8] is None else str(r[8])], tl=d["tl"].get(r[6])))
    if isinstance(d.get("uf"), list):
        move.append("w_unifire")
        hits += unifire_hits(d["uf"], d["ips"] if isinstance(d["ips"], list) else [], length)
    tracks = {k: dict(name=DRAWN_NAMES[k], color=build.WHOLE_TRACKS[k]["color"], shape=build.WHOLE_TRACKS[k]["shape"],
                      cat=build.WHOLE_TRACKS[k]["cat"], p=k) for k in ("diamond", "unifire") if "w_" + k in move}
    return dict(tracks=tracks, hits=hits, names=DRAWN_SCORE_NAMES, nohit=DRAWN_NOHIT, move=move, how=DRAWN_HOW, ss=d.get("ss"),
                twelve="; DIAMOND's hit is realigned, and UniFIRE is placed where its rules' signatures lie")


def pack_drawn(html, add, ver):
    """A drawn page with its data packed and its scripts started by assets/gz.js: the rows first, Mol* after."""
    i = html.index("<script>const PAGE=")
    j = html.index("</body>")
    m = re.match(r"<script>const PAGE=(.*?);window\.WHOLE=(.*?);window\.NOTRUN=(.*?);window\.HOW=(.*?);</script>", html[i:j], re.S)
    if not m:
        sys.exit("could not read the data of a drawn page")
    data = dict(zip(("page", "whole", "notrun", "how"), (json.loads(g) for g in m.groups())), add=add)
    srcs = re.findall(r'<script src="([^"]+)"></script>', html[i + m.end():j])
    mol = [x for x in srcs if "molstar" in x or x.endswith(".pdb.js") or "structure.js" in x]
    css = re.search(r'<link rel="stylesheet" href="([^"]*molstar[^"]*)">', html)
    head = html[:i].replace(css.group(0), "") if css else html[:i]
    cfg = dict(z=packed(data), base="../", dense="../assets/dense-bhf.js", molstar=[css.group(1)] + mol if css else None)
    return '%s<script src="../assets/gz.js%s"></script><script>bootPage(%s)</script>%s' % (head, ver, dumps(cfg), html[j:])


def dumps(x):
    return json.dumps(x, separators=(",", ":")).replace("</", "<\\/")


DRAWN_NAMES = {"diamond": "DIAMOND, best hit", "unifire": "UniFIRE (UniProt's rules)"}
DRAWN_SCORE_NAMES = {"diamond": ["evalue", "% identity (DIAMOND's own)", "% of your protein (DIAMOND's own)"],
                     "unifire": ["source", "rule", "what it predicts", "where it points"]}
DRAWN_NOHIT = {"diamond": "E ≤ 0.001: its table keeps nothing weaker", "unifire": "no rule fired"}
DRAWN_HOW = {"diamond": "best hit only: run by Biomni; its table has no positions, so this pair was aligned again (BLOSUM62, local, gap open -11, extend -1)",
             "unifire": "run by Biomni on InterProScan's matches. Rules work on signatures, so each prediction is placed where the rule's signatures lie on this protein; a site uses the rule's own position"}


NAME_CHECK = [
    ("", "no name, existing or proposed"),
    ("agree", "proposed name is the existing name"),
    ("conflict", "proposed name differs from the existing name"),
    ("unconfirmed", "existing name, nothing proposed"),
    ("new", "proposed name, no existing name"),
]


# How a gene card's name_verdict reads on the page (gene_cards/prompt.md defines the choices)
VERDICT_TEXT = {
    "existing": "the header name is better supported",
    "proposed": "the proposed name is better supported",
    "same_gene": "both names mean the same gene",
    "family_only": "both are in one family and the evidence cannot tell which member",
    "neither": "the evidence supports neither name",
}


def opentext(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def inp(folder, name):
    p = folder / name
    return p if p.exists() else folder / (name + ".gz")


def shard_of(pid):
    """File number of a protein's data; protein.js computes the same (Math.imul(h, 31) + code)."""
    h = 0
    for ch in pid:
        h = (h * 31 + ord(ch)) % 2 ** 32
    return h % N_SHARDS


def split_list(s, sep=";"):
    return [x.strip() for x in s.split(sep) if x.strip()]


def pfam_list(s):
    """'PF05729.19|NACHT;PF00560.38|LRR_1' -> [[accession, name], ...], repeats kept in table order."""
    out = []
    for x in split_list(s):
        acc, _, name = x.partition("|")
        out.append([acc, name])
    return out


def num(s):
    try:
        return float(s)
    except ValueError:
        return None


def read_fasta_md5(path):
    seqs, cur = {}, None
    for line in opentext(path):
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif line:
            seqs[cur].append(line)
    return {k: hashlib.md5("".join(v).rstrip("*").encode()).hexdigest() for k, v in seqs.items()}


KMERSEEK_MAX_E = 10.0  # a kmerseek best hit is shown when its E-value is at most this
STRUCT_MAX_AA = 200  # Biomni predicted 3-D structures only for proteins of at most this length
FOLDSEEK_MAX_E = 1e-3  # a Foldseek best hit is shown when its E-value is at most this, the cut the Biomni report calls confident
UNIFIRE_MAX_AA = 5000  # longer proteins were not run through InterProScan and UniFIRE (botryllus-mhc PR 13)


# Biomni also searched every protein against 95 Fester and 1 BHF protein, all from Botryllus itself
# (botryllus_allorecognition_hits.tsv, allorecognition_references.tsv). Matches within the same
# species are not shown on the protein pages, and the immune category those hits set is dropped from
# the immune table. The index page's allorecognition gene group (read_allorecognition) is the one
# place the two files are read.
NOT_READ = set()  # tables in the folder that the build does not read
DROPPED_CATEGORIES = {"allorecognition"}


def without_dropped_categories(immune):
    """The immune table without DROPPED_CATEGORIES and their evidence; a protein left with no
    category is no longer an immune-gene candidate."""
    out = {}
    for pid, r in immune.items():
        cats = [c for c in split_list(r["categories"]) if c not in DROPPED_CATEGORIES]
        if cats:
            ev = [e for e in r["evidence"].split(" | ") if e.partition(":")[0].strip() not in DROPPED_CATEGORIES]
            out[pid] = dict(r, categories=";".join(cats), evidence=" | ".join(ev))
    return out


def load(folder):
    tsv = lambda name: list(csv.DictReader(opentext(inp(folder, name)), delimiter="\t"))
    v2 = inp(folder, "botryllus_protein_annotations_v2.tsv").exists()  # v2 adds the FANTASIA columns
    ann = tsv("botryllus_protein_annotations_v2.tsv" if v2 else "botryllus_protein_annotations.tsv")
    lengths = {r["protein_id"]: (int(r["length"]), r["md5"]) for r in tsv("proteome_lengths.tsv")}
    ids = [r["protein_id"] for r in ann]
    if len(set(ids)) != len(ids):
        sys.exit("botryllus_protein_annotations.tsv names a protein twice")
    if set(ids) != set(lengths):
        sys.exit("proteome_lengths.tsv and the annotation table name different proteins")
    immune = without_dropped_categories({r["protein_id"]: r for r in tsv("botryllus_immune_genes_v2.tsv" if v2 else "botryllus_immune_genes.tsv")})
    jh = defaultdict(list)
    for r in tsv("botryllus_jackhmmer_immune.tsv"):
        jh[r["target"]].append([r["seed_family"], r["seed"], int(r["iteration"]), r["e_value"], r["score"]])
    rep_of, members = {}, defaultdict(list)
    for f in csv.reader(opentext(inp(folder, "botryllus_mmseqs_clusters.tsv")), delimiter="\t"):
        rep_of[f[1]] = f[0]
        members[f[0]].append(f[1])
    missing = set(lengths) - set(rep_of)
    if missing:
        sys.exit("%d proteins have no row in botryllus_mmseqs_clusters.tsv, e.g. %s" % (len(missing), sorted(missing)[0]))
    for k in (set(immune) | set(jh) | set(rep_of)) - set(lengths):
        sys.exit("%s is in a table but not in the proteome" % k)
    # Names fetched by gene_cards/fetch_names.py; optional, the pages show bare IDs without them
    uni = {r["accession"]: r for r in tsv("uniprot_names.tsv")} if inp(folder, "uniprot_names.tsv").exists() else {}
    go = {r["go_id"]: r["name"] for r in tsv("go_names.tsv")} if inp(folder, "go_names.tsv").exists() else {}
    # UniFIRE predictions (full_inputs/make_unifire_table.py); optional
    uf = None
    if inp(folder, "unifire_predictions.tsv").exists():
        uf = defaultdict(list)
        for r in tsv("unifire_predictions.tsv"):
            if r["protein_id"] not in lengths:
                sys.exit("%s is in unifire_predictions.tsv but not in the proteome" % r["protein_id"])
            v = r["value"]
            if r["type"] == "xref.GO" and go.get(v):
                v = "%s %s" % (v, go[v])
            uf[r["protein_id"]].append([r["source"], r["rule"], r["type"], v, int(r["start"]) if r["start"] else None,
                                        int(r["end"]) if r["end"] else None])
    # Best kmerseek hit against Swiss-Prot (full_inputs/make_kmerseek_table.py); optional
    km = None
    if inp(folder, "kmerseek_best_per_protein.tsv").exists():
        km = {}
        for r in tsv("kmerseek_best_per_protein.tsv"):
            if r["protein_id"] not in lengths:
                sys.exit("%s is in kmerseek_best_per_protein.tsv but not in the proteome" % r["protein_id"])
            if r["evalue"] and float(r["evalue"]) <= KMERSEEK_MAX_E:
                km[r["protein_id"]] = ["%s, k = %s" % (r["alphabet"], r["ksize"]), float(r["evalue"]), r["target"],
                                       r["target_name"], r["target_organism"], int(r["n_settings"])]
                if r.get("qstart"):  # where it lies, and the Swiss-Prot residues (make_kmerseek_table.py)
                    km[r["protein_id"]] += [r["alphabet"], int(r["ksize"]), int(r["qstart"]), int(r["qend"]), int(r["tstart"]), int(r["tend"]), r["target_subseq"]]
    # Whole-proteome InterProScan matches and best MMseqs2 hit per protein (make_interproscan_table.py,
    # mmseqs2_best_per_protein.tsv); optional
    ips = None
    if inp(folder, "interproscan.tsv").exists():
        ips = defaultdict(list)
        for f in csv.reader(opentext(inp(folder, "interproscan.tsv")), delimiter="\t"):
            ips[f[0]].append([f[3], f[4], f[5] if f[5] != "-" else "", int(f[6]), int(f[7]), f[8] if f[8] != "-" else "",
                              f[11] if f[11] != "-" else "", f[12] if f[12] != "-" else ""])
    mm = None
    if inp(folder, "mmseqs2_best_per_protein.tsv").exists():
        mm = {r["protein_id"]: [r["target"], float(r["fident"]), int(r["qstart"]), int(r["qend"]), r["evalue"], r["bits"],
                                uni.get(r["target"], {}).get("protein_name", ""), uni.get(r["target"], {}).get("organism", ""),
                                int(r["tstart"]), int(r["tend"])]  # last two: the range on the other protein
              for r in tsv("mmseqs2_best_per_protein.tsv")}
    # Best Foldseek hit per protein (make_foldseek_tables.py): ProstT5 + Foldseek on the whole
    # proteome, and Foldseek on predicted structures of proteins up to STRUCT_MAX_AA; optional
    def foldseek(name):
        if not inp(folder, name).exists():
            return None
        out = {}
        for r in tsv(name):
            if r["protein_id"] not in lengths:
                sys.exit("%s is in %s but not in the proteome" % (r["protein_id"], name))
            if name.startswith("foldseek_struct") and lengths[r["protein_id"]][0] > STRUCT_MAX_AA:
                sys.exit("%s is in %s but longer than %d aa" % (r["protein_id"], name, STRUCT_MAX_AA))
            if float(r["evalue"]) <= 10:
                acc = build.af_accession(r["target"])  # AF-P12345-F1-model_v6, AF-C3VPR6-2-F1-model_v6
                out[r["protein_id"]] = [acc, float(r["fident"]), int(r["qstart"]), int(r["qend"]), r["evalue"], r["bits"],
                                        uni.get(acc, {}).get("protein_name", ""), uni.get(acc, {}).get("organism", "")]
                if "lddt" in r:  # Foldseek on structures: TM-score and LDDT of the alignment
                    out[r["protein_id"]] += [float(r["alntmscore"]), float(r["lddt"])]
                out[r["protein_id"]] += [int(r["tstart"]), int(r["tend"])]  # last two: the range on the other protein
        return out
    fs = foldseek("foldseek_best_per_protein.tsv")
    fst = foldseek("foldseek_struct_best_per_protein.tsv")
    # DIAMOND's best hit realigned, for its position (make_diamond_positions.py): skipped pairs have no entry
    dm = {}
    if inp(folder, "diamond_positions.tsv").exists():
        dm = {r["protein_id"]: [int(r[c]) for c in ("qstart", "qend", "tstart", "tend", "identical", "columns")]
              for r in tsv("diamond_positions.tsv") if r["qstart"]}
    return ann, lengths, immune, jh, rep_of, members, uni, go, uf, km, ips, mm, fs, fst, dm


def packed(obj):
    """obj as gzip-compressed JSON in base64, for a data file or a drawn page; assets/gz.js reads it.
    A protein page then downloads a few hundred KB rather than every protein's table."""
    return base64.b64encode(gzip.compress(json.dumps(obj, separators=(",", ":")).encode(), 9, mtime=0)).decode()


def add_targets(out, tseqs):
    """Merge the best hits' Swiss-Prot sequences into the target files (targets/NNN.js, build.py's
    write_target_shards), keeping build.py's file for each accession, so no file is added. Returns
    the number of target files (0 if there are none)."""
    files = sorted((out / "targets").glob("*.js")) if (out / "targets").is_dir() else []
    if not files:
        if tseqs:
            (out / "targets").mkdir(exist_ok=True)
            files = [out / "targets" / ("%03d.js" % i) for i in range(max(1, len(tseqs) // 300))]
            for f in files:
                f.write_text("(window.TSEQ=window.TSEQ||{}),Object.assign(TSEQ,{});")
        else:
            return 0
    n = len(files)
    pre, post = "(window.TSEQ=window.TSEQ||{}),Object.assign(TSEQ,", ");"
    by = {}
    for acc, s in tseqs.items():
        by.setdefault(build.shard_of(acc, n), {})[acc] = s
    for i, f in enumerate(files):
        if i not in by:
            continue
        text = f.read_text()
        if not (text.startswith(pre) and text.endswith(post)):
            sys.exit("%s is not a target file build.py wrote" % f)
        d = json.loads(text[len(pre):-len(post)])
        d.update({k: v for k, v in by[i].items() if k not in d})
        f.write_text(pre + json.dumps(d, separators=(",", ":"), sort_keys=True) + post)
    return n


def read_cards(paths):
    """Gene cards from gene_cards/run_cards.py; a later file wins for a protein in two."""
    cards = {}
    for p in paths or []:
        for line in open(p):
            c = json.loads(line)
            cards[c["protein_id"]] = c
    return cards


def evidence_list(s):
    """'GO immune process:GO:0045087,GO:0051607 | NLR:Pfam:PF05729.19' -> [[category, detail], ...]"""
    out = []
    for part in s.split(" | "):
        cat, _, det = part.partition(":")
        if cat.strip():
            out.append([cat.strip(), det.strip()])
    return out


def read_proteome(path, lengths):
    """The proteome FASTA's sequences under the Biomni IDs (make_proteome_lengths.py's naming: the
    first copy of a gene's header is FUN015821, later ones FUN015821.2, ...), each checked against
    proteome_lengths.tsv's length and MD5, so a sequence can only land on its own protein."""
    seen, recs = Counter(), []
    for line in opentext(path):
        line = line.strip()
        if line.startswith(">"):
            gene = line[1:].split()[0].split("_")[0]
            seen[gene] += 1
            recs.append([gene if seen[gene] == 1 else "%s.%d" % (gene, seen[gene]), []])
        elif line:
            recs[-1][1].append(line)
    out = {}
    for pid, parts in recs:
        s = "".join(parts).rstrip("*")
        if pid not in lengths or lengths[pid] != (len(s), hashlib.md5(s.encode()).hexdigest()):
            sys.exit("%s in %s does not match proteome_lengths.tsv (length and MD5)" % (pid, path))
        out[pid] = s
    if len(out) != len(lengths):
        sys.exit("%s has %d proteins, proteome_lengths.tsv %d" % (path, len(out), len(lengths)))
    return out


def read_swissprot(path, want):
    """Sequences of the accessions in want, from a UniProt FASTA (headers sp|P12345|NAME_SPECIES ...)."""
    out, cur = {}, None
    for line in opentext(path):
        if line.startswith(">"):
            acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
            cur = acc if acc in want else None
            if cur:
                out[cur] = []
        elif cur:
            out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def build_site(folder, out, drawn=None, drawn_fasta=None, cards=None, proteome=None, swissprot=None):
    t0 = time.time()
    ann, lengths, immune, jh, rep_of, members, uni, go, uf, km, ips, mm, fs, fst, dm = load(folder)
    seqs = read_proteome(proteome, lengths) if proteome else {}
    ss, ss_summary = read_s4pred(folder, lengths)
    # the Swiss-Prot proteins of the best hits with a position (MMseqs2, ProstT5 + Foldseek, Foldseek on
    # structure), so protein.html can align the two proteins; isoforms (P12345-2) are not in uniprot_sprot.fasta
    best_targets = {r["sp_acc"] for r in ann if r["sp_acc"] and r["protein_id"] in dm} | {v[0] for d in (mm, fs, fst) if d for v in d.values() if v and float(v[4]) <= FOLDSEEK_MAX_E}
    tseqs = read_swissprot(swissprot, best_targets | {v[2] for v in (km or {}).values()}) if swissprot else {}
    tlen = {a: len(s) for a, s in tseqs.items()}  # a kmerseek hit keeps only its region's residues, so the sequence is left out
    tseqs = {a: s for a, s in tseqs.items() if a in best_targets}
    rules = read_rules(folder)
    aging, allo = read_genage(folder), read_allorecognition(folder)
    cards = read_cards(cards)
    cats = sorted({c for r in immune.values() for c in split_list(r["categories"])},
                  key=lambda c: (-sum(c in split_list(r["categories"]) for r in immune.values()), c))
    if len(cats) > 31:
        sys.exit("more immune categories than fit in one number")
    status_code = {s: i for i, (s, _) in enumerate(NAME_CHECK)}

    # Proteins also drawn by build.py, matched by sequence
    drawn_page = {}
    if drawn:
        by_md5 = defaultdict(list)
        for pid, (_, md5) in lengths.items():
            by_md5[md5].append(pid)
        for did, md5 in read_fasta_md5(drawn_fasta).items():
            page = "proteins/%s.html" % build.safe(did)
            if not (drawn / page).exists():
                continue
            for pid in by_md5.get(md5, []):
                drawn_page[pid] = page

    name = {r["protein_id"]: r["proposed_name"] or r["existing_name"] for r in ann}
    # a gene's proteins in FASTA order: FUN1, FUN1.2, FUN1.3 (make_proteome_lengths.py names them)
    gene_proteins = defaultdict(list)
    for r in ann:
        gene_proteins[r["protein_id"].split(".")[0]].append(r["protein_id"])
    rows, shards, n_no_other = [], defaultdict(dict), 0
    cache = {}

    def strong(d):  # the Foldseek hits the counts use: E-value FOLDSEEK_MAX_E or below; the pages show every best hit
        if d is None:
            return None
        if id(d) not in cache:
            cache[id(d)] = {k: v for k, v in d.items() if float(v[4]) <= FOLDSEEK_MAX_E}
        return cache[id(d)]

    pfam_ix = {}  # Pfam accession -> its index in the data file's pfam list: [accession, InterProScan's description]
    for r in ann:
        pid = r["protein_id"]
        im = immune.get(pid)
        mask = 0
        for c in split_list(im["categories"]) if im else []:
            mask |= 1 << cats.index(c)
        rep = rep_of[pid]
        size = len(members[rep])
        pf = pfam_list(r["pfam_domains"])
        found = other_program_found(r, pid, mm, strong(fs), strong(fst), jh, ips, uf)
        k = None if found or not km else km.get(pid)
        rows.append([pid, lengths[pid][0], r["existing_name"], r["proposed_name"], status_code[r["name_status"]],
                     r["sp_name"], r["sp_acc"], num(r["sp_pident"]), num(r["sp_qcov"]), r["sp_evalue"],
                     ", ".join(dict.fromkeys(n for _, n in pf)), mask,
                     "", size, drawn_page.get(pid, ""),  # 12 is empty: it held the Fester or BHF match
                     # 15: InterProScan's Pfam matches, [start, end, index into pfam]; 16: MMseqs2's best
                     # hit, [start, end, identity %, accession], for the index page's position track
                     [[x[3], x[4], pfam_ix.setdefault(x[1], len(pfam_ix))] for x in sorted(ips.get(pid, []), key=lambda x: x[3]) if x[0] == "Pfam"] if ips else [],
                     [mm[pid][2], mm[pid][3], round(100 * mm[pid][1], 1), mm[pid][0]] if mm and mm.get(pid) else None,
                     # 17: the other whole-proteome programs with a hit, as bits (OTHER_HITS)
                     sum(1 << i for i, d in enumerate((mm, strong(fs), strong(fst), km, jh)) if d and d.get(pid)),
                     1 if r["name_status"] == "conflict" and names_conflict(r["existing_name"], r["proposed_name"]) else 0,  # 18
                     # 19: kmerseek's best hit when no other program found anything: [E-value, setting, accession, name], or 0
                     [k[1], k[0], k[2], k[3]] if k and k[1] <= ONLY_KMERSEEK_E else 0,
                     # 20: the gene groups it is in (gene_groups), or 0
                     gene_groups(r, pid, uni, aging, mm, fs, fst, km, allo) or 0,
                     # 21: Foldseek best hits weaker than FOLDSEEK_MAX_E, which field 17 leaves out but the protein's page
                     # shows: [[index into OTHER_HITS, E-value], ...], or 0
                     [[i, float(d[pid][4])] for i, d in ((1, fs), (2, fst)) if d and d.get(pid) and float(d[pid][4]) > FOLDSEEK_MAX_E] or 0])
        n_no_other += not found
        mem = [m for m in members[rep] if m != pid]
        u = uni.get(r["sp_acc"], {})
        fa = None
        if r.get("fantasia_uniprot"):
            fu = uni.get(r["fantasia_uniprot"], {})
            fa = [r["fantasia_uniprot"], fu.get("protein_name", ""), fu.get("gene", ""), fu.get("organism", ""),
                  round(float(r["fantasia_distance"]), 3), round(float(r["fantasia_ri"]), 3),
                  [g.split("|", 1) for g in split_list(r["fantasia_go_terms"]) if "|" in g]]
        card = cards.get(pid)
        shards[shard_of(pid)][pid] = dict(
            go=[[g, go.get(g, "")] for g in split_list(r["go_terms"])], kw=split_list(r["sp_keywords"]), pfam=pf,
            sp_gene=u.get("gene", ""), sp_org=u.get("organism", ""), fa=fa,
            # UniFIRE: its predictions; "not run" for proteins over UNIFIRE_MAX_AA, which it skipped
            km=None if km is None else km.get(pid, 0),
            dm=dm.get(pid),
            ss=ss.get(pid) if ss else None,  # s4pred's helix, strand and coil as runs, drawn under the ruler
            mm=None if mm is None else mm.get(pid, 0),
            fs=None if fs is None else fs.get(pid, 0),  # 0: no hit at FOLDSEEK_MAX_E or below
            fst=None if fst is None else ("not run" if lengths[pid][0] > STRUCT_MAX_AA else fst.get(pid, 0)),
            # the sequence (only for a protein without a drawn page; a drawn page carries its own), and the
            # lengths of its best hits' Swiss-Prot proteins, for the line under each hit
            q=seqs.get(pid) if not drawn_page.get(pid) else None,
            tl={a: tlen[a] for a in (*(v[0] for v in (mm and mm.get(pid), fs and fs.get(pid), fst and fst.get(pid)) if v), km and km.get(pid) and km[pid][2], r["sp_acc"]) if a in tlen},
            ips=None if ips is None else ("not run" if lengths[pid][0] > UNIFIRE_MAX_AA and not ips.get(pid) else sorted(ips.get(pid, []), key=lambda x: (x[3], x[4]))),  # 0: searched, no hit at E-value KMERSEEK_MAX_E or below
            uf=None if uf is None else ("not run" if lengths[pid][0] > UNIFIRE_MAX_AA and not uf.get(pid) else rule_places(uf.get(pid, []), sorted(ips.get(pid, []), key=lambda x: (x[3], x[4])) if ips else [], rules)),
            r=rows[-1][:15] + [0, 0, 0, rows[-1][18]],  # what the protein's pages need from its row in full/rows.js
            card={k: card[k] for k in ("card_name", "one_line", "summary", "confidence", "confidence_reason", "name_verdict", "model")} if card else None,
            sp_e=r["sp_evalue"], ev=evidence_list(im["evidence"]) if im else [],
            jh=sorted(jh.get(pid, []), key=lambda x: float(x[3])), rep=rep,
            # every protein of its gene, itself included, when the gene has more than one: [ID, length, name]
            iso=[[m, lengths[m][0], name[m]] for m in gene_proteins[pid.split(".")[0]]] if len(gene_proteins[pid.split(".")[0]]) > 1 else None,
            mem=[[m, name[m]] for m in mem[:MEMBERS_MAX]], more=max(0, len(mem) - MEMBERS_MAX))

    pfam_desc = {x[1]: x[2] for v in (ips or {}).values() for x in v if x[0] == "Pfam"}
    (out / "full").mkdir(parents=True, exist_ok=True)
    asset_ver = "?v=" + hashlib.md5(b"".join(p.read_bytes() for p in sorted((Path(__file__).parent / "assets").glob("*.*")))).hexdigest()[:8]
    if drawn:
        shutil.copytree(drawn, out, dirs_exist_ok=True)
        (out / "index.html").rename(out / "drawn.html")
        # The drawn page is the protein's page: add the Biomni tables to it (assets/viewer.js
        # fills #biomni and #biomni-more under the drawing)
        first = {}
        for pid, page in drawn_page.items():
            first.setdefault(page, pid)
        row_of = {r[0]: r for r in rows}
        for page, pid in first.items():
            f = out / page
            html = f.read_text()
            html = html.replace("</main>", '<section id="biomni" data-id="%s"></section><section id="biomni-more"></section></main>' % pid, 1)
            if html.count('id="biomni') != 2:
                sys.exit("could not add the Biomni tables to %s" % page)
            f.write_text(pack_drawn(html, drawn_additions(row_of[pid], shards[shard_of(pid)][pid], lengths[pid][0]), asset_ver))
    shutil.copytree(Path(__file__).parent / "assets", out / "assets", dirs_exist_ok=True)
    n_tseq = add_targets(out, tseqs)
    for k, entries in shards.items():
        (out / "full" / ("%03d.js" % k)).write_text('(window.FULL_GZ=window.FULL_GZ||{})[%d]="%s";' % (k, packed(entries)))
    for k in range(N_SHARDS):
        if k not in shards:
            (out / "full" / ("%03d.js" % k)).write_text('(window.FULL_GZ=window.FULL_GZ||{})[%d]="%s";' % (k, packed({})))
    # isoforms: a gene's first protein has the bare ID, the others end in .2, .3, ... (make_proteome_lengths.py)
    n_iso = Counter(r[0].split(".")[0] for r in rows)
    iso = {g: k for g, k in n_iso.items() if k > 1}
    full = dict(cats=cats, check=[t for _, t in NAME_CHECK], other=OTHER_HITS, total=len(rows), rows=rows, groups=[g for g, _ in GENE_GROUPS], iso=iso,
                pfam=[[a, pfam_desc.get(a, "")] for a in sorted(pfam_ix, key=pfam_ix.get)])
    (out / "full" / "rows.js").write_text('window.FULL_Z="%s";' % packed(full))
    (out / "full" / "meta.js").write_text("window.FULL_META=%s;" % json.dumps(
        dict(cats=cats, check=full["check"], check_differs=status_code["conflict"], verdicts=VERDICT_TEXT), separators=(",", ":")))
    # What index.html needs to show its first screen before full/rows.js arrives
    used = {x[2] for r in rows[:100] for x in r[15]}
    first = dict(full, rows=rows[:100], pfam={i: full["pfam"][i] for i in used}, iso={g: k for g, k in iso.items() if g in {r[0].split(".")[0] for r in rows[:100]}},
                 ladder=[sum(1 for r in rows if r[6] and r[9] and float(r[9]) <= c) for c in (1e-3, 1e-10, 1e-50, 1e-100)])

    n = len(rows)
    summ = dict(
        total=n,
        swissprot_hit=sum(1 for r in rows if r[6]),
        pfam_domain=sum(1 for r in rows if r[10]),
        neither=sum(1 for r in rows if not r[6] and not r[10]),
        immune_candidates=sum(1 for r in rows if r[11]),
        jackhmmer_hit=sum(1 for pid in jh if pid in lengths),
        jackhmmer_not_immune=sum(1 for pid in jh if pid not in immune),
        in_cluster_of_2_or_more=sum(1 for r in rows if r[13] > 1),
        clusters=len(members),
        foldseek_hit=None if fs is None else len(strong(fs)),
        foldseek_struct_hit=None if fst is None else len(strong(fst)),
        drawn=sum(1 for r in rows if r[14]),
        drawn_pages=len({r[14] for r in rows if r[14]}),
        names_conflict=sum(1 for r in rows if r[18]),
        no_other_program=n_no_other,
        only_kmerseek_3=sum(1 for r in rows if r[19]),
        s4pred=ss_summary,  # None when s4pred_states.tsv is not there
        only_kmerseek_5=sum(1 for r in rows if r[19] and r[19][0] <= 1e-5),
        name_check={t: sum(1 for r in rows if r[4] == i) for i, (_, t) in enumerate(NAME_CHECK)},
        immune_categories={c: sum(1 for r in rows if r[11] >> i & 1) for i, c in enumerate(cats)},
        genes=len(n_iso), genes_with_isoforms=len(iso), most_isoforms=max(n_iso.values()),
        genes_by_isoforms=dict(sorted(Counter(n_iso.values()).items())),
        gene_groups={g: sum(1 for r in rows if r[20] and any(x[0] == i for x in r[20])) for i, (g, _) in enumerate(GENE_GROUPS)},
    )
    (out / "full_summary.json").write_text(json.dumps(summ, indent=1))
    inputs = sorted(p.name for p in folder.iterdir() if p.suffix in (".tsv", ".gz") and p.name not in NOT_READ)
    (out / "index.html").write_text(index_html(summ, cats, inputs, asset_ver, known_html(drawn_page), first))
    (out / "protein.html").write_text(protein_html(asset_ver, n_tseq))
    print("%d proteins, %d drawn, %d data files, %.1f s" % (n, summ["drawn"], N_SHARDS, time.time() - t0))
    return summ


def head(title, ver, dense=False):
    """The page head; dense: the first page's stylesheets (dense-bhf.css, a copy of the BHF page's dense.css; index.css) after site.css."""
    more = '<link rel="stylesheet" href="assets/dense-bhf.css%s"><link rel="stylesheet" href="assets/index.css%s">' % (ver, ver) if dense else ""
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1"><title>%s</title>%s'
            '<link rel="stylesheet" href="assets/site.css%s">%s</head><body>' % (build.esc(title), build.CSS_LINKS, ver, more))


def bars(items, total, key):
    """Counts as horizontal bars of one ink colour; each line filters the table when clicked."""
    top = max(v for _, v, _ in items) or 1
    rows = "".join(
        '<button class="bar" data-%s="%s"><span class="blab">%s</span><span class="bt"><i style="width:%.1f%%"></i></span>'
        '<span class="bn">%s</span></button>' % (key, k, build.esc(t), 100 * v / top, "{:,}".format(v))
        for k, v, t in items)
    return '<div class="bars">%s</div>' % rows


# Proteins the lab already knows, linked at the top of the index when they have a drawn page:
# (protein ID, short name, full name)
KNOWN = [("FUN008084", "BHF", "Botryllus histocompatibility factor")]


def known_html(drawn_page):
    links = ['<a class="known" href="%s"><b>%s</b> %s <span class="mono">%s</span></a>'
             % (build.esc(drawn_page[pid]), build.esc(short), build.esc(full), build.esc(pid))
             for pid, short, full in KNOWN if pid in drawn_page]
    return '<span class="knownrow">go straight to %s</span>' % " ".join(links) if links else ""


AXIS_AA = 1000  # the index page's position track runs from residue 1 to this; 95% of proteins are shorter
TRACK_PCT = 88  # the share of the track the scale takes; the rest holds a longer protein's "› 16,045 aa" (full.js TW)
# The other whole-proteome programs, in the order of row field 17's bits, with the cut each hit passed
OTHER_HITS = [["MMseqs2", "E-value 0.001 or below"], ["ProstT5 + Foldseek", "E-value 0.001 or below"],
              ["Foldseek on structure", "E-value 0.001 or below, proteins of 200 aa or less"],
              ["kmerseek", "region E-value %g or below" % KMERSEEK_MAX_E], ["jackhmmer, immune seeds", "any hit in Biomni's table"]]


# The paper of each program and database the pages name, for the index page's reference list. Every
# entry was read from its PubMed record on 2026-11-10 (authors, journal, volume, pages, DOI); none is
# written from memory. (what it is on the pages, authors, title, journal, volume and issue, year, pages, DOI)
CITATIONS = [
    ("DIAMOND", "Buchfink, Benjamin, Klaus Reuter, and Hajk-Georg Drost", "Sensitive protein alignments at tree-of-life scale using DIAMOND",
     "Nature Methods", "18, no. 4", 2021, "366–68", "10.1038/s41592-021-01101-x"),
    ("MMseqs2", "Steinegger, Martin, and Johannes Söding", "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets",
     "Nature Biotechnology", "35, no. 11", 2017, "1026–28", "10.1038/nbt.3988"),
    ("Foldseek", "van Kempen, Michel, et al.", "Fast and accurate protein structure search with Foldseek",
     "Nature Biotechnology", "42, no. 2", 2024, "243–46", "10.1038/s41587-023-01773-0"),
    ("ProstT5", "Heinzinger, Michael, et al.", "Bilingual language model for protein sequence and structure",
     "NAR Genomics and Bioinformatics", "6, no. 4", 2024, "lqae150", "10.1093/nargab/lqae150"),
    ("HMMER: hmmsearch, hmmscan, phmmer, jackhmmer", "Eddy, Sean R.", "Accelerated profile HMM searches",
     "PLoS Computational Biology", "7, no. 10", 2011, "e1002195", "10.1371/journal.pcbi.1002195"),
    ("pyhmmer", "Larralde, Martin, and Georg Zeller", "PyHMMER: a Python library binding to HMMER for efficient sequence analysis",
     "Bioinformatics", "39, no. 5", 2023, "btad214", "10.1093/bioinformatics/btad214"),
    ("HH-suite (hhsearch)", "Steinegger, Martin, et al.", "HH-suite3 for fast remote homology detection and deep protein annotation",
     "BMC Bioinformatics", "20, no. 1", 2019, "473", "10.1186/s12859-019-3019-7"),
    ("InterProScan", "Jones, Philip, et al.", "InterProScan 5: genome-scale protein function classification",
     "Bioinformatics", "30, no. 9", 2014, "1236–40", "10.1093/bioinformatics/btu031"),
    ("Pfam", "Paysan-Lafosse, Typhaine, et al.", "The Pfam protein families database: embracing AI/ML",
     "Nucleic Acids Research", "53, no. D1", 2025, "D523–34", "10.1093/nar/gkae997"),
    ("UniFIRE's UniRule rules", "MacDougall, Alistair, et al.", "UniRule: a unified rule resource for automatic annotation in the UniProt Knowledgebase",
     "Bioinformatics", "36, no. 17", 2020, "4643–48", "10.1093/bioinformatics/btaa485"),
    ("FANTASIA, the method across about 1,000 animal proteomes", "Martínez-Redondo, Gemma I., et al.", "FANTASIA leverages language models to decode the functional dark proteome across the animal tree of life",
     "Communications Biology", "8, no. 1", 2025, "1227", "10.1038/s42003-025-08651-2"),
    ("FANTASIA, the software", "Pérez-Canales, Francisco M., et al.", "FANTASIA suite: a reproducible and configurable framework for embedding-based functional annotation of proteins",
     "NAR Genomics and Bioinformatics", "8, no. 3", 2026, "lqag106", "10.1093/nargab/lqag106"),
    ("FANTASIA, the earlier test on model organisms", "Barrios-Núñez, Israel, et al.", "Decoding functional proteome information in model organisms using protein language models",
     "NAR Genomics and Bioinformatics", "6, no. 3", 2024, "lqae078", "10.1093/nargab/lqae078"),
    ("ESMFold, the model before ESMFold2", "Lin, Zeming, et al.", "Evolutionary-scale prediction of atomic-level protein structure with a language model",
     "Science", "379, no. 6637", 2023, "1123–30", "10.1126/science.ade2574"),
    ("Swiss-Prot (UniProt)", "The UniProt Consortium", "UniProt: the Universal Protein Knowledgebase in 2025",
     "Nucleic Acids Research", "53, no. D1", 2025, "D609–17", "10.1093/nar/gkae1010"),
    ("AlphaFold's models of Swiss-Prot, which Foldseek searched", "Bertoni, Damian, et al.", "AlphaFold Protein Structure Database 2025: a redesigned interface and updated structural coverage",
     "Nucleic Acids Research", "54, no. D1", 2026, "D358–62", "10.1093/nar/gkaf1226"),
    ("Mol*, the 3-D viewer", "Sehnal, David, et al.", "Mol* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures",
     "Nucleic Acids Research", "49, no. W1", 2021, "W431–37", "10.1093/nar/gkab314"),
    ("GenAge, the aging gene group", "Tacutu, Robi, et al.", "Human Ageing Genomic Resources: new and updated databases",
     "Nucleic Acids Research", "46, no. D1", 2018, "D1083–90", "10.1093/nar/gkx1042"),
    ("Biomni, the agent that ran the whole-proteome programs", "Huang, Kexin, et al.", "Biomni: a general-purpose biomedical AI agent",
     "bioRxiv", "preprint", 2025, "", "10.1101/2025.05.30.656746"),
    ("s4pred, the helix and strand prediction", "Moffat, Lewis, and David T. Jones", "Increasing the accuracy of single sequence prediction methods using a deep semi-supervised learning framework",
     "Bioinformatics", "37, no. 21", 2021, "3744–51", "10.1093/bioinformatics/btab491"),
]
# Named on the pages with no paper to cite
NO_PAPER = [("kmerseek", 'no paper yet; the code is at <a href="https://github.com/seanome/kmerseek">github.com/seanome/kmerseek</a>'),
            ("ESMFold2", "we found no PubMed record for this model; the paper of the first ESMFold is listed above"),
            ("UniFIRE's ARBA rules", "we found no PubMed record for ARBA; UniProt describes it in its online help")]


def citations_html():
    """The reference list of the index page, Chicago style, one line per program or database."""
    e = build.esc
    items = "".join(
        '<li><b>%s.</b> %s. "%s." <i>%s</i> %s (%d)%s. <a href="https://doi.org/%s">https://doi.org/%s</a></li>'
        % (e(what), e(authors).rstrip("."), e(title), e(journal), e(vol), year, ": " + e(pages) if pages else "", doi, doi)
        for what, authors, title, journal, vol, year, pages, doi in CITATIONS)
    none = "".join("<li><b>%s.</b> %s.</li>" % (e(what), why[0].upper() + why[1:]) for what, why in NO_PAPER)
    return '<ol class="refs">%s</ol><p class="note">No paper to cite:</p><ul class="refs">%s</ul>' % (items, none)


def index_html(s, cats, inputs, ver, known="", first=None):
    """The first page, after the dense-results-page skill: the answer in one sentence, then a sticky
    header (filters, search, a three-line key, column names, the residue ruler) over one grid row per
    protein, each with its Pfam matches and MMseqs2's best hit drawn along the protein."""
    f = "{:,}".format
    n = s["total"]
    pc = lambda v: "%d%%" % round(100 * v / n)
    chips = [("all", n, "all"), ("sp", s["swissprot_hit"], "Swiss-Prot match"), ("pfam", s["pfam_domain"], "Pfam domain"),
             ("none", s["neither"], "neither"), ("immune", s["immune_candidates"], "immune-gene candidate"),
             ("drawn", s["drawn"], "hits along the protein"), ("names", s["names_conflict"], "names conflict"),
             ("km3", s["only_kmerseek_3"], "only kmerseek, E ≤ 0.001"), ("km5", s["only_kmerseek_5"], "only kmerseek, E ≤ 0.00001")]
    only_tip = ("no other program found anything on this protein (DIAMOND, Pfam, InterProScan apart from its disorder, coil, membrane and signal "
                "predictions, MMseqs2, Foldseek, jackhmmer, FANTASIA, UniFIRE), and kmerseek's best hit over 35 settings passes this E-value. "
                "Of the %s proteins nothing else found, about %s would pass by chance if the E-values are right and the settings are independent")
    tips = {"all": "every protein", "sp": "DIAMOND found a Swiss-Prot match at the cutoff", "pfam": "hmmsearch found a Pfam domain",
            "none": "neither a Swiss-Prot match at the cutoff nor a Pfam domain; another program may still have a hit",
            "immune": "in at least one immune-gene category", "drawn": "has a page that shows every program's hits along the protein, as one of our test proteins",
            "names": "the FASTA header name and the proposed name are different names; capitals, punctuation and a trailing -like are ignored",
            "km3": only_tip % (f(s["no_other_program"]), f(round(s["no_other_program"] * 35 * 1e-3))),
            "km5": only_tip % (f(s["no_other_program"]), f(round(s["no_other_program"] * 35 * 1e-5)))}
    chip_html = "".join('<button class="chip%s" data-t="%s" title="%s">%s <b class="mono">%s</b></button>' % (" on" if k == "all" else "", k, tips[k], build.esc(w), f(v))
                        for k, v, w in chips)
    cat_opts = "".join('<option value="%d">%s (%s)</option>' % (i, build.esc(c), f(s["immune_categories"][c])) for i, c in enumerate(cats))
    chk_opts = "".join('<option value="%d">%s (%s)</option>' % (i, build.esc(t), f(s["name_check"][t])) for i, (_, t) in enumerate(NAME_CHECK))
    grp_opts = "".join('<option value="%d">%s (%s)</option>' % (i, build.esc(g), f(s["gene_groups"][g])) for i, (g, _) in enumerate(GENE_GROUPS))
    cat_bars = bars([(i, s["immune_categories"][c], c) for i, c in enumerate(cats)], n, "cat")
    chk_bars = bars([(i, s["name_check"][t], t) for i, (_, t) in enumerate(NAME_CHECK)], n, "chk")
    iso_bars = bars([(k, v, "1 protein, bare ID" if int(k) == 1 else "%s proteins, up to .%s" % (k, k)) for k, v in s["genes_by_isoforms"].items()], n, "iso")
    seq, dom = build.CATEGORIES["seq"]["color"], build.CATEGORIES["pfam"]["color"]
    ramp = "".join('<i style="opacity:%s;background:%s"></i>' % (o, seq) for o in (".3", ".5", ".75", "1"))
    # the same scale as full.js's track: residue 1 to AXIS_AA over the first TRACK_PCT% of the width
    at = lambda v: TRACK_PCT * (v - 1) / (AXIS_AA - 1)
    ticks = "".join('<text x="%.2f%%" y="10" class="ax" text-anchor="%s">%s</text><line x1="%.2f%%" x2="%.2f%%" y1="12" y2="16" stroke="#c3c2b7"/>'
                    % (at(v), "start" if v == 1 else "middle", f(v), at(v), at(v))
                    for v in (1, AXIS_AA // 4, AXIS_AA // 2, 3 * AXIS_AA // 4, AXIS_AA))
    return head("Botryllus proteome annotations", ver, dense=True) + """
<div class="titleblock">
 <div class="titlebar"><h1><i>Botryllus schlosseri</i> proteome</h1><span class="pn">%(n)s proteins predicted from the new genome: what each one matches</span>
 <span class="pn" title="A gene's first protein has the bare ID; its other isoforms end in .2, .3 and so on, in the order of the proteome FASTA. Each row says which isoform it is.">%(genes)s genes · %(gi)s have 2 to %(mi)s isoforms</span>%(known)s</div>
 <div class="sentence"><p id="sent">%(sp)s of %(n)s proteins (%(spp)s) match a Swiss-Prot protein, %(pf)s (%(pfp)s) carry a Pfam domain, and %(ne)s (%(nep)s) have neither.</p>
 <div id="sentnote"><p class="note">Swiss-Prot match = the best hit of DIAMOND (a fast BLAST-like sequence search) at E-value 0.001 or below,
 in Swiss-Prot release 2026_03 (575,748 proteins). Pfam domain = a Pfam family model that pyhmmer's hmmsearch finds above the score Pfam sets
 for that family. Biomni, the agent that ran the whole-proteome programs except kmerseek, made both tables. A protein with neither can still have a hit
 from another program: the "other programs" column says which. %(dp)s test proteins have a page that shows every program's hits along the protein
 (<a href="drawn.html">their own index</a>).</p></div></div>
 <details class="overview-fold"><summary>Counts by immune-gene category, by how the proposed name compares with the existing one, and by isoforms per gene</summary>
 <div class="overview">
 <section><h2>Immune-gene candidates, by category</h2><p class="note">%(im)s proteins; one protein can be in several categories. A category
 comes from a Pfam domain, the gene name, or an immune keyword or GO term of the best Swiss-Prot match; each protein page says which.
 A jackhmmer hit from an immune seed protein does not give a category: %(jn)s proteins have such a hit and no category. Bar length: share of
 the %(im)s; click a bar to list its proteins.</p>%(cb)s</section>
 <section><h2>Proposed name compared with the existing name</h2><p class="note">Existing = the name in the FASTA header. Proposed = from
 the best Swiss-Prot match or the Pfam domains. "Differs" includes a proposed name that only adds "-like" (NLRC5 and Nlrc5-like). Bars: of all
 %(n)s proteins.</p>%(kb)s</section>
 <section><h2>Genes by number of proteins (isoforms)</h2><p class="note">%(genes)s genes. A gene's first protein has the bare ID; its other
 isoforms end in .2, .3 and so on, in the order of the proteome FASTA. %(gi)s genes have more than one protein, which adds %(extra)s proteins.
 Bar length: number of genes; click a bar to list those genes' proteins.</p>%(ib)s</section>
 </div></details>
</div>
<header class="top ix">
 <div class="filters">
  <span class="lead2">DIAMOND E-value cutoff:</span><span class="ladder" id="eladder" title="Which DIAMOND hits count as a Swiss-Prot match. Click one."></span>
  <span class="lead2">show:</span><span class="ladder" id="tiles">%(chips)s</span>
  <label>immune category <select id="cat"><option value="-1">any</option>%(cat_opts)s</select></label>
  <label>names <select id="chk"><option value="-1">any</option>%(chk_opts)s</select></label>
  <label title="Proteins whose best hit from DIAMOND, MMseqs2, Foldseek or kmerseek is a Swiss-Prot protein of this group, at any E-value, weak hits included: the last column then names the program, its E-value and the hit. Aging: the hit's gene symbol is one of the 307 human genes of GenAge build 21.">gene group <select id="grp"><option value="-1">any</option>%(grp_opts)s</select></label>
  <span class="search"><input id="q" type="search" placeholder="search: every word must match" autocomplete="off" spellcheck="false" aria-label="Search the ID, names, Swiss-Prot match, Pfam domains and categories: every word must match"></span>
  <button class="btn" id="dl" title="Every protein the filters keep, not only the rows shown, one line each, as a tab-separated file">&#8595; TSV of these proteins</button>
  <button class="btn" id="help-btn" aria-expanded="false" title="What the columns and marks mean; keys">? Help</button>
  <button class="btn" id="clear" title="Show all proteins again: cutoff 0.001, no filter, no search, no sort">&#215; clear all</button>
  <span class="count" id="count" aria-live="polite"></span>
 </div>
 <div class="key"><span class="lead">along each protein:</span>
  <span><svg width="26" height="10"><line x1="0" x2="26" y1="5" y2="5" stroke="#c8c8c8" stroke-width="2"/></svg>the protein on one scale, residue 1 to %(ax)s; a longer one ends in &#8250; and its length</span>
  <span><svg width="16" height="10"><rect x="1" y="1" width="14" height="8" fill="#fff" stroke="%(dom)s" stroke-width="1.2"/></svg>Pfam match by InterProScan (can differ from the Pfam column)</span>
  <span class="kind"><b>bar: MMseqs2's best Swiss-Prot hit</b>, whose positions are known (DIAMOND's are not; the row names MMseqs2's hit when it is another protein); fainter = fewer identical residues <i class="ramp">%(ramp)s</i> &lt;20%% · 20–40%% · 40–60%% · 60%%+</span>
 </div>
 <div class="g ixg colhead">
  <div><button data-k="0" data-t="s">protein · existing name &#8594; proposed <i></i></button></div>
  <div class="num"><button data-k="1">length (aa) <i></i></button></div>
  <div><button data-k="5" data-t="s" title="the best hit of DIAMOND against Swiss-Prot 2026_03">Swiss-Prot match (DIAMOND) <i></i></button></div>
  <div class="num"><button data-k="7" title="identity: identical residues in DIAMOND's alignment; coverage: the share of your protein it spans">identity · coverage <i></i></button></div>
  <div class="num"><button data-k="9" data-t="e">E-value <i></i></button></div>
  <div><button data-k="10" data-t="s" title="Pfam models found by pyhmmer's hmmsearch (Biomni); under it, InterProScan's Pfam matches">Pfam domains <i></i></button></div>
  <div>immune categories</div>
  <div><button data-k="17" title="other whole-proteome programs with a hit on this protein; Help lists each one's cutoff">other programs with a hit <i></i></button></div>
  <div class="num"><button data-k="13" title="proteins in its MMseqs2 cluster (Biomni's clustering of the proteome), itself included">cluster <i></i></button></div>
  <div class="track">position along the protein (residue)</div>
 </div>
 <div class="g ixg ruler"><div class="lbl">%(cnt)s</div><div class="track"><svg width="100%%" height="16">%(ticks)s</svg></div></div>
</header>
<main><div id="grid" role="list" data-ax="%(axn)d" data-tw="%(tw)d" data-seq="%(seq)s" data-dom="%(dom)s"></div><p class="morebar"><button class="btn" id="more" hidden>Show 100 more</button></p></main>
<div class="help" id="help" role="dialog" aria-label="Help">
 <h2>Rows</h2><p>One row per protein, 100 at a time. Click a row, or press Enter on it, to open its page.</p>
 <h2>Columns</h2><dl>
 <dt>protein</dt><dd>the ID in the proteome FASTA (gene ID, underscore, transcript ID; a second isoform of the same gene gets .2), then the
 name in the FASTA header and the proposed name. "isoform 2 of 3" = the gene has three proteins and this is the second in the FASTA; a
 protein that is its gene's only one says nothing</dd>
 <dt>best Swiss-Prot match</dt><dd>DIAMOND's best hit (E-value 0.001 or below) in Swiss-Prot 2026_03, with its accession</dd>
 <dt>identity · coverage</dt><dd>identical residues in DIAMOND's alignment, and the share of your protein the alignment spans</dd>
 <dt>Pfam domains</dt><dd>Pfam models found by pyhmmer's hmmsearch at Pfam's gathering thresholds (Biomni). The boxes in the track are
 InterProScan's Pfam matches, which can differ</dd>
 <dt>cluster</dt><dd>proteins in the same MMseqs2 cluster (Biomni's clustering of the proteome), itself included</dd>
 <dt>other programs with a hit</dt><dd>%(other)s</dd>
 <dt>gene group</dt><dd>a filter, not a finding: it keeps a protein when the best hit of DIAMOND, MMseqs2, either Foldseek search or kmerseek
 is a Swiss-Prot protein whose name is in the group (%(grp)s), at any E-value, so weak hits from one program are kept. While it is on, the
 "other programs" column names the program, its E-value and the hit, best first. MHC class I and class II are the two halves of
 "HLA and other MHC". Aging keeps a protein when the gene symbol of such a hit is one of the 307 human genes of GenAge build 21
 (genomics.senescence.info/genes, Tacutu et al. 2018, CC BY 3.0), capitals ignored, in whichever species the hit is; a hit with no gene
 symbol in our names table cannot count. Allorecognition is the one group that is not a match to another species: it keeps the %(allo_n)s
 proteins that DIAMOND matched to one of 96 <i>Botryllus schlosseri</i> proteins from NCBI (Fester family, Fester co-receptor and one partial
 BHF), so a hit can be the protein itself; the column names the best of those</dd>
 <dt>position along the protein</dt><dd>the protein from residue 1 to %(ax)s on one scale for every row; a longer protein ends in &#8250; and
 its length</dd></dl>
 <h2>Programs</h2><dl><dt>DIAMOND, MMseqs2</dt><dd>fast sequence searches, like BLAST, against Swiss-Prot</dd>
 <dt>hmmsearch, InterProScan</dt><dd>find Pfam family models (profiles of a domain's sequence) in the protein</dd>
 <dt>jackhmmer</dt><dd>repeats a profile search to reach distant relatives; here run from immune-gene seed proteins against this proteome</dd>
 <dt>ProstT5 + Foldseek</dt><dd>predicts a 3-D shape alphabet from the sequence and searches AlphaFold's models of Swiss-Prot</dd>
 <dt>Foldseek on structure</dt><dd>the same search from a structure ESMFold2 predicted, for proteins of 200 aa or less</dd>
 <dt>kmerseek</dt><dd>finds shared runs of residues in a reduced amino-acid alphabet</dd></dl>
 <h2>Sorting</h2><p>Click a column name with &#8597; to sort; the first click puts the best first (longest, most identical, smallest E-value,
 biggest cluster, A to Z).</p>
 <h2>Keys</h2><p><kbd>j</kbd>/<kbd>k</kbd> next / previous row · <kbd>Enter</kbd> or <kbd>Space</kbd> open its page · <kbd>/</kbd> search ·
 <kbd>Esc</kbd> close this, clear the search, or leave the row</p>
</div>
<section class="below"><h2>Which Swiss-Prot</h2><p class="note">"Swiss-Prot" on this page is release 2026_03 (575,748 proteins), the copy
DIAMOND searched; 533,805 of its proteins are among DIAMOND's hits (at most 25 per protein). Each protein's page says which copy each other
program searched.</p></section>
<section class="below"><h2>Papers of the programs and databases</h2><p class="note">One paper for each program and database these pages
name. Each entry was read from its PubMed record on 2026-11-10.</p>%(refs)s</section>
<footer><h2>Input files read</h2><p class="note">Built by build_full.py from %(inputs)s. These tables come from Biomni's runs and our kmerseek run on Sherlock
(full_inputs/README.md in the repository says how each was made). This page only reads those tables; building it ran no search.</p></footer>
<script>window.FULL_FIRST=%(first)s;</script><script src="assets/gz.js%(ver)s"></script><script src="assets/full.js%(ver)s"></script></body></html>
""" % dict(first=dumps(first), n=f(n), known=known, sp=f(s["swissprot_hit"]), spp=pc(s["swissprot_hit"]), pf=f(s["pfam_domain"]), pfp=pc(s["pfam_domain"]),
           ne=f(s["neither"]), nep=pc(s["neither"]), dp=f(s["drawn_pages"]), dr=f(s["drawn"]), im=f(s["immune_candidates"]),
           jn=f(s["jackhmmer_not_immune"]), cb=cat_bars, kb=chk_bars, chips=chip_html, cat_opts=cat_opts, chk_opts=chk_opts, grp_opts=grp_opts,
           refs=citations_html(), genes=f(s["genes"]), gi=f(s["genes_with_isoforms"]), mi=s["most_isoforms"], ib=iso_bars, extra=f(n - s["genes"]),
           allo_n=f(next(v for g, v in s["gene_groups"].items() if g.startswith(ALLO))),
           grp=build.esc("; ".join(g for g, p in GENE_GROUPS if p not in (GENAGE, ALLO))),
           other="; ".join("%s: %s" % (a, c) for a, c in OTHER_HITS), ax=f(AXIS_AA), axn=AXIS_AA, tw=TRACK_PCT, seq=seq, dom=dom, ramp=ramp, ticks=ticks, cnt="", inputs=build.esc(", ".join(inputs)), ver=ver)


def protein_html(ver, n_tseq=0):
    """protein.html: the dense page (assets/dense-bhf.js, as on a drawn page) for any protein, from
    the whole-proteome tables; assets/protein-dense.js builds its data from full/NNN.js. The key lists
    the programs with positions in those tables; the Swiss-Prot table lists the programs run on every protein."""
    run = {t["key"]: ("run" if t["key"] in ("mmseqs2", "kmerseek", "prostt5", "foldseek", "interproscan") else "not run") for t in build.TRACKS}
    ctx = dict(run_status=run)
    whole = [["w_diamond"], ["w_jackhmmer"]]
    const = dict(
        tracks={t["key"]: dict(name=t["name"], color=build.PROGRAMS[t["program"]]["color"], shape=build.PROGRAMS[t["program"]]["shape"],
                               cat=build.PROGRAMS[t["program"]]["cat"], p=t["program"]) for t in build.TRACKS},
        program_names={k: v["name"] for k, v in build.PROGRAMS.items()}, ips_drawn=build.IPS_DRAWN, ips_evalue=sorted(build.IPS_EVALUE),
        fantasia_max=build.FANTASIA_MAX_DIST, kmerseek_max_e=KMERSEEK_MAX_E, struct_max_aa=STRUCT_MAX_AA,
        pfam_color=build.CATEGORIES["pfam"]["color"], ver=ver, tseq_files=n_tseq, alpha=build.parse_alphabets(Path(__file__).parent / "test_inputs" / "data" / "kmerseek_alphabets.rs"))
    for k, t in (("diamond", "DIAMOND, best hit"), ("unifire", "UniFIRE (UniProt's rules)")):
        const["tracks"][k] = dict(name=t, color=build.WHOLE_TRACKS[k]["color"], shape=build.WHOLE_TRACKS[k]["shape"], cat=build.WHOLE_TRACKS[k]["cat"], p=k)
    const["program_names"].update(diamond="DIAMOND", unifire="UniFIRE")
    const["tracks"]["kmerseek"]["name"] = "kmerseek, best hit"
    head = build.page_head("Botryllus protein: what every annotation program found", 0).replace(
        "</head>", '<link rel="stylesheet" href="assets/dense-bhf.css%s"></head>' % ver, 1)
    help_ = ('<aside class="help" id="help"><h2>What this page is</h2><p>Every program\'s best hit on this protein, from Biomni\'s runs on '
             'the whole proteome (Biomni is the agent that ran them, except kmerseek, which ran on Sherlock). MMseqs2, ProstT5 + Foldseek and Foldseek on the '
             'predicted structure are placed where they lie; InterProScan places every match; kmerseek, DIAMOND, jackhmmer with immune seeds, '
             'FANTASIA and UniFIRE name a hit with no position. The 51 test proteins also have phmmer, jackhmmer against '
             'Swiss-Prot, hmmscan, hhsearch and MMseqs2 three rounds, every hit, and the residues of each.</p><h2>Rows</h2><p>A <b>program '
             'row</b> folds its hits into one row; click it for one row per hit. A hit row gives the residues of your protein it covers, its '
             'E-value and score, the identity where the program gives one, and the matched protein; click it for its numbers.</p>'
             '<h2>Controls</h2><p>The <b>E-value cutoff</b> chips keep the hits at or below that value; the sentence and every count follow '
             'it. The best hits of MMseqs2 and the Foldseek runs are all at 0.001 or below, since the tables keep only those. <b>Search</b> '
             'looks through the matched protein\'s name, accession and species; every word typed must appear.</p><h2>Keyboard</h2><dl>'
             '<dt><kbd>j</kbd> <kbd>k</kbd></dt><dd>next and previous row</dd><dt><kbd>Enter</kbd> <kbd>Space</kbd></dt><dd>unfold a program, '
             'or open a hit</dd><dt><kbd>Esc</kbd></dt><dd>close every open hit and this help</dd><dt><kbd>?</kbd></dt><dd>this help</dd></dl></aside>')
    return head + """<div class="titleblock"><div class="titlebar"><a class="back" href="index.html" title="back to the table of all proteins"><span aria-hidden="true">←</span> All proteins</a><h1 title="the ID as written in the proteome FASTA">loading…</h1><span class="pn"></span><span class="len"></span></div><div class="sentence"><p id="sent"></p><div id="sentnote"></div></div></div>
<header class="top">%s%s%s</header>
<main><div class="withviewer noviewer"><div id="grid"></div></div>
<section class="below" id="biomni" data-base="" data-undrawn="1"></section><section class="below" id="biomni-more"></section>%s</main>
<footer><h2>Input files read</h2><p class="note">Built by build_full.py from tables of the whole proteome, Biomni's and our kmerseek run's (full_inputs/README.md in the
repository says how each was made). This page only reads them; building it ran no search.</p></footer>
<div id="tip" role="tooltip" hidden></div>%s
<script>window.DENSE_CONST=%s;</script><script src="assets/gz.js%s"></script><script src="full/meta.js%s"></script><script src="assets/protein-dense.js%s"></script></body></html>
""" % (build.DENSE_FILTERS % (round(build.math.log10(build.SHOW_E)), round(build.math.log10(build.COUNT_E))),
       build.dense_key(ctx, False, False, (("diamond", "DIAMOND"), ("unifire", "UniFIRE"))), build.DENSE_COLHEAD, build.dense_swissprot(ctx, whole), help_,
       json.dumps(const, separators=(",", ":")), ver, ver, ver)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--drawn", type=Path, help="a site made by build.py, to link to and copy in")
    ap.add_argument("--drawn-fasta", type=Path, help="the proteins.fasta that site was built from")
    ap.add_argument("--cards", type=Path, nargs="*", help="gene cards from gene_cards/run_cards.py (.jsonl)")
    ap.add_argument("--proteome-fasta", type=Path, help="the proteome FASTA (not in the repository: the proteome is unpublished); "
                    "each protein's sequence goes into its data file, so protein.html can show its residues")
    ap.add_argument("--swissprot-fasta", type=Path, help="uniprot_sprot.fasta(.gz) of the release the best hits searched (2026_03); "
                    "their Swiss-Prot proteins go into the target files, so protein.html can align the two")
    a = ap.parse_args()
    if bool(a.drawn) != bool(a.drawn_fasta):
        sys.exit("--drawn and --drawn-fasta go together")
    if a.out.exists() and any(a.out.iterdir()) and not (a.out / "full_summary.json").exists():
        print("not replacing %s: it is not empty and holds no full_summary.json, so build_full.py did not make it" % a.out)
        sys.exit(2)
    if a.out.exists():
        shutil.rmtree(a.out)
    a.out.mkdir(parents=True)
    (a.out / "full_summary.json").write_text('{"finished": false}\n')
    build_site(a.inputs, a.out, a.drawn, a.drawn_fasta, a.cards, a.proteome_fasta, a.swissprot_fasta)


if __name__ == "__main__":
    main()
