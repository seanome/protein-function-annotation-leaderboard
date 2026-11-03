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
    full/rows.js        what the index shows for each protein
    full/NNN.js         the rest of each protein's data, in 96 files; a protein's file number is a
                        hash of its ID (shard_of), so protein.html loads one file
    full_summary.json   the counts on the index
    drawn.html, proteins/, targets/, assets/, search.js  copied from --drawn; drawn.html is
                        build.py's index of the proteins drawn there
"""
import argparse
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
# species are not shown, so those two files are not read, and the immune category those hits set
# is dropped from the immune table.
NOT_READ = {"botryllus_allorecognition_hits.tsv", "allorecognition_references.tsv"}
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
            if float(r["evalue"]) <= FOLDSEEK_MAX_E:
                acc = build.af_accession(r["target"])  # AF-P12345-F1-model_v6, AF-C3VPR6-2-F1-model_v6
                out[r["protein_id"]] = [acc, float(r["fident"]), int(r["qstart"]), int(r["qend"]), r["evalue"], r["bits"],
                                        uni.get(acc, {}).get("protein_name", ""), uni.get(acc, {}).get("organism", "")]
                if "lddt" in r:  # Foldseek on structures: TM-score and LDDT of the alignment
                    out[r["protein_id"]] += [float(r["alntmscore"]), float(r["lddt"])]
                out[r["protein_id"]] += [int(r["tstart"]), int(r["tend"])]  # last two: the range on the other protein
        return out
    fs = foldseek("foldseek_best_per_protein.tsv")
    fst = foldseek("foldseek_struct_best_per_protein.tsv")
    return ann, lengths, immune, jh, rep_of, members, uni, go, uf, km, ips, mm, fs, fst


# The text that repeats across the data files (GO term names and IDs, UniProt keywords, InterProScan's
# member databases, signatures and descriptions, UniFIRE's rules and predictions) is stored once, in
# full/str.js, and the data files hold its number: this keeps the site under the host's 256 MiB. A
# data file calls FULL_X (defined in full/str.js) to turn the numbers back into text as it loads.
STR_X_JS = """window.FULL_X=function(d){var S=window.FULL_STR,s=function(i){return typeof i==="number"?S[i]:i;},p=function(g){return [s(g[0]),s(g[1])];};
for(var id in d){var e=d[id];if(e.go)e.go=e.go.map(p);if(e.kw)e.kw=e.kw.map(s);
if(Array.isArray(e.ips))e.ips=e.ips.map(function(m){return [s(m[0]),s(m[1]),s(m[2]),m[3],m[4],m[5],s(m[6]),s(m[7])];});
if(Array.isArray(e.uf))e.uf=e.uf.map(function(u){return [s(u[0]),s(u[1]),s(u[2]),s(u[3]),u[4],u[5]];});
if(e.fa&&e.fa[6])e.fa[6]=e.fa[6].map(p);}return d;};"""


def encode_strings(shards):
    """Replace the repeated strings of every shard entry by their number in one table (FULL_X undoes it)."""
    table = {}
    S = lambda x: table.setdefault(x, len(table)) if isinstance(x, str) else x
    for entries in shards.values():
        for e in entries.values():
            e["go"] = [[S(a), S(b)] for a, b in e["go"]]
            e["kw"] = [S(k) for k in e["kw"]]
            if isinstance(e.get("ips"), list):
                e["ips"] = [[S(m[0]), S(m[1]), S(m[2]), m[3], m[4], m[5], S(m[6]), S(m[7])] for m in e["ips"]]
            if isinstance(e.get("uf"), list):
                e["uf"] = [[S(u[0]), S(u[1]), S(u[2]), S(u[3]), u[4], u[5]] for u in e["uf"]]
            if e.get("fa"):
                e["fa"][6] = [[S(a), S(b)] for a, b in e["fa"][6]]
    return sorted(table, key=table.get)


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
    ann, lengths, immune, jh, rep_of, members, uni, go, uf, km, ips, mm, fs, fst = load(folder)
    seqs = read_proteome(proteome, lengths) if proteome else {}
    # the Swiss-Prot proteins of the best hits with a position (MMseqs2, ProstT5 + Foldseek, Foldseek on
    # structure), so protein.html can align the two proteins; isoforms (P12345-2) are not in uniprot_sprot.fasta
    best_targets = {v[0] for d in (mm, fs, fst) if d for v in d.values() if v}
    tseqs = read_swissprot(swissprot, best_targets) if swissprot else {}
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
    rows, shards = [], defaultdict(dict)
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
        rows.append([pid, lengths[pid][0], r["existing_name"], r["proposed_name"], status_code[r["name_status"]],
                     r["sp_name"], r["sp_acc"], num(r["sp_pident"]), num(r["sp_qcov"]), r["sp_evalue"],
                     ", ".join(dict.fromkeys(n for _, n in pf)), mask,
                     "", size, drawn_page.get(pid, ""),  # 12 is empty: it held the Fester or BHF match
                     # 15: InterProScan's Pfam matches, [start, end, index into pfam]; 16: MMseqs2's best
                     # hit, [start, end, identity %, accession], for the index page's position track
                     [[x[3], x[4], pfam_ix.setdefault(x[1], len(pfam_ix))] for x in sorted(ips.get(pid, []), key=lambda x: x[3]) if x[0] == "Pfam"] if ips else [],
                     [mm[pid][2], mm[pid][3], round(100 * mm[pid][1], 1), mm[pid][0]] if mm and mm.get(pid) else None,
                     # 17: the other whole-proteome programs with a hit, as bits (OTHER_HITS)
                     sum(1 << i for i, d in enumerate((mm, fs, fst, km, jh)) if d and d.get(pid)),
                     1 if r["name_status"] == "conflict" and names_conflict(r["existing_name"], r["proposed_name"]) else 0])  # 18
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
            mm=None if mm is None else mm.get(pid, 0),
            fs=None if fs is None else fs.get(pid, 0),  # 0: no hit at FOLDSEEK_MAX_E or below
            fst=None if fst is None else ("not run" if lengths[pid][0] > STRUCT_MAX_AA else fst.get(pid, 0)),
            # the sequence (only for a protein without a drawn page; a drawn page carries its own), and the
            # lengths of its best hits' Swiss-Prot proteins, for the line under each hit
            q=seqs.get(pid) if not drawn_page.get(pid) else None,
            tl={v[0]: len(tseqs[v[0]]) for v in (mm and mm.get(pid), fs and fs.get(pid), fst and fst.get(pid)) if v and v[0] in tseqs},
            ips=None if ips is None else ("not run" if lengths[pid][0] > UNIFIRE_MAX_AA and not ips.get(pid) else sorted(ips.get(pid, []), key=lambda x: (x[3], x[4]))),  # 0: searched, no hit at E-value KMERSEEK_MAX_E or below
            uf=None if uf is None else ("not run" if lengths[pid][0] > UNIFIRE_MAX_AA and not uf.get(pid) else uf.get(pid, [])),
            card={k: card[k] for k in ("card_name", "one_line", "summary", "confidence", "confidence_reason", "name_verdict", "model")} if card else None,
            sp_e=r["sp_evalue"], ev=evidence_list(im["evidence"]) if im else [],
            jh=sorted(jh.get(pid, []), key=lambda x: float(x[3])), rep=rep,
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
        for page, pid in first.items():
            f = out / page
            html = f.read_text()
            html = html.replace("</main>", '<section id="biomni" data-id="%s"></section><section id="biomni-more"></section></main>' % pid, 1)
            html = html.replace("</body>", '<script src="../full/rows.js%s"></script><script src="../full/str.js%s"></script>'
                                '<script src="../assets/viewer.js%s"></script></body>' % (asset_ver, asset_ver, asset_ver), 1)
            if html.count('id="biomni') != 2:
                sys.exit("could not add the Biomni tables to %s" % page)
            f.write_text(html)
    shutil.copytree(Path(__file__).parent / "assets", out / "assets", dirs_exist_ok=True)
    n_tseq = add_targets(out, tseqs)
    strings = encode_strings(shards)
    (out / "full" / "str.js").write_text("window.FULL_STR=%s;%s" % (json.dumps(strings, separators=(",", ":")), STR_X_JS))
    for s in range(N_SHARDS):
        (out / "full" / ("%03d.js" % s)).write_text(
            "window.FULL_SHARD=Object.assign(window.FULL_SHARD||{},FULL_X(%s));" % json.dumps(shards.get(s, {}), separators=(",", ":")))
    (out / "full" / "rows.js").write_text("window.FULL=%s;" % json.dumps(
        dict(cats=cats, check=[t for _, t in NAME_CHECK], check_differs=status_code["conflict"], verdicts=VERDICT_TEXT, rows=rows, other=OTHER_HITS,
             pfam=[[a, pfam_desc.get(a, "")] for a in sorted(pfam_ix, key=pfam_ix.get)]), separators=(",", ":")))

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
        foldseek_hit=None if fs is None else len(fs),
        foldseek_struct_hit=None if fst is None else len(fst),
        drawn=sum(1 for r in rows if r[14]),
        drawn_pages=len({r[14] for r in rows if r[14]}),
        names_conflict=sum(1 for r in rows if r[18]),
        name_check={t: sum(1 for r in rows if r[4] == i) for i, (_, t) in enumerate(NAME_CHECK)},
        immune_categories={c: sum(1 for r in rows if r[11] >> i & 1) for i, c in enumerate(cats)},
    )
    (out / "full_summary.json").write_text(json.dumps(summ, indent=1))
    inputs = sorted(p.name for p in folder.iterdir() if p.suffix in (".tsv", ".gz") and p.name not in NOT_READ)
    (out / "index.html").write_text(index_html(summ, cats, inputs, asset_ver, known_html(drawn_page)))
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


def index_html(s, cats, inputs, ver, known=""):
    """The first page, after the dense-results-page skill: the answer in one sentence, then a sticky
    header (filters, search, a three-line key, column names, the residue ruler) over one grid row per
    protein, each with its Pfam matches and MMseqs2's best hit drawn along the protein."""
    f = "{:,}".format
    n = s["total"]
    pc = lambda v: "%d%%" % round(100 * v / n)
    chips = [("all", n, "all"), ("sp", s["swissprot_hit"], "Swiss-Prot match"), ("pfam", s["pfam_domain"], "Pfam domain"),
             ("none", s["neither"], "neither"), ("immune", s["immune_candidates"], "immune-gene candidate"),
             ("drawn", s["drawn"], "every program drawn"), ("names", s["names_conflict"], "names conflict")]
    tips = {"all": "every protein", "sp": "DIAMOND found a Swiss-Prot match at the cutoff", "pfam": "hmmsearch found a Pfam domain",
            "none": "neither a Swiss-Prot match at the cutoff nor a Pfam domain; another program may still have a hit",
            "immune": "in at least one immune-gene category", "drawn": "has a page that draws every program's hits along the protein",
            "names": "the FASTA header name and the proposed name are different names; capitals, punctuation and a trailing -like are ignored"}
    chip_html = "".join('<button class="chip%s" data-t="%s" title="%s">%s <b class="mono">%s</b></button>' % (" on" if k == "all" else "", k, tips[k], build.esc(w), f(v))
                        for k, v, w in chips)
    cat_opts = "".join('<option value="%d">%s (%s)</option>' % (i, build.esc(c), f(s["immune_categories"][c])) for i, c in enumerate(cats))
    chk_opts = "".join('<option value="%d">%s (%s)</option>' % (i, build.esc(t), f(s["name_check"][t])) for i, (_, t) in enumerate(NAME_CHECK))
    cat_bars = bars([(i, s["immune_categories"][c], c) for i, c in enumerate(cats)], n, "cat")
    chk_bars = bars([(i, s["name_check"][t], t) for i, (_, t) in enumerate(NAME_CHECK)], n, "chk")
    seq, dom = build.CATEGORIES["seq"]["color"], build.CATEGORIES["pfam"]["color"]
    ramp = "".join('<i style="opacity:%s;background:%s"></i>' % (o, seq) for o in (".3", ".5", ".75", "1"))
    # the same scale as full.js's track: residue 1 to AXIS_AA over the first TRACK_PCT% of the width
    at = lambda v: TRACK_PCT * (v - 1) / (AXIS_AA - 1)
    ticks = "".join('<text x="%.2f%%" y="10" class="ax" text-anchor="%s">%s</text><line x1="%.2f%%" x2="%.2f%%" y1="12" y2="16" stroke="#c3c2b7"/>'
                    % (at(v), "start" if v == 1 else "middle", f(v), at(v), at(v))
                    for v in (1, AXIS_AA // 4, AXIS_AA // 2, 3 * AXIS_AA // 4, AXIS_AA))
    return head("Botryllus proteome annotations", ver, dense=True) + """
<div class="titleblock">
 <div class="titlebar"><h1><i>Botryllus schlosseri</i> proteome</h1><span class="pn">%(n)s proteins predicted from the new genome: what each one matches</span>%(known)s</div>
 <div class="sentence"><p id="sent">%(sp)s of %(n)s proteins (%(spp)s) match a Swiss-Prot protein, %(pf)s (%(pfp)s) carry a Pfam domain, and %(ne)s (%(nep)s) have neither.</p>
 <div id="sentnote"><p class="note">Swiss-Prot match = the best hit of DIAMOND (a fast BLAST-like sequence search) at E-value 0.001 or below,
 in Swiss-Prot release 2026_03 (575,748 proteins). Pfam domain = a Pfam family model that pyhmmer's hmmsearch finds above the score Pfam sets
 for that family. Biomni, the agent that ran the whole-proteome programs except kmerseek, made both tables. A protein with neither can still have a hit
 from another program: the "other programs" column says which. %(dp)s proteins have a page that draws every program's hits along the protein
 (<a href="drawn.html">their own index</a>).</p></div></div>
 <details class="overview-fold"><summary>Counts by immune-gene category and by how the proposed name compares with the existing one</summary>
 <div class="overview">
 <section><h2>Immune-gene candidates, by category</h2><p class="note">%(im)s proteins; one protein can be in several categories. A category
 comes from a Pfam domain, the gene name, or an immune keyword or GO term of the best Swiss-Prot match; each protein page says which.
 A jackhmmer hit from an immune seed protein does not give a category: %(jn)s proteins have such a hit and no category. Bar length: share of
 the %(im)s; click a bar to list its proteins.</p>%(cb)s</section>
 <section><h2>Proposed name compared with the existing name</h2><p class="note">Existing = the name in the FASTA header. Proposed = from
 the best Swiss-Prot match or the Pfam domains. "Differs" includes a proposed name that only adds "-like" (NLRC5 and Nlrc5-like). Bars: of all
 %(n)s proteins.</p>%(kb)s</section>
 </div></details>
</div>
<header class="top ix">
 <div class="filters">
  <span class="lead2">DIAMOND E-value cutoff:</span><span class="ladder" id="eladder" title="Which DIAMOND hits count as a Swiss-Prot match. Click one."></span>
  <span class="lead2">show:</span><span class="ladder" id="tiles">%(chips)s</span>
  <label>immune category <select id="cat"><option value="-1">any</option>%(cat_opts)s</select></label>
  <label>names <select id="chk"><option value="-1">any</option>%(chk_opts)s</select></label>
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
 name in the FASTA header and the proposed name</dd>
 <dt>best Swiss-Prot match</dt><dd>DIAMOND's best hit (E-value 0.001 or below) in Swiss-Prot 2026_03, with its accession</dd>
 <dt>identity · coverage</dt><dd>identical residues in DIAMOND's alignment, and the share of your protein the alignment spans</dd>
 <dt>Pfam domains</dt><dd>Pfam models found by pyhmmer's hmmsearch at Pfam's gathering thresholds (Biomni). The boxes in the track are
 InterProScan's Pfam matches, which can differ</dd>
 <dt>cluster</dt><dd>proteins in the same MMseqs2 cluster (Biomni's clustering of the proteome), itself included</dd>
 <dt>other programs with a hit</dt><dd>%(other)s</dd>
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
<footer><h2>Input files read</h2><p class="note">Built by build_full.py from %(inputs)s. These tables come from Biomni's runs and our kmerseek run on Sherlock
(full_inputs/README.md in the repository says how each was made). This page only reads those tables; building it ran no search.</p></footer>
<script src="full/rows.js%(ver)s"></script><script src="assets/full.js%(ver)s"></script></body></html>
""" % dict(n=f(n), known=known, sp=f(s["swissprot_hit"]), spp=pc(s["swissprot_hit"]), pf=f(s["pfam_domain"]), pfp=pc(s["pfam_domain"]),
           ne=f(s["neither"]), nep=pc(s["neither"]), dp=f(s["drawn_pages"]), dr=f(s["drawn"]), im=f(s["immune_candidates"]),
           jn=f(s["jackhmmer_not_immune"]), cb=cat_bars, kb=chk_bars, chips=chip_html, cat_opts=cat_opts, chk_opts=chk_opts,
           other="; ".join("%s: %s" % (a, c) for a, c in OTHER_HITS), ax=f(AXIS_AA), axn=AXIS_AA, tw=TRACK_PCT, seq=seq, dom=dom, ramp=ramp, ticks=ticks, cnt="", inputs=build.esc(", ".join(inputs)), ver=ver)


def protein_html(ver, n_tseq=0):
    """protein.html: the dense page (assets/dense-bhf.js, as on a drawn page) for any protein, from
    the whole-proteome tables; assets/protein-dense.js builds its data from full/NNN.js. The key lists
    the programs with positions in those tables; the Swiss-Prot table lists the programs run on every protein."""
    run = {t["key"]: ("run" if t["key"] in ("mmseqs2", "prostt5", "foldseek", "interproscan") else "not run") for t in build.TRACKS}
    ctx = dict(run_status=run)
    whole = [["w_diamond"], ["w_jackhmmer"]]
    const = dict(
        tracks={t["key"]: dict(name=t["name"], color=build.PROGRAMS[t["program"]]["color"], shape=build.PROGRAMS[t["program"]]["shape"],
                               cat=build.PROGRAMS[t["program"]]["cat"], p=t["program"]) for t in build.TRACKS},
        program_names={k: v["name"] for k, v in build.PROGRAMS.items()}, ips_drawn=build.IPS_DRAWN, ips_evalue=sorted(build.IPS_EVALUE),
        fantasia_max=build.FANTASIA_MAX_DIST, kmerseek_max_e=KMERSEEK_MAX_E, struct_max_aa=STRUCT_MAX_AA,
        pfam_color=build.CATEGORIES["pfam"]["color"], ver=ver, tseq_files=n_tseq)
    head = build.page_head("Botryllus protein: what every annotation program found", 0).replace(
        "</head>", '<link rel="stylesheet" href="assets/dense-bhf.css%s"></head>' % ver, 1)
    help_ = ('<aside class="help" id="help"><h2>What this page is</h2><p>Every program\'s best hit on this protein, from Biomni\'s runs on '
             'the whole proteome (Biomni is the agent that ran them, except kmerseek, which ran on Sherlock). MMseqs2, ProstT5 + Foldseek and Foldseek on the '
             'predicted structure are drawn where they lie; InterProScan draws every match; kmerseek, DIAMOND, jackhmmer with immune seeds, '
             'FANTASIA and UniFIRE name a hit with no position. The 51 proteins with a drawn page also have phmmer, jackhmmer against '
             'Swiss-Prot, hmmscan, hhsearch and MMseqs2 three rounds, every hit, and the residues of each.</p><h2>Rows</h2><p>A <b>program '
             'row</b> folds its hits into one row; click it for one row per hit. A hit row gives the residues of your protein it covers, its '
             'E-value and score, the identity where the program gives one, and the matched protein; click it for its numbers.</p>'
             '<h2>Controls</h2><p>The <b>E-value cutoff</b> chips keep the hits at or below that value; the sentence and every count follow '
             'it. The best hits of MMseqs2 and the Foldseek runs are all at 0.001 or below, since the tables keep only those. <b>Search</b> '
             'looks through the matched protein\'s name, accession and species; every word typed must appear.</p><h2>Keyboard</h2><dl>'
             '<dt><kbd>j</kbd> <kbd>k</kbd></dt><dd>next and previous row</dd><dt><kbd>Enter</kbd> <kbd>Space</kbd></dt><dd>unfold a program, '
             'or open a hit</dd><dt><kbd>Esc</kbd></dt><dd>close every open hit and this help</dd><dt><kbd>?</kbd></dt><dd>this help</dd></dl></aside>')
    return head + """<div class="titleblock"><div class="titlebar"><h1 title="the ID as written in the proteome FASTA">loading…</h1><span class="pn"></span><span class="len"></span><a class="back" href="index.html">All proteins</a></div><div class="sentence"><p id="sent"></p><div id="sentnote"></div></div></div>
<header class="top">%s%s%s</header>
<main><div class="withviewer noviewer"><div id="grid"></div></div>
<section class="below" id="biomni" data-base="" data-undrawn="1"></section><section class="below" id="biomni-more"></section>%s</main>
<footer><h2>Input files read</h2><p class="note">Built by build_full.py from tables of the whole proteome, Biomni's and our kmerseek run's (full_inputs/README.md in the
repository says how each was made). This page only reads them; building it ran no search.</p></footer>
<div id="tip" role="tooltip" hidden></div>%s
<script>window.DENSE_CONST=%s;</script><script src="full/rows.js%s"></script><script src="full/str.js%s"></script><script src="assets/protein-dense.js%s"></script></body></html>
""" % (build.DENSE_FILTERS % (round(build.math.log10(build.SHOW_E)), round(build.math.log10(build.COUNT_E))),
       build.dense_key(ctx, False, False), build.DENSE_COLHEAD, build.dense_swissprot(ctx, whole), help_,
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
