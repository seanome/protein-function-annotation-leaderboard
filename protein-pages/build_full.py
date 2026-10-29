#!/usr/bin/env python3
"""Pages for every protein of the Botryllus proteome (45,339), from the Biomni annotation tables.

    python build_full.py full_inputs --drawn site --drawn-fasta test_inputs/data/proteins.fasta --out site_full

Reads, from the input folder (see full_inputs/README.md for where each file came from):
    botryllus_protein_annotations.tsv[.gz]   one row per protein: names, best Swiss-Prot hit, GO, Pfam
    botryllus_immune_genes.tsv                immune-gene candidates, their categories and evidence
    botryllus_jackhmmer_immune.tsv            jackhmmer hits of 20 immune seed proteins
    botryllus_allorecognition_hits.tsv        DIAMOND hits to Fester and BHF reference proteins
    allorecognition_references.tsv            NCBI titles of those reference proteins
    botryllus_mmseqs_clusters.tsv             MMseqs2 clusters (representative, member)
    proteome_lengths.tsv[.gz]                 each protein's length and the MD5 of its sequence

None of these tables gives where on the protein a hit lies, so these pages list hits in tables
and draw nothing on a ruler. Proteins also in a site made by build.py (--drawn) link to that
page, where every program's hits are drawn; they are matched by sequence (--drawn-fasta), since
the two FASTA files number isoforms differently.

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
NAME_CHECK = [
    ("", "no name, existing or proposed"),
    ("agree", "proposed name is the existing name"),
    ("conflict", "proposed name differs from the existing name"),
    ("unconfirmed", "existing name, nothing proposed"),
    ("new", "proposed name, no existing name"),
]


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


def load(folder):
    tsv = lambda name: list(csv.DictReader(opentext(inp(folder, name)), delimiter="\t"))
    ann = tsv("botryllus_protein_annotations.tsv")
    lengths = {r["protein_id"]: (int(r["length"]), r["md5"]) for r in tsv("proteome_lengths.tsv")}
    ids = [r["protein_id"] for r in ann]
    if len(set(ids)) != len(ids):
        sys.exit("botryllus_protein_annotations.tsv names a protein twice")
    if set(ids) != set(lengths):
        sys.exit("proteome_lengths.tsv and the annotation table name different proteins")
    immune = {r["protein_id"]: r for r in tsv("botryllus_immune_genes.tsv")}
    jh = defaultdict(list)
    for r in tsv("botryllus_jackhmmer_immune.tsv"):
        jh[r["target"]].append([r["seed_family"], r["seed"], int(r["iteration"]), r["e_value"], r["score"]])
    titles = {r["accession"]: re.sub(r"\s*\[Botryllus schlosseri\]$", "", r["title"])
              for r in tsv("allorecognition_references.tsv")}
    allo = defaultdict(list)
    # No header; DIAMOND's qseqid sseqid pident length evalue bitscore qlen slen (the qlen column
    # equals the FASTA length on every row, tests/test_build_full.py)
    for f in csv.reader(opentext(inp(folder, "botryllus_allorecognition_hits.tsv")), delimiter="\t"):
        allo[f[0]].append([f[1], titles.get(f[1], ""), float(f[2]), int(f[3]), f[4], float(f[5]), int(f[6]), int(f[7])])
    rep_of, members = {}, defaultdict(list)
    for f in csv.reader(opentext(inp(folder, "botryllus_mmseqs_clusters.tsv")), delimiter="\t"):
        rep_of[f[1]] = f[0]
        members[f[0]].append(f[1])
    missing = set(lengths) - set(rep_of)
    if missing:
        sys.exit("%d proteins have no row in botryllus_mmseqs_clusters.tsv, e.g. %s" % (len(missing), sorted(missing)[0]))
    for k in (set(immune) | set(jh) | set(allo) | set(rep_of)) - set(lengths):
        sys.exit("%s is in a table but not in the proteome" % k)
    return ann, lengths, immune, jh, allo, rep_of, members


def evidence_list(s):
    """'GO immune process:GO:0045087,GO:0051607 | NLR:Pfam:PF05729.19' -> [[category, detail], ...]"""
    out = []
    for part in s.split(" | "):
        cat, _, det = part.partition(":")
        if cat.strip():
            out.append([cat.strip(), det.strip()])
    return out


def build_site(folder, out, drawn=None, drawn_fasta=None):
    t0 = time.time()
    ann, lengths, immune, jh, allo, rep_of, members = load(folder)
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
    for r in ann:
        pid = r["protein_id"]
        im = immune.get(pid)
        mask = 0
        for c in split_list(im["categories"]) if im else []:
            mask |= 1 << cats.index(c)
        al = sorted(allo.get(pid, []), key=lambda a: -a[5])
        rep = rep_of[pid]
        size = len(members[rep])
        pf = pfam_list(r["pfam_domains"])
        rows.append([pid, lengths[pid][0], r["existing_name"], r["proposed_name"], status_code[r["name_status"]],
                     r["sp_name"], r["sp_acc"], num(r["sp_pident"]), num(r["sp_qcov"]), r["sp_evalue"],
                     ", ".join(dict.fromkeys(n for _, n in pf)), mask,
                     "%s, %.1f%%" % (al[0][1] or al[0][0], al[0][2]) if al else "", size, drawn_page.get(pid, "")])
        mem = [m for m in members[rep] if m != pid]
        shards[shard_of(pid)][pid] = dict(
            go=split_list(r["go_terms"]), kw=split_list(r["sp_keywords"]), pfam=pf,
            sp_e=r["sp_evalue"], ev=evidence_list(im["evidence"]) if im else [],
            jh=sorted(jh.get(pid, []), key=lambda x: float(x[3])), al=al, rep=rep,
            mem=[[m, name[m]] for m in mem[:MEMBERS_MAX]], more=max(0, len(mem) - MEMBERS_MAX))

    if drawn:
        shutil.copytree(drawn, out, dirs_exist_ok=True)
        (out / "index.html").rename(out / "drawn.html")
    (out / "full").mkdir(parents=True, exist_ok=True)
    asset_ver = "?v=" + hashlib.md5(b"".join(p.read_bytes() for p in sorted((Path(__file__).parent / "assets").glob("*.*")))).hexdigest()[:8]
    shutil.copytree(Path(__file__).parent / "assets", out / "assets", dirs_exist_ok=True)
    for s in range(N_SHARDS):
        (out / "full" / ("%03d.js" % s)).write_text(
            "window.FULL_SHARD=Object.assign(window.FULL_SHARD||{},%s);" % json.dumps(shards.get(s, {}), separators=(",", ":")))
    (out / "full" / "rows.js").write_text("window.FULL=%s;" % json.dumps(
        dict(cats=cats, check=[t for _, t in NAME_CHECK], rows=rows), separators=(",", ":")))

    n = len(rows)
    summ = dict(
        total=n,
        swissprot_hit=sum(1 for r in rows if r[6]),
        pfam_domain=sum(1 for r in rows if r[10]),
        neither=sum(1 for r in rows if not r[6] and not r[10]),
        immune_candidates=sum(1 for r in rows if r[11]),
        jackhmmer_hit=sum(1 for pid in jh if pid in lengths),
        jackhmmer_not_immune=sum(1 for pid in jh if pid not in immune),
        allorecognition_hit=sum(1 for r in rows if r[12]),
        in_cluster_of_2_or_more=sum(1 for r in rows if r[13] > 1),
        clusters=len(members),
        drawn=sum(1 for r in rows if r[14]),
        drawn_pages=len({r[14] for r in rows if r[14]}),
        name_check={t: sum(1 for r in rows if r[4] == i) for i, (_, t) in enumerate(NAME_CHECK)},
        immune_categories={c: sum(1 for r in rows if r[11] >> i & 1) for i, c in enumerate(cats)},
    )
    (out / "full_summary.json").write_text(json.dumps(summ, indent=1))
    inputs = sorted(p.name for p in folder.iterdir() if p.suffix in (".tsv", ".gz"))
    (out / "index.html").write_text(index_html(summ, cats, inputs, asset_ver))
    (out / "protein.html").write_text(protein_html(asset_ver))
    print("%d proteins, %d drawn, %d data files, %.1f s" % (n, summ["drawn"], N_SHARDS, time.time() - t0))
    return summ


def head(title, ver):
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1"><title>%s</title>%s'
            '<link rel="stylesheet" href="assets/site.css%s"></head><body>' % (build.esc(title), build.CSS_LINKS, ver))


def bars(items, total, key):
    """Counts as horizontal bars of one ink colour; each line filters the table when clicked."""
    top = max(v for _, v, _ in items) or 1
    rows = "".join(
        '<button class="bar" data-%s="%s"><span class="blab">%s</span><span class="bt"><i style="width:%.1f%%"></i></span>'
        '<span class="bn">%s</span></button>' % (key, k, build.esc(t), 100 * v / top, "{:,}".format(v))
        for k, v, t in items)
    return '<div class="bars">%s</div>' % rows


def index_html(s, cats, inputs, ver):
    f = "{:,}".format
    nums = [("all", s["total"], "proteins"), ("sp", s["swissprot_hit"], "with a Swiss-Prot hit"),
            ("pfam", s["pfam_domain"], "with a Pfam domain"), ("none", s["neither"], "with neither"),
            ("immune", s["immune_candidates"], "immune-gene candidates"),
            ("allo", s["allorecognition_hit"], "with a Fester or BHF hit"),
            ("drawn", s["drawn"], "with every program's hits drawn")]
    tiles = "".join('<button class="num" data-t="%s"><span class="big">%s</span><span>%s</span></button>' % (k, f(v), build.esc(t))
                    for k, v, t in nums)
    cat_bars = bars([(i, s["immune_categories"][c], c) for i, c in enumerate(cats)], s["total"], "cat")
    chk_bars = bars([(i, s["name_check"][t], t) for i, (_, t) in enumerate(NAME_CHECK)], s["total"], "chk")
    return head("Botryllus proteome annotations", ver) + """
<header class="wide"><h1>Every <i>Botryllus schlosseri</i> protein, and what is known about it</h1>
<p class="summary">%s proteins of the new genome's predicted proteome. For each: its best Swiss-Prot match, its Pfam domains,
whether it looks like an immune gene, and its matches to the allorecognition proteins Fester and BHF.</p></header>
<main class="wide">
<section class="nums">%s</section>
<p class="note">Click a number to list those proteins. %s pages draw every annotation program's hits along the protein
(<a href="drawn.html">their own index</a>); %s proteins link to one, since some proteins have the same sequence.</p>
<div class="twocol">
<section><h2>Immune-gene candidates, by category</h2>
<p class="note">%s proteins; one protein can be in several categories. A category comes from a Pfam domain, the gene
name, an immune keyword or GO term of the best Swiss-Prot match, or a Fester or BHF hit; each protein page says which.
jackhmmer hits from immune seed proteins do not set a category: %s proteins have one and are in no category.
Click a bar to list its proteins.</p>%s</section>
<section><h2>Proposed name compared with the existing name</h2>
<p class="note">The existing name is the one in the FASTA header. The proposed name comes from the best Swiss-Prot match
or the Pfam domains. "Differs" includes a proposed name that only adds "-like" (NLRC5 and Nlrc5-like).</p>%s</section>
</div>
<h2>Proteins</h2>
<div class="search"><label for="q">Search</label><input id="q" type="search" autocomplete="off"
 placeholder="an ID, gene, protein name, Pfam domain or category, e.g. NLRC5, lectin, TIR, Fester"></div>
<p class="note">The search reads the ID, the existing and proposed names, the best Swiss-Prot match's name and accession,
the Pfam domain names, the immune categories and the Fester or BHF match. Every word typed must appear.</p>
<p id="count" class="note"></p>
<div class="tablewrap"><table id="t" class="full"><thead><tr>
<th data-k="0" data-t="s">Protein</th><th data-k="1">Length (aa)</th><th data-k="2" data-t="s">Existing name</th>
<th data-k="3" data-t="s">Proposed name</th><th data-k="5" data-t="s">Best Swiss-Prot match</th>
<th data-k="7">Identity (%%)</th><th data-k="8">Coverage of your protein (%%)</th>
<th data-k="10" data-t="s">Pfam domains</th><th>Immune categories</th><th data-k="12" data-t="s">Fester or BHF match (identity)</th>
<th data-k="13">Proteins in its MMseqs2 cluster</th><th>Drawn</th></tr></thead><tbody></tbody></table></div>
<p><button id="more" hidden>Show 100 more</button></p>
</main>
<footer><h2>Input files read</h2><p class="note">Built by build_full.py from %s. These tables come from a Biomni run
(full_inputs/README.md in the repository says how each was made). Identity and coverage are DIAMOND's, against
Swiss-Prot. No tool was run to make these pages.</p></footer>
<script src="full/rows.js%s"></script><script src="assets/full.js%s"></script></body></html>
""" % (f(s["total"]), tiles, f(s["drawn_pages"]), f(s["drawn"]), f(s["immune_candidates"]), f(s["jackhmmer_not_immune"]), cat_bars, chk_bars,
       build.esc(", ".join(inputs)), ver, ver)


def protein_html(ver):
    return head("Botryllus protein", ver) + """
<header class="wide"><a class="back" href="index.html">All proteins</a><h1 class="mono" id="pid"></h1>
<p class="len" id="sub"></p></header>
<main class="wide" id="main"><p class="note">Loading…</p></main>
<script src="full/rows.js%s"></script><script src="assets/viewer.js%s"></script></body></html>
""" % (ver, ver)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--drawn", type=Path, help="a site made by build.py, to link to and copy in")
    ap.add_argument("--drawn-fasta", type=Path, help="the proteins.fasta that site was built from")
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
    build_site(a.inputs, a.out, a.drawn, a.drawn_fasta)


if __name__ == "__main__":
    main()
