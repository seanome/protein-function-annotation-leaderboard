#!/usr/bin/env python3
"""Build a static website that shows, one page per protein, what each annotation tool found.

    python build.py INPUT_DIR [--out site] [--jobs 8]

INPUT_DIR holds the files listed in README.md. The script only reads them; it runs no search
tool and makes no network call. It writes INPUT_DIR's results as plain HTML into --out.
Needs Python 3.9 or newer and Biopython (for the recomputed alignments).
"""
import argparse
import csv
import datetime
import gc
import gzip
import hashlib
import html
import json
import math
import multiprocessing
import os
import re
import shutil
import sys
import time
from pathlib import Path

E_MAX = 10.0  # a hit is drawn when its E-value is at most this, the -E 10 the searches were run with
COUNT_E = 1e-3  # a hit counts (sentence, index, summary.json) when its E-value is at most this; --count-at
TARGET_PFAM_E = 1e-3  # a Swiss-Prot protein's Pfam families: its hmmscan or InterProScan domains at or below this
HH_PROB_MIN = 50.0  # hhsearch hits also need Prob of at least this
OVERLAP = 0.5  # two intervals overlap when they share at least half of the shorter one
TOP_N = 5  # hits drawn per track before "show all"

# Colour is the kind of program (three hues and a dark grey); inside a kind, the box's shape
# tells the programs apart. Checked in greyscale and with a deuteranopia simulation (Machado
# et al. 2009): the four colours differ by at least 7 in lightness (CIE L*) and by at least 23
# in colour difference (CIE76) under deuteranopia. None is a pLDDT colour, and none is the
# orange or blue of the residue panel's hydrophobic and polar classes.
CATEGORIES = {
    "seq": {"name": "sequence search", "color": "#008B7A"},
    "struct": {"name": "structure", "color": "#6A3D9A"},
    "kmer": {"name": "k-mer", "color": "#C0392B"},
    "pfam": {"name": "Pfam-model scan", "color": "#3D3D3D"},
}
SHAPES = ("filled", "open", "hatched", "rounded")
PROGRAMS = {
    "mmseqs2": {"name": "MMseqs2", "cat": "seq", "shape": "filled"},
    "kmerseek": {"name": "kmerseek", "cat": "kmer", "shape": "filled"},
    "foldseek": {"name": "Foldseek", "cat": "struct", "shape": "filled"},
    "phmmer": {"name": "phmmer", "cat": "seq", "shape": "open"},
    "jackhmmer": {"name": "jackhmmer", "cat": "seq", "shape": "hatched"},
    "hhsearch": {"name": "HH-suite (hhsearch)", "cat": "seq", "shape": "rounded"},
    "hmmscan": {"name": "hmmscan with Pfam", "cat": "pfam", "shape": "filled"},
    "interproscan": {"name": "InterProScan", "cat": "pfam", "shape": "open"},
}
for _p in PROGRAMS.values():
    _p["color"] = CATEGORIES[_p["cat"]]["color"]
PROTEIN_LINE = "#C8C8C8"  # the protein's backbone line, and only that
AGREEMENT = "#000000"
EXPRESSION = "#1a1a1a"  # expression bars are thin outlines in text ink; their shape, not a grey, marks them
PLDDT = [  # AlphaFold's own colours and bins
    (90, 101, "#0053D6", "very high (above 90)"),
    (70, 90, "#65CBF3", "confident (70 to 90)"),
    (50, 70, "#FFDB13", "low (50 to 70)"),
    (0, 50, "#FF7D45", "very low (below 50)"),
]

# The ten tracks, in drawing order. "second" marks the second run of a program (MMseqs2 three
# rounds, ProstT5 + Foldseek); it takes its program's colour and shape and has its own row.
# Settings are the leaderboard's setting lines; run_settings.tsv can replace them.
TRACKS = [
    dict(key="mmseqs2", file="mmseqs2.m8", program="mmseqs2", second=False, name="MMseqs2",
         setting="sensitivity 7 (-s 7), one round", database="Swiss-Prot",
         score="evalue", kind="m8"),
    dict(key="mmseqs2_3rounds", file="mmseqs2_3rounds.m8", program="mmseqs2", second=True,
         name="MMseqs2, three rounds", setting="sensitivity 7 (-s 7), three rounds",
         database="Swiss-Prot", score="evalue", kind="m8"),
    dict(key="kmerseek", file="kmerseek.tsv", program="kmerseek", second=False, name="kmerseek",
         setting="", database="Swiss-Prot", score="region_evalue", kind="kmerseek"),
    dict(key="foldseek", file="foldseek.m8", program="foldseek", second=False, name="Foldseek",
         setting="3Di shape letters and amino acids together",
         database="AlphaFold models of Swiss-Prot", score="evalue", kind="m8"),
    dict(key="prostt5", file="prostt5_foldseek.m8", program="foldseek", second=True,
         name="ProstT5 + Foldseek", setting="3Di letters predicted from sequence by ProstT5",
         database="AlphaFold models of Swiss-Prot", score="evalue", kind="m8"),
    dict(key="phmmer", file="phmmer.domtblout", program="phmmer", second=False, name="phmmer",
         setting="phmmer's default options, E-value cutoff 10", database="Swiss-Prot",
         score="i-Evalue", kind="domtbl_seq"),
    dict(key="jackhmmer", file="jackhmmer.domtblout", program="jackhmmer", second=False,
         name="jackhmmer", setting="three rounds of search (-N 3); the last round is drawn",
         database="Swiss-Prot", score="i-Evalue", kind="domtbl_seq"),
    dict(key="hhsearch", file="hhsearch", program="hhsearch", second=False,
         name="HH-suite (hhsearch)", setting="one .hhr file per protein",
         database="Swiss-Prot or PDB70 profiles", score="Prob and E-value", kind="hhr"),
    dict(key="hmmscan", file="hmmscan.domtblout", program="hmmscan", second=False,
         name="hmmscan with Pfam", setting="scans Pfam's models; no other protein searched",
         database="Pfam-A", score="i-Evalue", kind="domtbl_pfam"),
    dict(key="interproscan", file="interproscan.tsv", program="interproscan", second=False,
         name="InterProScan", setting="its member databases' own cutoffs", database="InterPro",
         score="score column", kind="interproscan"),
]
TRACK = {t["key"]: t for t in TRACKS}
PROGRAM_ORDER = list(PROGRAMS)

# InterProScan analyses that predict a feature of the protein (a coil, a membrane span, a
# signal peptide, disorder) rather than a family or domain; they do not count as a hit.
IPS_FEATURES = {"Coils", "MobiDBLite", "Phobius", "TMHMM", "DeepTMHMM", "SignalP_EUK",
                "SignalP_GRAM_POSITIVE", "SignalP_GRAM_NEGATIVE", "SignalP-noTM", "SignalP-TM",
                "SignalP", "AntiFam"}
# InterProScan analyses whose score column is an E-value; for the others the score is not an
# E-value and the match counts if the member database reported it.
IPS_EVALUE = {"Pfam", "SMART", "SUPERFAMILY", "Gene3D", "CDD", "NCBIfam", "TIGRFAM", "PIRSF",
              "PANTHER", "HAMAP", "SFLD", "FunFam", "PIRSR", "AntiFam"}


# ----------------------------------------------------------------------------- reading ----

def opentext(path):
    """A text file, or its gzipped copy when the name ends in .gz."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def accession(t):
    """Swiss-Prot accession from a target name: sp|P12345|NAME, AF-P12345-F1-model_v4, P12345."""
    if t.startswith("AF-"):
        return t.split("-")[1]
    if t.startswith(("sp|", "tr|")):
        parts = t.split("|")
        return parts[1] if len(parts) > 1 else t
    return t


def counts(h):
    """A drawn hit counts when its E-value is at most COUNT_E. A hit with no E-value (an
    InterProScan match from a member database whose score is not an E-value) counts if drawn."""
    return h["passes"] and (h["evalue"] is None or h["evalue"] <= COUNT_E)


def fmt_thr(e):
    """0.001, 10, 1e-05: how a threshold is written in page text."""
    return "%g" % e


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def hit(track, q0, q1, target, t0, t1, evalue, scores, passes, aln=None):
    return dict(track=track, q0=int(q0), q1=int(q1), target=target, t0=t0, t1=t1,
                evalue=evalue, scores=scores, passes=passes, aln=aln)


M8_COLS = "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits".split(",")


# The two big kinds of file (MMseqs2/Foldseek tables, phmmer/jackhmmer tables) can hold tens of
# millions of hits for a proteome. While reading, only the raw text line of each hit at the
# threshold is kept (as (track, line)); the worker drawing a protein's page turns its lines
# into hit records with m8_hit / domtbl_hit. Holding every hit as a record took over 32 GB for
# 7,216 Botryllus proteins.
USED_TARGETS = set()  # Swiss-Prot accessions of the search hits at the threshold
PFAM_DESC = {}  # Pfam accession -> its one-line description, from hmmscan's last column
# The tools' own alignments, from HMMER's main output (phmmer.out, jackhmmer.out, hmmscan.out):
# (track, query, target, hmm from, hmm to, ali from, ali to) -> (model row, target row), the
# two aligned rows exactly as HMMER printed them.
HMMER_ALN = {}


def m8_hit(line, key):
    """One MMseqs2 or Foldseek row. Two optional last columns, qaln and taln (the aligned
    sequences, from --format-output ...,qaln,taln), give the tool's own alignment."""
    cols = M8_COLS + (["alntmscore", "theader"] if key in ("foldseek", "prostt5") else ["theader"]) + ["qaln", "taln"]
    d = dict(zip(cols, line.rstrip("\n").split("\t")))
    e = fnum(d["evalue"])
    sc = [("evalue", d["evalue"]), ("fident", d.get("fident", "")), ("bits", d.get("bits", ""))]
    if key in ("foldseek", "prostt5"):
        sc.append(("alntmscore", d.get("alntmscore", "")))
    sc = [(k, v if v not in ("", None) else "not in the file") for k, v in sc]
    aln = None
    if d.get("qaln") and d.get("taln"):
        aln = dict(q=d["qaln"], t=d["taln"], qs=int(d["qstart"]), ts=int(d["tstart"]),
                   src="the alignment %s wrote (the qaln and taln columns of %s)" % (TRACK[key]["name"], TRACK[key]["file"]))
    return hit(key, d["qstart"], d["qend"], accession(d["target"]), int(d["tstart"]),
               int(d["tend"]), e, sc, e is not None and e <= E_MAX, aln)


def parse_m8(path, key):
    out = {}
    with opentext(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("query\t"):
                continue
            f = line.split("\t", 9)
            e = fnum(f[8])
            if e is None or e > E_MAX:  # only hits at the threshold are drawn
                continue
            out.setdefault(f[0], []).append((key, line))
            USED_TARGETS.add(accession(f[1]))
    return out


def domtbl_rows(path):
    with opentext(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            yield line, line.split(None, 22)


def parse_domtbl_seq(path, key):
    """phmmer / jackhmmer: the query (your protein) became the model, so its positions are the
    'hmm coord' columns; 'ali coord' and 'env coord' are positions on the Swiss-Prot protein.
    For jackhmmer only the last round is kept (see last_round).
    Rows come grouped by query (HMMER finishes one query before the next), so the file is read
    one query at a time. Only hits at the threshold are kept: they are the only ones drawn."""
    out = {}

    def flush(q, rows):
        if key == "jackhmmer":
            rows = last_round(rows)
        for line, f in rows:
            e = fnum(f[12])
            if e is None or e > E_MAX:
                continue
            out.setdefault(q, []).append((key, line))
            USED_TARGETS.add(accession(f[0]))

    cur_q, rows = None, []
    for line, f in domtbl_rows(path):
        if f[3] != cur_q:
            if rows:
                flush(cur_q, rows)
            cur_q, rows = f[3], []
        rows.append((line, f))
    if rows:
        flush(cur_q, rows)
    return out


def domtbl_hit(line, key):
    f = line.split(None, 22)
    e = fnum(f[12])
    sc = [("i-Evalue", f[12]), ("score", f[13]), ("target env from-to", f[19] + "-" + f[20])]
    h = hit(key, f[15], f[16], accession(f[0]), int(f[17]), int(f[18]), e, sc, e is not None and e <= E_MAX)
    raw = HMMER_ALN.get((key, f[3], f[0], int(f[15]), int(f[16]), int(f[17]), int(f[18])))
    if raw:
        h["raw_aln"] = raw  # turned into rows in build_protein, which has the query's residues
    return h


def parse_hmmer_out(path, key, wanted):
    """The domain alignments in HMMER's main output (-o, best with --notextw), for queries in
    `wanted`. For jackhmmer, a query's later round replaces its earlier ones, so the last
    round is kept, as for the table. Each block after '== domain N' holds, in order: optional
    RF/CS lines, the model row ('name from residues to'), the match line, the target row
    ('name from residues to') and a PP line."""
    row = re.compile(r"^\s*(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s*$")
    q, tgt, per_q, block = None, None, {}, None

    def finish_block():
        if block and len(block) >= 2:
            (mn, mf, ms, mt), (tn, tf, ts, tt) = block[0], block[1]
            per_q[(key, q, tgt, int(mf), int(mt), int(tf), int(tt))] = (ms, ts, int(mf), int(tf))

    with opentext(path) as fh:
        for line in fh:
            if line.startswith("Query:"):
                finish_block(); block = None
                if q in wanted:
                    HMMER_ALN.update(per_q)
                q, per_q = line.split()[1], {}
            elif line.startswith("@@ Round:"):
                finish_block(); block = None
                per_q = {}
            elif line.startswith(">> "):
                finish_block(); block = None
                tgt = line[3:].split()[0]
            elif line.lstrip().startswith("== domain"):
                finish_block(); block = []
            elif block is not None:
                st = line.strip()
                if not st:
                    if block:
                        finish_block(); block = None
                    continue
                if st.endswith(" RF") or st.endswith(" CS"):
                    continue
                m = row.match(line)
                if m and (len(block) == 0 or m.group(1) == tgt or len(block) == 1):
                    if len(block) < 2:
                        block.append(m.groups())
    finish_block()
    if q in wanted:
        HMMER_ALN.update(per_q)


def materialize(x):
    """A hit record from what the reader kept: a record already, or (track, raw line)."""
    if isinstance(x, dict):
        return x
    key, line = x
    return m8_hit(line, key) if TRACK[key]["kind"] == "m8" else domtbl_hit(line, key)


def last_round(rows):
    """rows: (line, fields) pairs of one query."""
    """Keep the last jackhmmer round if the table holds more than one. HMMER 3.4's jackhmmer
    writes --domtblout for the last round only (on the 50-protein test set every protein has
    one block), so this is a guard for other versions. A new round is taken to start when a
    target comes back that was already seen, or when a new target has a better full-sequence
    E-value than the one before it (inside a round, targets are best first)."""
    blocks, cur, seen, prev_t, prev_e = [], [], set(), None, None
    for item in rows:
        f = item[1] if isinstance(item, tuple) else item
        t, e = f[0], fnum(f[6])
        if t != prev_t and cur and (t in seen or (e is not None and prev_e is not None and e < prev_e)):
            blocks.append(cur)
            cur, seen = [], set()
        cur.append(item)
        seen.add(t)
        prev_t, prev_e = t, e
    if cur:
        blocks.append(cur)
    return blocks[-1] if blocks else []


def parse_domtbl_pfam(path, proteins):
    """hmmscan: query = protein sequence, target = Pfam model; env coords are on the query.
    Rows whose query is not one of your proteins are Swiss-Prot proteins that were also
    scanned; they give the Pfam domains of a search hit's target."""
    out, target_pfam, pfam_names = {}, {}, {}
    for _, f in domtbl_rows(path):
        pf = f[1].split(".")[0] if f[1] != "-" else f[0]
        pfam_names[pf] = f[0]
        if len(f) > 22 and f[22].strip():
            PFAM_DESC[pf] = f[22].strip()
        e = fnum(f[12])
        if f[3] in proteins:
            sc = [("i-Evalue", f[12]), ("score", f[13]), ("Pfam", pf)]
            h = hit("hmmscan", f[19], f[20], f[0], None, None, e, sc, e is not None and e <= E_MAX)
            h["names"] = {pf}
            h["hmm"] = [int(f[15]), int(f[16]), int(f[2])]  # model from, to, model length
            raw = HMMER_ALN.get(("hmmscan", f[3], f[0], int(f[15]), int(f[16]), int(f[17]), int(f[18])))
            if raw:
                h["raw_aln"] = raw
            out.setdefault(f[3], []).append(h)
        elif e is not None and e <= TARGET_PFAM_E:
            target_pfam.setdefault(accession(f[3]), []).append((int(f[19]), int(f[20]), {pf}))
    return out, target_pfam, pfam_names


def parse_interproscan(path, proteins):
    out, target_ann, names = {}, {}, {}
    with opentext(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            prot, analysis, sig, desc = f[0], f[3], f[4], f[5]
            ipr = f[11] if len(f) > 11 and f[11] not in ("", "-") else None
            ipr_desc = f[12] if len(f) > 12 else ""
            score = f[8]
            e = fnum(score)
            if analysis in IPS_FEATURES:
                passes = False
            elif analysis in IPS_EVALUE and e is not None:
                passes = e <= E_MAX
            else:
                passes = True
            nm = set()
            if analysis == "Pfam":
                nm.add(sig.split(".")[0])
                names[sig.split(".")[0]] = desc or sig
            if ipr:
                nm.add(ipr)
                names[ipr] = ipr_desc or ipr
            if prot in proteins:
                sc = [("member database", analysis), ("signature", sig), ("score", score),
                      ("InterPro entry", (ipr + " " + ipr_desc) if ipr else "none")]
                # the score is kept as an E-value only where it is one
                h = hit("interproscan", f[6], f[7], desc or sig, None, None,
                        e if analysis in IPS_EVALUE else None, sc, passes)
                h["names"] = nm
                h["db"] = analysis
                out.setdefault(prot, []).append(h)
            elif passes and nm and (e is None or e <= TARGET_PFAM_E):
                target_ann.setdefault(accession(prot), []).append((int(f[6]), int(f[7]), nm))
    return out, target_ann, names


HHR_HIT = re.compile(r"^\s*(\d+)\s+(.*?)\s+(\d+(?:\.\d+)?)\s+(\S+)\s+(\S+)\s+(-?\d+(?:\.\d+)?)\s+"
                     r"(-?\d+(?:\.\d+)?)\s+(\d+)\s+(\d+)-(\d+)\s+(\d+)-(\d+)\s*\(\s*(\d+)\)\s*$")


def parse_hhr(path):
    """One .hhr file: the hit table (query and template ranges, Prob, E-value) and, for each
    hit, the alignment of the query and template residues as hhsearch printed it."""
    query, hits, alns, cur = None, {}, {}, None
    with opentext(path) as fh:
        lines = fh.read().splitlines()
    in_table = False
    for line in lines:
        if line.startswith("Query ") and query is None:
            query = line.split()[1]
        if line.startswith(" No Hit"):
            in_table = True
            continue
        if in_table:
            if not line.strip():
                in_table = False
                continue
            m = HHR_HIT.match(line)
            if m:
                no = int(m.group(1))
                hits[no] = dict(name=m.group(2).split()[0], prob=float(m.group(3)),
                                evalue=fnum(m.group(4)), q0=int(m.group(9)), q1=int(m.group(10)),
                                t0=int(m.group(11)), t1=int(m.group(12)))
            continue
        if line.startswith("No "):
            cur = int(line.split()[1])
            alns[cur] = {"q": [], "t": [], "qs": None, "ts": None}
        elif cur is not None and (line.startswith("Q ") or line.startswith("T ")):
            parts = line.split()
            if len(parts) >= 5 and parts[1] not in ("Consensus",) and not parts[1].startswith("ss_"):
                side = "q" if line[0] == "Q" else "t"
                # Q name start SEQ end (len)
                try:
                    start = int(parts[2])
                except ValueError:
                    continue
                alns[cur][side].append(parts[3])
                if alns[cur][side + "s"] is None:
                    alns[cur][side + "s"] = start
    out = []
    for no, h in sorted(hits.items()):
        e = h["evalue"]
        a = alns.get(no)
        aln = None
        if a and a["q"] and a["t"]:
            aln = dict(q="".join(a["q"]), t="".join(a["t"]), qs=a["qs"], ts=a["ts"])
        sc = [("Prob", "%g" % h["prob"]), ("E-value", "%g" % e if e is not None else "")]
        passes = e is not None and e <= E_MAX and h["prob"] >= HH_PROB_MIN
        out.append(hit("hhsearch", h["q0"], h["q1"], accession(h["name"]), h["t0"], h["t1"], e,
                       sc, passes, aln=aln))
    return query, out


def parse_kmerseek(path):
    """kmerseek regions. Column names are used as they are in the file. Region coordinates in
    kmerseek's output are 0-based and end-exclusive, so 1 is added to the start."""
    out, meta = {}, {}
    with opentext(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        col = {c: i for i, c in enumerate(head)}

        def pick(*names):
            for n in names:
                if n in col:
                    return n
            raise SystemExit("kmerseek.tsv has none of the columns %s" % (names,))
        cq, ct = pick("query", "query_name"), pick("target", "target_name")
        cqs, cqe = pick("query_start", "region_start"), pick("query_end", "region_end")
        cts, cte = pick("target_start"), pick("target_end")
        ce = pick("region_evalue")
        go_cols = [c for c in head if re.search(r"(^|_)go($|_)|pfam2go", c, re.I)]
        meta["go_cols"] = go_cols
        meta["moltype"] = set()
        meta["ksize"] = set()
        for f in rd:
            d = dict(zip(head, f))
            e = fnum(d[ce])
            sc = [(ce, d[ce])]
            for g in go_cols:
                sc.append((g, d[g] or "none"))
            if not go_cols:
                sc.append(("GO terms", "no GO column in the file"))
            if "moltype" in d:
                meta["moltype"].add(d["moltype"])
            if "ksize" in d:
                meta["ksize"].add(d["ksize"])
            h = hit("kmerseek", int(d[cqs]) + 1, int(d[cqe]), accession(d[ct]), int(d[cts]) + 1,
                    int(d[cte]), e, sc, e is not None and e <= E_MAX)
            h["moltype"] = d.get("moltype")
            if h["passes"]:  # only hits at the threshold are drawn
                out.setdefault(d[cq], []).append(h)
                USED_TARGETS.add(h["target"])
    return out, meta


def parse_alphabets(path):
    """Residue classes of every kmerseek alphabet, read from kmerseek's own source file
    src/rust/alphabets.rs: {moltype: {residue: class letter}}."""
    src = Path(path).read_text()
    variant_to_moltype = dict(re.findall(r'Self::(\w+)\s*=>\s*"([a-z0-9_]+)"', src))
    variant_to_static = dict(re.findall(r"Self::(\w+)\s*=>\s*&([A-Z0-9_]+)", src))
    tables = {}
    for m in re.finditer(r"static ([A-Z0-9_]+):\s*LazyLock<HashMap<u8, u8>>\s*=\s*LazyLock::new\(\|\|\s*(\w+)\(([^)]*)\)\)", src, re.S):
        name, fn, args = m.groups()
        if fn in ("build_hp", "build_hpc"):
            groups = re.findall(r'b"([A-Z]+)"', args)
            letters = ["h", "p", "c"][: len(groups)]
            tables[name] = {r: letters[i] for i, g in enumerate(groups) for r in g}
        elif fn == "build_reduced":
            const = args.strip()
            cm = re.search(r"(?:const|static) %s:[^=]*=\s*&\[(.*?)\];" % re.escape(const), src, re.S)
            if cm:
                clusters = re.findall(r'"([A-Z]+)"', cm.group(1))
                tables[name] = {r: c[0].lower() for c in clusters for r in c}
    out = {}
    for v, mt in variant_to_moltype.items():
        st = variant_to_static.get(v)
        if st in tables:
            out[mt] = tables[st]
    out["protein20"] = {r: r.lower() for r in "ACDEFGHIKLMNPQRSTVWY"}
    return out


def read_tsv_table(path):
    with opentext(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        return head, [r for r in rd if r]


# ------------------------------------------------------------------------ agreement ----

def overlaps(a0, a1, b0, b1):
    inter = min(a1, b1) - max(a0, b0) + 1
    if inter <= 0:
        return False
    return inter >= OVERLAP * min(a1 - a0 + 1, b1 - b0 + 1)


def target_names(h, target_ann):
    """Pfam families and InterPro entries on the matched part of a search hit's target, or None
    when nothing is known about that target."""
    if "names" in h:
        return h["names"]
    doms = target_ann.get(h["target"])
    if doms is None or h["t0"] is None:
        return None
    lo, hi = min(h["t0"], h["t1"]), max(h["t0"], h["t1"])
    return set().union(*[n for (d0, d1, n) in doms if overlaps(lo, hi, d0, d1)] or [set()])


AGREE_OTHERS = [[]]  # the other agreeing Pfam families of the last protein agreement() looked at
# For a 'disagree' protein: [(Pfam family, [programs that named it]), ...] from the overlapping
# pairs of hits, best E-value first.
DISAGREE_NAMES = [[]]


def agreement(hits, target_ann):
    """Returns (programs with a counted hit, status, agreed name id).
    status, for two or more programs:
      'agree': hits of two programs overlap and name the same Pfam family or InterPro entry;
      'disagree': hits of two programs overlap, both name a family, and none is shared;
      'overlap_unknown': hits of two programs overlap, but no family is known for one of them;
      'apart': no hit of one program overlaps a hit of another.
    'disagree' is checked before 'overlap_unknown'; for one program, 'one'; for none, 'none'.

    Proteins with many hits (a kinase can have thousands) make comparing every pair slow, so
    hits are grouped by the family name they carry, repeated ranges are dropped, and the search
    for an overlapping pair between two programs stops at the first one found, best
    E-value first. This gives the same answer as comparing every pair."""
    passing = [h for h in hits if counts(h)]
    progs = sorted({TRACK[h["track"]]["program"] for h in passing}, key=PROGRAM_ORDER.index)
    AGREE_OTHERS[0] = []
    DISAGREE_NAMES[0] = []
    if not progs:
        return progs, "none", None
    if len(progs) == 1:
        return progs, "one", None
    # by_name: family -> program -> {range: best E-value}; known: program -> {range: (best
    # E-value, families)} for hits that name a family; unknown_hits: the same for hits that name
    # none (target not annotated, or no family on the matched part).
    by_name, unknown_hits, known, all_by_prog = {}, {}, {}, {}
    for h in passing:
        prog = TRACK[h["track"]]["program"]
        e = h["evalue"] if h["evalue"] is not None else math.inf
        iv = (h["q0"], h["q1"])
        all_by_prog.setdefault(prog, {})
        all_by_prog[prog][iv] = min(e, all_by_prog[prog].get(iv, math.inf))
        if not h["_names"]:
            unknown_hits.setdefault(prog, {})[iv] = min(e, unknown_hits.get(prog, {}).get(iv, math.inf))
            continue
        ke, kn = known.setdefault(prog, {}).get(iv, (math.inf, set()))
        known[prog][iv] = (min(e, ke), kn | h["_names"])
        for n in h["_names"]:
            d = by_name.setdefault(n, {}).setdefault(prog, {})
            d[iv] = min(e, d.get(iv, math.inf))

    def best_pair(xa, xb, ranges=False):
        """Lowest E-value of an overlapping pair, one range from each dict, or None. With
        ranges=True, returns (E-value, range from xa, range from xb)."""
        cand = sorted([(e, iv, xb, 0) for iv, e in xa.items()] + [(e, iv, xa, 1) for iv, e in xb.items()],
                      key=lambda c: (c[0], c[1]))
        for e, (a0, a1), other, side in cand:
            for (b0, b1) in other:
                if overlaps(a0, a1, b0, b1):
                    if not ranges:
                        return e
                    return (e, (a0, a1), (b0, b1)) if side == 0 else (e, (b0, b1), (a0, a1))
        return None

    support = {}
    for n, per_prog in by_name.items():
        fs = list(per_prog)
        for i, fa in enumerate(fs):
            for fb in fs[i + 1:]:
                e = best_pair(per_prog[fa], per_prog[fb])
                if e is not None:
                    got, best = support.get(n, (set(), math.inf))
                    support[n] = (got | {fa, fb}, min(best, e))
    if support:
        # The Pfam family named by the most programs. On a tie, the one whose own domain on
        # this protein (an hmmscan or InterProScan hit) has the best E-value, then the best
        # E-value of an agreeing pair.
        own = {}
        for h in passing:
            if h["track"] in ("hmmscan", "interproscan") and h["_names"]:
                for n in h["_names"]:
                    own[n] = min(own.get(n, math.inf), h["evalue"] if h["evalue"] is not None else math.inf)
        order = sorted(support, key=lambda n: (-len(support[n][0]), own.get(n, math.inf), support[n][1], n))
        AGREE_OTHERS[0] = order[1:]
        return progs, "agree", order[0]
    # No family is shared. Overlapping hits that both name a family: 'disagree'. The first
    # overlapping pair (best E-value first) of each two programs gives the families listed.
    named, kp = {}, sorted(known, key=PROGRAM_ORDER.index)
    for i, fa in enumerate(kp):
        for fb in kp[i + 1:]:
            ea = {iv: v[0] for iv, v in known[fa].items()}
            eb = {iv: v[0] for iv, v in known[fb].items()}
            got = best_pair(ea, eb, ranges=True)
            if got:
                e, ra, rb = got
                for prog, n_set in ((fa, known[fa][ra][1]), (fb, known[fb][rb][1])):
                    for n in sorted(n_set):
                        best, ps = named.get(n, (math.inf, []))
                        named[n] = (min(best, e), ps + [prog] if prog not in ps else ps)
    if named:
        DISAGREE_NAMES[0] = [(n, sorted(named[n][1], key=PROGRAM_ORDER.index))
                             for n in sorted(named, key=lambda n: (named[n][0], n))]
        return progs, "disagree", None
    for fa, xa in unknown_hits.items():
        for fb, xb in all_by_prog.items():
            if fb != fa and best_pair(xa, xb) is not None:
                return progs, "overlap_unknown", None
    return progs, "apart", None


def coverage(hits, length):
    """Number of programs with a counted hit over each residue (list of length `length`)."""
    per_prog = {}
    for h in hits:
        if counts(h):
            per_prog.setdefault(TRACK[h["track"]]["program"], []).append((h["q0"], h["q1"]))
    diff = [0] * (length + 2)
    for ivs in per_prog.values():
        ivs = sorted((max(1, a), min(length, b)) for a, b in ivs if b >= a)
        merged = []
        for a, b in ivs:
            if merged and a <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        for a, b in merged:
            diff[a] += 1
            diff[b + 1] -= 1
    out, c = [], 0
    for i in range(1, length + 1):
        c += diff[i]
        out.append(c)
    return out


# ------------------------------------------------------------------------ alignment ----

_ALIGNER = None


def aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio.Align import PairwiseAligner, substitution_matrices
        a = PairwiseAligner()
        a.mode = "local"
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score = -11
        a.extend_gap_score = -1
        _ALIGNER = a
    return _ALIGNER


def clean(s):
    return re.sub(r"[^A-Z]", "X", s.upper())


def recompute(qseq, q0, tseq, t0):
    """Local alignment of two segments; returns gapped strings and 1-based starts."""
    if not qseq or not tseq:
        return None
    try:
        alns = aligner().align(clean(qseq), clean(tseq))
        a = alns[0]
    except (ValueError, IndexError, OverflowError):
        return None
    # a.aligned: the matched blocks of each sequence; what lies between two blocks is a gap
    # in one sequence (Biopython 1.79 and later).
    qblocks, tblocks = a.aligned
    qa, ta = [], []
    for k, ((i0, i1), (j0, j1)) in enumerate(zip(qblocks, tblocks)):
        if k:
            pi, pj = qblocks[k - 1][1], tblocks[k - 1][1]
            if i0 > pi:
                qa.append(qseq[pi:i0])
                ta.append("-" * (i0 - pi))
            if j0 > pj:
                qa.append("-" * (j0 - pj))
                ta.append(tseq[pj:j0])
        qa.append(qseq[i0:i1])
        ta.append(tseq[j0:j1])
    c0q, c0t = int(qblocks[0][0]), int(tblocks[0][0])
    return dict(q="".join(qa), t="".join(ta), qs=q0 + c0q, ts=t0 + c0t)


# ------------------------------------------------------------------------- drawing ----

def esc(s):
    return html.escape(str(s), quote=True)


def fmt_e(e):
    if e is None:
        return "none"
    if e == 0:
        return "0"
    return "%.1e" % e if e < 0.01 else "%.3g" % e


def text_w(s, px=11):
    """Width of a box label in the monospace font, measured: 0.58 of the font size per character."""
    return len(s) * px * 0.58


LABEL_GAP = 8  # px between a box and a label placed beside it


def place_labels(hits, px, W):
    """Where each box's label goes, and the lane (row inside the tool's row) of each box.
    A label that fits inside its box goes inside. Otherwise it goes to the right of the box,
    joined to it by a short leader line, or to the left when the right would leave the drawing;
    the font is never made smaller and the name is never cut. Boxes stack into lanes by the
    room they take with their label, so no two labels or boxes overlap. Returns, per hit,
    (lane, side, label x) with side 'in', 'right' or 'left'; and the number of lanes."""
    lanes, out = [], []
    for h in hits:
        x0, x1 = (h["q0"] - 1) * px, h["q1"] * px
        w = max(2.0, x1 - x0)
        lw = text_w(h["label"])
        if lw + 6 <= w:
            side, lx, a, b = "in", x0 + 3, x0, x0 + w
        elif x1 + LABEL_GAP + lw <= W or x0 - LABEL_GAP - lw < 0:
            lx = min(x1 + LABEL_GAP, W - lw)
            side, a, b = "right", x0, max(x0 + w, lx + lw)
        else:
            lx = x0 - LABEL_GAP - lw
            side, a, b = "left", lx, x0 + w
        for li, ends in enumerate(lanes):
            if all(a > e1 + 4 or b < e0 - 4 for e0, e1 in ends):
                ends.append((a, b))
                break
        else:
            lanes.append([(a, b)])
            li = len(lanes) - 1
        out.append((li, side, lx))
    return out, max(1, len(lanes))


LANE_H, BOX_H, ROW_PAD = 18, 14, 8


def shape_svg(x, y, w, h, prog, extra=""):
    """A program's box: its kind's colour, its own shape (filled, open, hatched, rounded)."""
    col, shape = prog["color"], prog["shape"]
    if shape == "open":
        return ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="2" fill="#fff" stroke="%s" stroke-width="2"%s/>'
                % (x + 1, y + 1, max(1.0, w - 2), h - 2, col, extra))
    if shape == "hatched":
        return ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="2" fill="url(#hatch-%s)" stroke="%s" '
                'stroke-width="1"%s/>' % (x + 0.5, y + 0.5, max(1.0, w - 1), h - 1, prog["cat"], col, extra))
    rx = h / 2 if shape == "rounded" else 2
    return '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%.1f" fill="%s"%s/>' % (x, y, w, h, rx, col, extra)


def e_opacity(e):
    """E-value as box opacity, on a log scale: 0.25 at E-value 10, full at 1e-30 or below (and
    for an E-value of 0 or none). The eye ranks hits by strength without reading tooltips."""
    if e is None or e <= 0:
        return 1.0
    return round(max(0.25, min(1.0, 0.25 + 0.75 * (1 - math.log10(e)) / 31)), 2)


def hatch_defs():
    """The hatch patterns, one per colour, in a hidden SVG that every drawing on the page uses."""
    pats = "".join('<pattern id="hatch-%s" width="4" height="4" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
                   '<rect width="4" height="4" fill="#fff"/><line x1="1" y1="0" x2="1" y2="4" stroke="%s" stroke-width="1.6"/>'
                   '</pattern>' % (k, c["color"]) for k, c in CATEGORIES.items())
    return '<svg width="0" height="0" style="position:absolute" aria-hidden="true"><defs>%s</defs></svg>' % pats


def box_svg(h, idx, x0, x1, y, track, side="in", lx=None):
    prog = PROGRAMS[track["program"]]
    w = max(2.0, x1 - x0)
    op = e_opacity(h["evalue"])
    rect = shape_svg(x0, y, w, BOX_H, prog, ' opacity="%.2f"' % op if op < 1 else "")
    ty = y + BOX_H - 3
    if side == "in":
        solid = prog["shape"] in ("filled", "rounded") and op >= 0.6  # white text needs a dark box
        tcol, cls = ("#ffffff", "bl") if solid else ("#1a1a1a", "bl halo")
        lab = '<text x="%.1f" y="%d" fill="%s" class="%s">%s</text>' % (x0 + 3, ty, tcol, cls, esc(h["label"]))
    else:  # beside the box, with a leader line in the box's colour
        lw = text_w(h["label"])
        l0, l1 = (x0 + w, lx - 1) if side == "right" else (lx + lw + 1, x0)
        lead = ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>'
                % (l0, y + BOX_H / 2, l1, y + BOX_H / 2, prog["color"])) if l1 > l0 else ""
        lab = lead + '<text x="%.1f" y="%d" fill="#1a1a1a" class="bl out">%s</text>' % (lx, ty, esc(h["label"]))
    return '<g class="hit" tabindex="0" data-i="%d">%s%s</g>' % (idx, rect, lab)


def page_svg(p, hits_by_track, cov, plddt, W, px, run_status, m=8):
    """The track area as two SVGs (names on the left, the drawing on the right), with the best
    TOP_N hits of each tool. Each row is a <g class="row"> in its own coordinates, so the
    page's script can redraw a tool's row with all of its hits ("Show all")."""
    L = p["length"]
    rows = []  # (key, height, left svg, right svg)
    r = ['<line x1="0" y1="18" x2="%.1f" y2="18" stroke="%s" stroke-width="2"/>' % (W, PROTEIN_LINE)]
    for pos in range(50, L + 1, 50):
        x = (pos - 0.5) * px
        r.append('<line x1="%.1f" y1="18" x2="%.1f" y2="%d" stroke="%s"/>' % (x, x, 22 if pos % 100 else 24, PROTEIN_LINE))
        if pos % 100 == 0:
            r.append('<text x="%.1f" y="10" class="rl" text-anchor="middle">%d</text>' % (x, pos))
    r.append('<text x="1" y="10" class="rl">1</text>')
    r.append('<text x="%.1f" y="30" class="rl" text-anchor="end">%d aa</text>' % (W - 1, L))
    rows.append(("ruler", 32, '<text x="4" y="20" class="tn">residue</text>', "".join(r)))
    for t in TRACKS:
        k = t["key"]
        hs = hits_by_track.get(k, [])
        shown = hs[:TOP_N]
        place, n_l = place_labels(shown, px, W)
        h_row = ROW_PAD + n_l * LANE_H + 4
        if run_status[k] == "not run":
            g = ['<text x="4" y="%d" class="nr">not run: no %s in the input folder</text>'
                 % (ROW_PAD + 11, esc(t["file"]))]
        elif not hs:
            g = ['<text x="4" y="%d" class="nh">no hit with E-value at most %s</text>' % (ROW_PAD + 11, fmt_thr(COUNT_E))]
        else:
            g = ['<line x1="0" y1="%d" x2="%.1f" y2="%d" stroke="%s" stroke-width="1.5"/>'
                 % (ROW_PAD + BOX_H // 2, W, ROW_PAD + BOX_H // 2, PROTEIN_LINE)]
        for h, (li, side, lx) in zip(shown, place):
            g.append(box_svg(h, h["i"], (h["q0"] - 1) * px, h["q1"] * px, ROW_PAD + li * LANE_H, t, side, lx))
        more = len(hs) - len(shown)
        if more > 0:
            h_row = max(h_row, ROW_PAD + 34)  # room for the "+N more" line under the name
        extra = ('<text x="26" y="%d" class="tm more">+%d more hit%s</text>' % (ROW_PAD + 26, more, "" if more == 1 else "s")) if more > 0 else ""
        left = ('<g transform="translate(4,%d)">%s</g><text x="26" y="%d" class="tn">%s</text>%s'
                % (ROW_PAD, legend_mark(t), ROW_PAD + 11, esc(t["name"]), extra))
        rows.append((k, h_row, left, "".join(g)))
    if plddt is not None:
        g = []
        for i, v in enumerate(plddt[:L]):
            c = next((col for lo, hi, col, _ in PLDDT if lo <= v < hi), PLDDT[-1][2])
            g.append('<rect x="%.2f" y="%d" width="%.2f" height="%d" fill="%s"/>' % (i * px, ROW_PAD, px + 0.05, BOX_H, c))
        rows.append(("plddt", ROW_PAD + BOX_H + 8, '<text x="4" y="%d" class="tn">pLDDT</text>' % (ROW_PAD + 11), "".join(g)))
    # programs with a hit: a step line from 0 to m, the programs run; 48 px tall whatever m is
    ah, top = 48, 10
    step_h = ah / max(1, m)
    pts, prev = [], None
    for i, c in enumerate(cov):
        if c != prev:
            if prev is not None:
                pts.append("L%.2f,%.1f" % (i * px, top + ah - prev * step_h))
            pts.append(("M" if prev is None else "L") + "%.2f,%.1f" % (i * px, top + ah - c * step_h))
            prev = c
    if prev is not None:
        pts.append("L%.2f,%.1f" % (L * px, top + ah - prev * step_h))
    ticks = "".join('<line x1="0" y1="%.1f" x2="5" y2="%.1f" stroke="%s"/>' % (top + ah - v * step_h, top + ah - v * step_h, AGREEMENT)
                    for v in (0, m))
    right = ticks + '<path d="%s" fill="none" stroke="%s" stroke-width="1.6"/>' % ("".join(pts), AGREEMENT)
    left = '<text x="4" y="%d" class="tn">programs with a hit</text>' % (top + 12) + "".join(
        '<text x="100%%" dx="-4" y="%.1f" class="rl" text-anchor="end">%s</text>' % (top + ah - v * step_h + 3, lab)
        for v, lab in ((0, "0"), (m, "%d programs" % m)))
    rows.append(("agree", top + ah + 8, left, right))
    y, L_out, R_out = 0, [], []
    for k, hgt, left, right in rows:
        sep_l = '<line x1="0" y1="0" x2="100%" y2="0" class="sep"/>' if y else ""
        sep_r = '<line x1="0" y1="0" x2="%.1f" y2="0" class="sep"/>' % W if y else ""
        L_out.append('<g class="row" data-k="%s" data-h="%d" transform="translate(0,%d)">%s%s</g>' % (k, hgt, y, sep_l, left))
        R_out.append('<g class="row" data-k="%s" data-h="%d" transform="translate(0,%d)">%s%s</g>' % (k, hgt, y, sep_r, right))
        y += hgt
    return y, "".join(L_out), "".join(R_out)


def legend_mark(t, w=16, h=BOX_H):
    return shape_svg(0, 0, w, h, PROGRAMS[t["program"]])


def evalue_key():
    """Boxes at four E-values, as drawn: fainter is weaker."""
    marks = "".join('<span class="kk">%s%s</span>' % (svg_inline(
        '<rect width="12" height="10" rx="2" fill="#1a1a1a" opacity="%.2f"/>' % e_opacity(e), 12, 10), lab)
        for e, lab in ((10, "10"), (1e-3, "0.001"), (1e-10, "1e-10"), (1e-30, "1e-30 or less")))
    return '<p class="kind">Fainter box: higher (weaker) E-value. %s</p>' % marks


def kind_key():
    """The two lines that open the legend: colour = kind of program, shape = program."""
    sw = "".join('<span class="kk"><i style="background:%s"></i>%s</span>' % (c["color"], esc(c["name"]))
                 for c in CATEGORIES.values())
    return ('<p class="kind">Colour: the kind of program. %s</p><p class="kind">Shape: which program, inside a '
            'kind: filled, open, hatched or rounded.</p>' % sw) + evalue_key()


def svg_inline(inner, w, h, cls=""):
    vb = "" if isinstance(w, str) else ' viewBox="0 0 %s %d"' % (w, h)
    return '<svg class="%s" width="%s" height="%d"%s xmlns="http://www.w3.org/2000/svg">%s</svg>' % (
        cls, w, h, vb, inner)


def sparkline(cov, length, m=8, w=90, h=16):
    """The "programs with a hit" line squeezed to w pixels, 0 to m programs high, the same
    scale for every protein. data-c holds the highest count in each pixel, for the hover."""
    if not cov:
        return svg_inline("", w, h)
    bins = []
    for b in range(w):
        a0 = int(b * length / w)
        a1 = max(a0 + 1, int((b + 1) * length / w))
        bins.append(max(cov[a0:a1]))
    d = []
    for b, c in enumerate(bins):
        yy = h - 1 - c * (h - 2) / max(1, m)
        d.append(("M" if b == 0 else "L") + "%d,%.1f L%d,%.1f" % (b, yy, b + 1, yy))
    base = '<line x1="0" y1="%d" x2="%d" y2="%d" stroke="%s"/>' % (h - 0.5, w, h - 0.5, PROTEIN_LINE)
    svg = svg_inline(base + '<path d="%s" fill="none" stroke="%s" stroke-width="1.2"/>' % (" ".join(d), AGREEMENT), w, h, "spark")
    return svg.replace("<svg ", '<svg data-c="%s" data-l="%d" data-m="%d" ' % ("".join("%x" % min(c, 15) for c in bins), length, m), 1)


# --------------------------------------------------------------------------- pages ----

CSS_LINKS = (
    '<link rel="preconnect" href="https://fonts.googleapis.com">'
    '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
    '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600'
    '&family=Source+Sans+3:wght@400;600;700&display=swap">'
)


ASSET_VER = [""]  # "?v=<hash of assets/>": a rebuilt site makes browsers fetch the new scripts


def page_head(title, depth):
    pre = "../" * depth
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            '<title>%s</title>%s<link rel="stylesheet" href="%sassets/site.css%s"></head><body>%s'
            % (esc(title), CSS_LINKS, pre, ASSET_VER[0], hatch_defs()))


def legend_html(ctx, has_plddt, has_expr):
    """The legend, in a narrow column left of the drawing (above it on a phone), so the reader
    meets the key before the marks: one line per tool, its mark then its name. The rule for a
    hit, the setting, the database and the columns that give the box open on a click."""
    items = []
    for t in TRACKS:
        k = t["key"]
        prog = PROGRAMS[t["program"]]
        mark = svg_inline(legend_mark(t, 12, 10), 12, 10)
        if ctx["run_status"][k] == "not run":
            items.append('<li class="off">%s<span>%s: not run</span></li>' % (mark, esc(t["name"])))
            continue
        c, d = fmt_thr(COUNT_E), fmt_thr(E_MAX)
        if k == "hhsearch":
            thr = "counts at Prob ≥ 50 and E-value at most %s; drawn up to E-value %s" % (c, d)
        elif k == "interproscan":
            thr = ("counts at E-value at most %s where there is one, drawn up to %s; coils, disorder, membrane "
                   "and signal-peptide predictions do not count" % (c, d))
        else:
            thr = "counts at %s at most %s; drawn up to %s" % (t["score"], c, d)
        coords = {"domtbl_seq": "hmm from-to (positions on your protein)",
                  "domtbl_pfam": "env from-to, the domain's envelope: the stretch of your protein hmmscan assigns "
                                 "to the domain, a little wider than the aligned part",
                  "m8": "qstart-qend", "kmerseek": "region_start-region_end, plus 1 at the start",
                  "hhr": "Query HMM range", "interproscan": "start-stop"}[t["kind"]]
        lb = ctx["leaderboard"].get(k)
        detail = [thr, ctx["settings"][k], ctx["databases"][k], "box: " + coords] + ([lb] if lb else [])
        items.append('<li><details><summary>%s<span>%s</span></summary><div class="sub">%s</div></details></li>' % (
            mark, esc(t["name"]), "<br>".join(esc(x) for x in detail)))
    line = svg_inline('<line x1="0" y1="5" x2="12" y2="5" stroke="%s" stroke-width="2"/>' % PROTEIN_LINE, 12, 10)
    items.append('<li>%s<span>your protein, residue 1 to the end</span></li>' % line)
    if has_plddt:
        sw = "".join('<span class="sw" style="background:%s"></span>%s ' % (c, esc(lbl)) for _, _, c, lbl in PLDDT)
        items.append('<li><span></span><span><b>pLDDT</b>, AlphaFold\'s confidence in each residue, 0 to 100: %s</span></li>' % sw)
    step = svg_inline('<path d="M0,9 L4,9 L4,2 L8,2 L8,6 L12,6" fill="none" stroke="%s" stroke-width="1.4"/>' % AGREEMENT, 12, 10)
    items.append('<li>%s<span>programs with a hit at each residue, 0 to %d (the programs run), each program counted once</span></li>'
                 % (step, ctx["n_run"]))
    if has_expr:
        bar = svg_inline('<rect x="0.5" y="2.5" width="11" height="5" fill="none" stroke="%s"/>' % EXPRESSION, 12, 10)
        items.append('<li>%s<span>expression of this protein\'s gene in each sample, from expression.tsv</span></li>' % bar)
    return ('<aside class="side legend"><h2>Legend</h2>%s<p class="note">Click a tool for its setting, database '
            'and what makes a hit.</p><ul>%s</ul></aside>' % (kind_key(), "".join(items)))


def of_programs(m):
    """'of 8 programs' when every program was run, else 'of the 5 programs run'."""
    return "of 8 programs" if m == 8 else "of the %d programs run" % m


def header_sentence(n, status, agreed, names, m=8):
    """The sentence at the top of a protein page. m: programs with at least one input file."""
    if n == 0:
        return "No tool reports a hit." if m == 8 else "None of the %d programs run reports a hit." % m
    s = "%d %s report%s a hit." % (n, of_programs(m), "s" if n == 1 else "")
    if n == 1:
        return s
    if status == "agree":
        return s + " They agree on %s." % family_label(agreed, names)
    if status == "overlap_unknown":
        return s + " Their hits overlap, but no Pfam family is known for the Swiss-Prot proteins they match."
    if status == "disagree":
        fams = DISAGREE_NAMES[0]
        # "Pkinase PF00069 (phmmer)": name and accession without a second pair of brackets
        listed = ", ".join("%s (%s)" % (" ".join(dict.fromkeys([names.get(n) or n, n])),
                                        ", ".join(PROGRAMS[p]["name"] for p in ps))

                           for n, ps in fams[:6])
        more = " and %d more" % (len(fams) - 6) if len(fams) > 6 else ""
        return s + " Their hits overlap but name different Pfam families: %s%s." % (listed, more)
    return s + " Their hits are on different parts of the protein and do not overlap."



def family_label(n, names):
    nm = names.get(n)
    return "%s (%s)" % (nm, n) if nm and nm != n else n


def family_full(n, names):
    """family_label plus Pfam's one-line description, e.g.
    'MHC_I (PF00129): Class I Histocompatibility antigen, domains alpha 1 and 2'."""
    lab, d = family_label(n, names), PFAM_DESC.get(n)
    return "%s: %s" % (lab, d) if d and d != names.get(n) else lab


def swissprot_header(h):
    """(entry name, description, gene, organism) from a UniProt FASTA header such as
    'sp|P04439|HLAA_HUMAN HLA class I histocompatibility antigen, A alpha chain OS=Homo sapiens
    OX=9606 GN=HLA-A PE=1 SV=2'. Missing parts are ''."""
    first, _, rest = h.partition(" ")
    parts = first.split("|")
    entry = parts[2] if len(parts) > 2 else ""
    desc = re.split(r" (?:OS|OX|GN|PE|SV)=", " " + rest, maxsplit=1)[0].strip()
    org = re.search(r" OS=(.*?)(?= (?:OX|GN|PE|SV)=|$)", " " + rest)
    gene = re.search(r" GN=(\S+)", " " + rest)
    return entry, desc, gene.group(1) if gene else "", org.group(1).strip() if org else ""


def hmmer_rows(h, seq):
    """HMMER's printed rows as an alignment of your protein against the other sequence.
    phmmer/jackhmmer: the model was built from your protein, so its row is written out with
    your protein's residues at the model positions ('.' in the model row: a residue the target
    has and the model does not, a gap in your protein). hmmscan: your protein is the target
    row, and the other row is the Pfam model's consensus ('.' becomes a gap)."""
    ms, ts, mf, tf = h["raw_aln"]
    tool = TRACK[h["track"]]["name"]
    if h["track"] == "hmmscan":
        return dict(q=ts, t=ms.replace(".", "-"), qs=tf, ts=mf,
                    src="the alignment hmmscan wrote (hmmscan.out); the lower row is the Pfam model's consensus: "
                        "capital letters where the family keeps one residue, small letters where it varies more")
    q, pos = [], mf
    for c in ms:
        if c == ".":
            q.append("-")
        else:
            q.append(seq[pos - 1] if pos - 1 < len(seq) else "X")
            pos += 1
    return dict(q="".join(q), t=ts.replace(".", "-"), qs=mf, ts=tf,
                src="the alignment %s wrote (%s.out), with its model row written as your protein's residues"
                    % (tool, h["track"]))


def non_overlapping(doms, label):
    """Domains to draw on a small map: in the given order, skipping any that overlap one kept."""
    kept = []
    for d0, d1, nm in doms:
        if all(d1 < a or d0 > b for a, b, _ in kept):
            kept.append((d0, d1, label(nm)))
    return [list(k) for k in sorted(kept)]


def id_gene_name(pid):
    """The gene name a gene model's ID carries after its number (RFC5 in FUN000798_RFC5,
    TTN in FUN019272_TTN.1), or None for an ID with no name (FUN017965_FUN017965)."""
    head, _, rest = pid.partition("_")
    rest = re.sub(r"\.\d+$", "", rest)
    if not rest or rest == head or re.fullmatch(r"FUN\d+", rest):
        return None
    return rest


def build_protein(task):
    """Runs in a worker: draws one protein page, returns its row for the index."""
    pid, ctx_key = task
    ctx = CTX
    p = ctx["proteins"][pid]
    hits = [materialize(x) for x in ctx["hits"].get(pid, [])]
    L = p["length"]
    seq = p["seq"]
    for h in hits:
        h["_names"] = target_names(h, ctx["target_ann"])
    progs, status, agreed = agreement(hits, ctx["target_ann"])
    cov = coverage(hits, L)
    best_e = min([h["evalue"] for h in hits if h["passes"] and h["evalue"] is not None], default=None)

    # every hit with E-value at most E_MAX goes into the page's data; per track, best E-value first.
    # The drawing as built shows the hits that count (E-value at most COUNT_E); the slider shows more.
    by_track = {}
    data = []
    score_names = {}
    terms = {}  # text the index page's search box looks through -> bit mask of the programs whose hits carry it
    pname = ctx["prot_names"].get(pid)
    if pname:
        terms[pname] = 0  # 0: the protein's own name from protein_names.tsv, not a hit
    gname = id_gene_name(pid)
    gmask, gexample = 0, None  # programs with a hit on a Swiss-Prot protein of that gene
    alns = {}  # hit number -> alignment rows, written next to the page and read on a click
    for t in TRACKS:
        hs = sorted([h for h in hits if h["track"] == t["key"] and h["passes"]],
                    key=lambda h: (h["evalue"] if h["evalue"] is not None else 1e300, h["q0"]))
        for rank, h in enumerate(hs):
            h["i"] = len(data)
            name = ctx["entry_names"].get(h["target"], h["target"]) if t["key"] not in ("hmmscan", "interproscan") else h["target"]
            h["label"] = name
            d = dict(k=t["key"], l=name, q=[h["q0"], h["q1"]], s=[v for _, v in h["scores"]], e=h["evalue"])
            # For the tooltip and the search box: the full name, gene and species of a Swiss-Prot
            # target, and the Pfam families on the matched part of it (or the family itself, for
            # hmmscan and InterProScan). f is None when no Pfam output covers the target.
            if t["key"] == "hmmscan":
                pf = next(iter(h["names"]))
                if PFAM_DESC.get(pf):
                    d["n"] = PFAM_DESC[pf]
            elif t["key"] != "interproscan":
                info = ctx["target_info"].get(h["target"])
                if info:
                    d["n"], d["g"], d["o"] = info
            d["f"] = None if h["_names"] is None else [family_full(n, ctx["names"]) for n in sorted(h["_names"])]
            bit = 1 << PROGRAM_ORDER.index(t["program"])
            for txt in [name, h["target"], d.get("n"), d.get("g")] + (d["f"] or []):
                if txt:
                    terms[txt] = terms.get(txt, 0) | bit
            if gname and counts(h):
                entry_gene = name.split("_")[0] if "_" in name else ""
                if gname.upper() in ((d.get("g") or "").upper(), entry_gene.upper()):
                    gmask |= bit
                    gexample = gexample or name
            score_names.setdefault(t["key"], [n for n, _ in h["scores"]])
            a = None
            if t["key"] == "hmmscan":
                d["r"] = h["hmm"][:2]
                d["tl"] = h["hmm"][2]
                if h.get("raw_aln"):
                    a = hmmer_rows(h, seq)
            elif t["key"] != "interproscan":
                d["t"] = h["target"]
                d["r"] = [h["t0"], h["t1"]]
                tseq = ctx["targets"].get(h["target"])
                if tseq:
                    d["tl"] = len(tseq)
                doms = ctx["target_ann"].get(h["target"])
                if doms:
                    d["td"] = non_overlapping(doms, lambda nm: family_label(sorted(nm)[0], ctx["names"]) if nm else "")
                if t["key"] == "hhsearch":
                    a = dict(h["aln"], src="the alignment hhsearch wrote (the .hhr file)") if h.get("aln") else None
                elif t["key"] == "kmerseek":
                    d["m"] = h.get("moltype")
                elif h.get("aln"):
                    a = h["aln"]
                elif h.get("raw_aln"):
                    a = hmmer_rows(h, seq)
                elif tseq is not None and rank < TOP_N:
                    lo, hi = min(h["t0"], h["t1"]), max(h["t0"], h["t1"])
                    a = recompute(seq[h["q0"] - 1:h["q1"]], h["q0"], tseq[lo - 1:hi], lo)
                    if a:
                        a["src"] = ("recomputed when the site was built, because %s was not given: the two ranges "
                                    "aligned again with Biopython's PairwiseAligner (local, BLOSUM62, gap open -11, "
                                    "gap extend -1)" % ("%s.out" % t["key"] if t["kind"] == "domtbl_seq" else "a qaln/taln column"))
            if a:
                alns[h["i"]] = a
                d["a"] = 1
            data.append(d)
        by_track[t["key"]] = hs

    plddt = ctx["plddt"].get(pid)
    W = int(max(640, min(L * 2.0, 3200)))
    px = W / L
    counted_by_track = {k: [h for h in v if counts(h)] for k, v in by_track.items()}
    hgt, left, right = page_svg(p, counted_by_track, cov, plddt, W, px, ctx["run_status"], ctx["n_run"])
    blocks = ['<div class="tracks" id="tracks"><div class="names">%s</div><div class="scroll">%s</div></div>'
              % (svg_inline(left, "100%", hgt), svg_inline(right, W, hgt, "trk"))]
    n_more = sum(max(0, len(v) - TOP_N) for v in counted_by_track.values())

    expr = ""
    if ctx["expression"]:
        samples, rows = ctx["expression"]
        vals = rows.get(pid)
        if vals:
            mx = max([v for v in vals if v is not None] + [1e-9])
            items = []
            for s, v in zip(samples, vals):
                wpct = 0 if v is None else 100 * v / mx
                items.append('<div class="er"><span class="es">%s</span><span class="eb"><span style="width:%.1f%%"></span></span><span class="ev">%s</span></div>'
                             % (esc(s), wpct, "missing" if v is None else "%g" % v))
            expr = '<section><h2>Expression</h2><p class="note">Values as in expression.tsv (count or TPM, whichever the file holds); bar length is relative to this protein\'s highest sample.</p>%s</section>' % "".join(items)

    sentence = header_sentence(len(progs), status, agreed, ctx["names"], ctx["n_run"])
    others = AGREE_OTHERS[0]
    other_note = ""
    if others:
        shown_o = [family_label(n, ctx["names"]) for n in others[:8]]
        other_note = ('<p class="note">%d other Pfam famil%s also meet%s the rule: %s%s. The one named has hits '
                      'from the most programs; on a tie, the best E-value of its own hmmscan or InterProScan '
                      'domain on this protein.</p>' % (len(others), "y" if len(others) == 1 else "ies",
                                                       "s" if len(others) == 1 else "", esc(", ".join(shown_o)),
                                                       " and %d more" % (len(others) - 8) if len(others) > 8 else ""))
    gnote = ""
    if gname:
        gps = [PROGRAMS[f]["name"] for b, f in enumerate(PROGRAM_ORDER) if gmask & (1 << b)]
        gnote = ('<p class="note">The ID names the gene %s. %s</p>' % (esc(gname), (
            "Hits on a Swiss-Prot protein of that gene (gene name or entry name %s_...): %s; for example %s."
            % (esc(gname), esc(", ".join(gps)), esc(gexample))) if gps else
            "No hit is on a Swiss-Prot protein with that gene name or entry name %s_..." % esc(gname)))
    out = [page_head("%s · tool comparison" % pid, 1)]
    out.append('<header class="wide"><a class="back" href="../index.html">All proteins</a><h1 class="mono">%s</h1>'
               '<p class="len">%s%d amino acids</p><p class="summary">%s</p>%s%s</header><main class="wide">' % (
        esc(pid), (esc(pname) + " · ") if pname else "", L, esc(sentence), other_note, gnote))
    # E-value slider, in steps of 10 times: from 10 (every drawn hit) down to the smallest
    # E-value above 0 on this page. It starts at the counting threshold.
    pos_e = [h["evalue"] for v in by_track.values() for h in v if h["evalue"]]
    c_pos = int(round(math.log10(COUNT_E)))
    lo = min(c_pos, max(-300, math.floor(math.log10(min(pos_e))))) if pos_e else c_pos
    slider = ('<div class="eslider"><label for="emax">Show hits with E-value at most <output id="ev">%s</output></label>'
              '<input id="emax" type="range" min="%d" max="1" step="1" value="%d" data-start="%d">'
              '<span class="note">The sentence at the top and the index count hits at E-value at most %s. Move the '
              'slider up to see weaker hits, down to see fewer.</span></div>'
              % (fmt_thr(COUNT_E), lo, c_pos, c_pos, fmt_thr(COUNT_E))) if by_track else ""
    toggle = ""
    if n_more:
        toggle = ('<button id="toggle" type="button">Show all %d hits (%d more)</button>'
                  % (sum(len(v) for v in counted_by_track.values()), n_more))
    # The drawing comes first; each of its rows starts with the tool's name and mark, and the
    # table under it gives each tool's setting, database and threshold.
    out.append('<div class="withside">' + legend_html(ctx, plddt is not None, bool(ctx["expression"] and ctx["expression"][1].get(pid))))
    out.append('<section><h2>What each tool found</h2><p class="lead">Each box is one hit, at the position the tool '
               'reported on your protein.</p>%s%s%s<p class="note">A box is drawn for each hit with '
               'E-value at most the slider\'s value; a Swiss-Prot protein matched in two places gives two boxes. Up to %d hits per tool are drawn, best '
               'E-value first (ties: the leftmost first); "+N more" under a tool\'s name counts the rest. Hover or '
               'tap a box for the protein\'s full name, gene, species and Pfam families, and the scores; click it '
               'to see the residues. A long protein\'s drawing scrolls sideways.</p></section>'
               % (slider, "".join(blocks), toggle, TOP_N))
    out.append('</div>')
    out.append(expr)
    out.append('<div id="tip" role="tooltip" hidden></div><aside id="panel" hidden><button id="close" type="button">Close</button><div id="pbody"></div></aside>')
    out.append('</main>%s' % ctx["footer"])
    page_data = dict(
        id=pid, q=seq, px=px, w=W, top=TOP_N, m=ctx["n_run"], hits=data, names=score_names, alpha=ctx["alphabets"],
        tracks={t["key"]: dict(name=t["name"], color=PROGRAMS[t["program"]]["color"], shape=PROGRAMS[t["program"]]["shape"],
                               cat=PROGRAMS[t["program"]]["cat"], p=t["program"]) for t in TRACKS},
        pfam_color=CATEGORIES["pfam"]["color"],
        tseq_files=ctx["n_shards"])
    out.append('<script>const PAGE=%s;</script><script src="../assets/protein.js%s"></script></body></html>'
               % (json.dumps(page_data, separators=(",", ":")).replace("</", "<\\/"), ASSET_VER[0]))
    (ctx["out"] / "proteins" / (safe(pid) + ".html")).write_text("".join(out))
    (ctx["out"] / "proteins" / (safe(pid) + ".aln.js")).write_text(
        "window.ALN=%s;" % json.dumps({str(k): v for k, v in alns.items()}, separators=(",", ":")))
    return dict(id=pid, length=L, n=len(progs), progs=progs, status=status,
                agreed=family_label(agreed, ctx["names"]) if agreed else None,
                best_e=best_e, spark=sparkline(cov, L, ctx["n_run"]), terms=terms, name=pname,
                gname=gname, gprogs=[PROGRAMS[f]["name"] for b, f in enumerate(PROGRAM_ORDER) if gmask & (1 << b)],
                gexample=gexample)


def shard_of(acc, n):
    """Which targets/NNN.js file holds a Swiss-Prot sequence. protein.js computes the same."""
    h = 0
    for ch in acc:
        h = (h * 31 + ord(ch)) & 0xFFFFFFFF
    return h % n


def write_target_shards(out, hits, targets):
    """Swiss-Prot sequences of every drawn search hit, split into files a page loads only when
    a box is clicked, so a sequence is stored once for the whole site, not once per page."""
    used = {t for t in USED_TARGETS if t in targets}
    n = max(1, min(256, len(used) // 200 + 1))
    shards = {}
    for acc in used:
        shards.setdefault(shard_of(acc, n), {})[acc] = targets[acc]
    (out / "targets").mkdir()
    for i in range(n):
        (out / "targets" / ("%03d.js" % i)).write_text(
            "(window.TSEQ=window.TSEQ||{}),Object.assign(TSEQ,%s);" % json.dumps(shards.get(i, {}), separators=(",", ":"), sort_keys=True))
    return n, len(used)


def safe(pid):
    return re.sub(r"[^A-Za-z0-9._-]", "_", pid)


CTX = {}


# The six classes of protein, in the order of the index tiles: (tile text, filter button text).
STATUSES = {
    "none": ("with no hit from any program", "No program has a hit"),
    "one": ("with a hit from exactly one program", "Exactly one program has a hit"),
    "agree": ("where hits of two programs overlap and name the same Pfam family", "Programs agree on a Pfam family"),
    "disagree": ("where hits of two programs overlap and name different Pfam families", "Programs name different families"),
    "overlap_unknown": ("where hits of two programs overlap but no Pfam family is known for one of them",
                        "Overlap, but no family known"),
    "apart": ("where two or more programs have hits, but on different parts of the protein", "Hits do not overlap"),
}


def index_page(rows, ctx):

    total = len(rows)
    pos_best = [r["best_e"] for r in rows if r["best_e"]]
    lo_idx = max(-300, math.floor(math.log10(min(pos_best)))) if pos_best else -1
    cats = {s: 0 for s in STATUSES}
    only = {f: 0 for f in PROGRAM_ORDER}
    for r in rows:
        cats[r["status"]] += 1
        if r["n"] == 1:
            only[r["progs"][0]] += 1
    assert sum(cats.values()) == total
    out = [page_head("Protein annotation tools, side by side", 0)]
    out.append('<header class="wide"><h1>Which tools found something in each protein</h1></header><main class="wide">')
    nums = [("all", total, "proteins in all")] + [(s, cats[s], STATUSES[s][0]) for s in STATUSES]
    out.append('<section class="nums">%s</section>' % "".join(
        '<a class="num" href="#%s"><span class="big">%s</span><span>%s</span></a>' % (k, "{:,}".format(v), esc(t)) for k, v, t in nums))
    # Sidebar: for each program, how many proteins only it has a hit on; each line filters the table.
    side = []
    for f in PROGRAM_ORDER:
        mark = svg_inline(shape_svg(0, 0, 12, 10, PROGRAMS[f]), 12, 10)
        if ctx["program_run"][f]:
            side.append('<li><a data-f="only-%s" href="#only-%s">%s<span class="pn">%s</span><span class="pc">%s</span></a></li>'
                        % (f, f, mark, esc(PROGRAMS[f]["name"]), "{:,}".format(only[f])))
        else:
            side.append('<li class="off">%s<span class="pn">%s</span><span class="pc">not run</span></li>'
                        % (mark, esc(PROGRAMS[f]["name"])))
    aside = ('<aside class="side"><h2>Only one program has a hit</h2><p class="note">Proteins on which this '
             'program, and no other, reports a hit. Click one to list them.</p><ul>%s</ul></aside>' % "".join(side))
    filt = ['<button data-f="all">All</button>'] + ['<button data-f="%s">%s</button>' % (s, esc(STATUSES[s][1])) for s in STATUSES]
    out.append('<div class="withside"><section><h2>Proteins</h2><p class="note">One page per protein of your proteome. On each page, every '
               'tool\'s hits are drawn on the same ruler of residue positions. Each program counts once, even when it was '
               'run twice (MMseqs2 one round and three rounds; Foldseek on predicted structures and ProstT5 + Foldseek), '
               'so there are 8 programs for 10 rows; %d of the 8 were run on this proteome. A hit counts when its E-value is at most %s '
               '(hhsearch: also Prob at least 50); the protein pages also draw weaker hits, up to E-value %s. "Agree" means two hits from different programs overlap by '
               'at least half of the shorter one and name the same Pfam family or InterPro entry. The small '
               'picture is the "programs with a hit" line of each protein page: its height at each point is the '
               'number of programs (0 to %d) with a hit over that part of the protein, squeezed to the same 90-pixel width '
               'and the same 0-to-%d height for every protein whatever its length; hover over it for the position and count. "Best E-value, any tool" is the lowest E-value any tool printed '
               '(each tool computes its own, so they are not on one scale); 0 means the tool printed 0.</p>'
               '<div class="search"><label for="q">Search</label><input id="q" type="search" autocomplete="off" '
               'placeholder="a gene, protein or Pfam family, e.g. HLA, immune, MHC, RING"></div>'
               '<p class="note">The search looks at every hit drawn on the protein pages (E-value at most %s): the name, gene and species of '
               'each Swiss-Prot protein a tool matched, the Pfam families on the matched part of it, and the '
               'Pfam families hmmscan or InterProScan found on your protein. Every word typed must appear; the '
               'last column shows what matched and which programs\' hits carry it.</p>'
               '<div class="filters" role="group" aria-label="Filter the table">%s</div>'
               '<div class="eslider"><label for="emax">Best E-value at most <output id="ev">10</output></label>'
               '<input id="emax" type="range" min="%d" max="1" step="1" value="1"><span class="note">At 10, every '
               'protein is listed, including those with no hit; moves in steps of 10 times.</span></div>'
               '<label class="unnamed"><input id="unnamed" type="checkbox"> Only proteins with no gene name in '
               'their ID (FUN…_FUN…), %s of %s</label>'
               '<p id="count" class="note"></p>' % (ctx["n_run"], fmt_thr(COUNT_E), fmt_thr(E_MAX), ctx["n_run"], ctx["n_run"], fmt_thr(E_MAX), "".join(filt), lo_idx,
 "{:,}".format(sum(1 for r in rows if not r["gname"])), "{:,}".format(total)))
    trs = []
    for r in rows:
        trs.append('<tr data-s="%s" data-f="%s"><td><a class="mono" href="proteins/%s.html">%s</a>%s</td><td class="r" data-v="%d">%s</td>'
                   '<td class="r" data-v="%d">%d</td><td>%s</td><td class="r" data-v="%s">%s</td><td>%s</td></tr>' % (
                       r["status"], " ".join(r["progs"]), safe(r["id"]), esc(r["id"]),
                       '<div class="pname">%s</div>' % esc(r["name"]) if r["name"] else "",
                       r["length"], "{:,}".format(r["length"]), r["n"], r["n"], esc(r["agreed"] or "none"),
                       "%g" % r["best_e"] if r["best_e"] is not None else "1e999", fmt_e(r["best_e"]) if r["best_e"] is not None else "no hit", r["spark"]))
        if r["gname"]:
            gcell = '<td data-v="%d"><span class="mono">%s</span> %s</td>' % (
                len(r["gprogs"]), esc(r["gname"]),
                ("found by " + esc(", ".join(r["gprogs"]))) if r["gprogs"] else '<span class="no">no hit has this gene</span>')
        else:
            gcell = '<td data-v="-1"></td>'
        trs[-1] = trs[-1][:-5] + gcell + '<td class="m"></td></tr>'
        trs[-1] = trs[-1].replace("<tr ", '<tr data-n="%d" ' % (1 if r["gname"] else 0), 1)
    out.append('<div class="tablewrap"><table id="t"><thead><tr><th data-k="0" data-t="s">Protein</th><th data-k="1">Length (aa)</th>'
               '<th data-k="2">Programs with a hit (%s)</th><th data-k="3" data-t="s">Pfam family they agree on</th><th data-k="4">Best E-value, any tool</th>'
               '<th>Programs with a hit, along the protein</th><th data-k="6">Gene name in the ID, and the programs whose hits are on that gene</th><th class="m">What the search matched</th></tr></thead><tbody>%s</tbody></table></div></section>%s</div>' % (esc(of_programs(ctx["n_run"])), "".join(trs), aside))
    # The search table: every distinct text once, then for each protein (in table order) pairs of
    # [text number, program bit mask]. Kept in its own file so the table itself stays small.
    vocab, pos, per_row = [], {}, []
    for r in rows:
        pairs = []
        for txt, mask in r["terms"].items():
            if txt not in pos:
                pos[txt] = len(vocab)
                vocab.append(txt)
            pairs.extend((pos[txt], mask))
        per_row.append(pairs)
    (ctx["out"] / "search.js").write_text("window.SEARCH=%s;" % json.dumps(
        dict(programs=[PROGRAMS[f]["name"] for f in PROGRAM_ORDER], text=vocab, rows=per_row), separators=(",", ":")))
    out.append('<div id="tip" role="tooltip" hidden></div></main>%s<script src="search.js%s"></script><script src="assets/index.js%s"></script></body></html>' % (ctx["footer"], ASSET_VER[0], ASSET_VER[0]))
    (ctx["out"] / "index.html").write_text("".join(out))
    return cats, only


def footer(inputs, built, secs=None):
    rows = "".join('<tr><td class="mono">%s</td><td class="r">%s</td><td>%s</td></tr>' % (
        esc(n), "{:,} bytes".format(s) if s is not None else "", esc(m) if m else "not in the input folder")
        for n, s, m in inputs)
    return ('<footer><h2>Input files read</h2><p class="note">Built %s by build.py%s. No tool was run '
            'to make these pages; they only read the files below.</p><table class="small">'
            '<tr><th>File</th><th>Size</th><th>Last modified</th></tr>%s</table></footer>' % (
                esc(built), "" if secs is None else " in %.0f seconds" % secs, rows))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", type=Path)
    ap.add_argument("--out", type=Path, default=Path("site"))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--count-at", type=float, default=COUNT_E,
                    help="a hit counts when its E-value is at most this (a power of 10; default 0.001)")
    a = ap.parse_args()
    lg = math.log10(a.count_at) if a.count_at > 0 else None
    if lg is None or abs(lg - round(lg)) > 1e-9 or a.count_at > E_MAX:
        print("--count-at must be a power of 10 no larger than %g (the slider moves in steps of 10 times)" % E_MAX)
        sys.exit(2)
    globals()["COUNT_E"] = a.count_at
    t_start = time.time()
    # Millions of small hit records make Python's cycle collector rescan the heap again and
    # again while the files are read; the records hold no cycles, so it is turned off.
    gc.disable()
    D = a.inputs

    def inp(name):
        """D/name, or D/name.gz when only the gzipped copy is there."""
        return D / name if (D / name).exists() or not (D / (name + ".gz")).exists() else D / (name + ".gz")
    out = a.out
    if out.exists():
        # Only an empty folder or a site this script built (it holds summary.json) is replaced.
        if not out.is_dir() or (any(out.iterdir()) and not (out / "summary.json").exists()):
            print("not replacing %s: it is not empty and holds no summary.json, so it is not a site build.py made" % out)
            sys.exit(2)
        shutil.rmtree(out)

    (out / "proteins").mkdir(parents=True)
    # marks the folder as this script's from the start, so a build that stops halfway can be redone
    (out / "summary.json").write_text('{"finished": false}\n')
    here = Path(__file__).resolve().parent
    shutil.copytree(here / "assets", out / "assets")
    hsh = hashlib.sha1()
    for f in sorted((here / "assets").rglob("*")):
        if f.is_file():
            hsh.update(f.read_bytes())
    ASSET_VER[0] = "?v=" + hsh.hexdigest()[:10] + datetime.datetime.now().strftime("%y%m%d%H%M%S")

    def stamp(name):
        p = inp(name)
        if not p.exists():
            return (name, None, None)
        if p.is_dir():
            files = sorted(p.glob("*.hhr"))
            size = sum(f.stat().st_size for f in files)
            mt = max((f.stat().st_mtime for f in files), default=p.stat().st_mtime)
            return (name + "/ (%d .hhr files)" % len(files), size, datetime.datetime.fromtimestamp(mt).strftime("%Y-%m-%d %H:%M"))
        st = p.stat()
        return (p.name, st.st_size, datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M"))

    names_read = ["proteins.fasta", "swissprot.fasta"] + [t["file"] for t in TRACKS] + [
        "phmmer.out", "jackhmmer.out", "hmmscan.out", "protein_names.tsv", "plddt.tsv", "expression.tsv", "leaderboard.tsv", "run_settings.tsv", "kmerseek_alphabets.rs"]
    inputs = [stamp(n) for n in names_read]

    prots = {}
    for h, s in read_fasta(inp("proteins.fasta")).items():
        s = s.rstrip("*")
        prots[h.split()[0]] = dict(seq=s, length=len(s))
    # Optional: a name for some of your proteins (ID, tab, name), shown on its page and searched.
    prot_names = {}
    if inp("protein_names.tsv").exists():
        for line in opentext(inp("protein_names.tsv")):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] in prots and not line.startswith("#"):
                prot_names[f[0]] = f[1]
    entry_names, targets, target_info = {}, {}, {}
    if inp("swissprot.fasta").exists():
        for h, s in read_fasta(inp("swissprot.fasta")).items():
            acc = accession(h.split()[0])
            entry, desc, gene, org = swissprot_header(h)
            if entry:
                entry_names[acc] = entry
            if desc or gene or org:
                target_info[acc] = (desc, gene, org)
            targets[acc] = s
    print("read %d proteins, %d Swiss-Prot sequences (%.0f s)" % (len(prots), len(targets), time.time() - t_start), flush=True)

    for key in ("phmmer", "jackhmmer", "hmmscan"):
        if inp(key + ".out").exists():
            n0 = len(HMMER_ALN)
            parse_hmmer_out(inp(key + ".out"), key, set(prots))
            print("  %s.out: %d alignments (%.0f s)" % (key, len(HMMER_ALN) - n0, time.time() - t_start), flush=True)
    hits, run_status, target_ann, names = {}, {}, {}, {}
    kmer_meta = {}
    hmmscan_target = {}
    for t in TRACKS:
        p = inp(t["file"])
        if not p.exists():
            run_status[t["key"]] = "not run"
            continue
        run_status[t["key"]] = "run"
        if t["kind"] == "m8":
            got = parse_m8(p, t["key"])
        elif t["kind"] == "domtbl_seq":
            got = parse_domtbl_seq(p, t["key"])
        elif t["kind"] == "domtbl_pfam":
            got, hmmscan_target, pn = parse_domtbl_pfam(p, prots)
            names.update(pn)
        elif t["kind"] == "interproscan":
            got, ips_target, nm = parse_interproscan(p, prots)
            names.update(nm)
            target_ann.update(ips_target)
        elif t["kind"] == "kmerseek":
            got, kmer_meta = parse_kmerseek(p)
        elif t["kind"] == "hhr":
            got = {}
            for f in sorted(p.glob("*.hhr")):
                q, hs = parse_hhr(f)
                USED_TARGETS.update(h["target"] for h in hs if h["passes"])
                got.setdefault(q or f.stem, []).extend(hs)
        unknown = [q for q in got if q not in prots]
        if unknown:
            print("  %s: %d query ids are not in proteins.fasta and are skipped, e.g. %s" % (t["file"], len(unknown), unknown[0]), flush=True)
        for q, hs in got.items():
            if q in prots:
                hits.setdefault(q, []).extend(hs)
        print("  %s: %d hits on %d proteins (%.0f s)" % (t["file"], sum(len(v) for v in got.values()), len(got), time.time() - t_start), flush=True)
    # Pfam families of targets: InterProScan's output on that target if there is one, else hmmscan's
    for acc, doms in hmmscan_target.items():
        target_ann.setdefault(acc, doms)

    settings = {t["key"]: t["setting"] for t in TRACKS}
    databases = {t["key"]: t["database"] for t in TRACKS}
    if kmer_meta.get("moltype"):
        settings["kmerseek"] = "alphabet %s, k = %s" % (", ".join(sorted(kmer_meta["moltype"])), ", ".join(sorted(kmer_meta["ksize"])))
    if inp("run_settings.tsv").exists():
        head, rws = read_tsv_table(inp("run_settings.tsv"))
        for r in rws:
            d = dict(zip(head, r))
            if d.get("track") in settings:
                if d.get("setting"):
                    settings[d["track"]] = d["setting"]
                if d.get("database"):
                    databases[d["track"]] = d["database"]

    leaderboard = {}
    if inp("leaderboard.tsv").exists():
        head, rws = read_tsv_table(inp("leaderboard.tsv"))
        for r in rws:
            d = dict(zip(head, r))
            k = d.get("track") or d.get("tool")
            key = next((t["key"] for t in TRACKS if k in (t["key"], t["name"])), None)
            if key:
                leaderboard[key] = "On the leaderboard, zebrafish, Pfam truth, read %s: Fmax %s, $%s per 1,000 queries" % (
                    d.get("date_read", "?"), d.get("fmax", "?"), d.get("cost_per_1000_queries", "?"))

    plddt = {}
    if inp("plddt.tsv").exists():
        with opentext(inp("plddt.tsv")) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[0] in prots:
                    plddt[f[0]] = [float(x) for x in f[1].split(",") if x]
    expression = None
    if inp("expression.tsv").exists():
        head, rws = read_tsv_table(inp("expression.tsv"))
        expression = (head[1:], {r[0]: [fnum(x) for x in r[1:]] for r in rws})

    alphabets = {}
    for cand in (inp("kmerseek_alphabets.rs"),):
        if cand.exists():
            alphabets = parse_alphabets(cand)
    used = kmer_meta.get("moltype", set())
    alphabets = {k: v for k, v in alphabets.items() if k in used}

    program_run = {f: any(run_status[t["key"]] == "run" for t in TRACKS if t["program"] == f) for f in PROGRAM_ORDER}
    n_run = sum(program_run.values())  # M in "N of the M programs run"
    built = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    CTX.update(dict(proteins=prots, prot_names=prot_names, hits=hits, target_ann=target_ann, names=names, entry_names=entry_names,
                    target_info=target_info,
                    targets=targets, run_status=run_status, settings=settings, databases=databases,
                    leaderboard=leaderboard, plddt=plddt, expression=expression, alphabets=alphabets,
                    out=out, footer="", program_run=program_run, n_run=n_run))
    CTX["footer"] = footer(inputs, built)
    CTX["n_shards"], n_used = write_target_shards(out, hits, targets)
    print("wrote %d Swiss-Prot sequences into %d files (%.0f s)" % (n_used, CTX["n_shards"], time.time() - t_start), flush=True)
    print("parsed all inputs (%.0f s); drawing %d pages with %d workers" % (time.time() - t_start, len(prots), a.jobs), flush=True)
    ids = list(prots)
    gc.freeze()  # workers share the parsed hits without copying them
    if a.jobs > 1:
        with multiprocessing.get_context("fork").Pool(a.jobs) as pool:
            rows = pool.map(build_protein, [(i, None) for i in ids], chunksize=max(1, len(ids) // (a.jobs * 20)))
    else:
        rows = [build_protein((i, None)) for i in ids]
    secs = time.time() - t_start
    CTX["footer"] = footer(inputs, built, secs)
    cats, only = index_page(rows, CTX)
    summary = dict(count_at=COUNT_E, drawn_up_to=E_MAX, programs_run=n_run, total=len(rows),
 **cats, only_program=only,
                   run_status=run_status, seconds=round(time.time() - t_start, 1))

    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
