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

E_MAX = 1e-3  # a hit counts when its E-value is at most this
HH_PROB_MIN = 50.0  # hhsearch hits also need Prob of at least this
OVERLAP = 0.5  # two intervals overlap when they share at least half of the shorter one
TOP_N = 5  # hits drawn per track before "show all"

# One colour per tool family. None of these is a hue the pLDDT scale uses (dark blue, light
# blue, yellow, orange), and the greys are kept for the protein line, expression bars and
# InterProScan, at clearly different lightness.
FAMILIES = {
    "mmseqs2": {"name": "MMseqs2", "color": "#5E57A3"},
    "kmerseek": {"name": "kmerseek", "color": "#B2182B"},
    "foldseek": {"name": "Foldseek", "color": "#3A7D1E"},
    "phmmer": {"name": "phmmer", "color": "#C46A1B"},
    "jackhmmer": {"name": "jackhmmer", "color": "#E7298A"},
    "hhsearch": {"name": "HH-suite (hhsearch)", "color": "#00796B"},
    "hmmscan": {"name": "hmmscan with Pfam", "color": "#6B6B00"},
    "interproscan": {"name": "InterProScan", "color": "#4A4A4A"},
}
PROTEIN_LINE = "#C8C8C8"
AGREEMENT = "#000000"
EXPRESSION = "#8C8C8C"
PLDDT = [  # AlphaFold's own colours and bins
    (90, 101, "#0053D6", "very high (above 90)"),
    (70, 90, "#65CBF3", "confident (70 to 90)"),
    (50, 70, "#FFDB13", "low (50 to 70)"),
    (0, 50, "#FF7D45", "very low (below 50)"),
]

# The ten tracks, in drawing order. "second" marks the second run of a family, drawn as an
# open box. Settings are the leaderboard's setting lines; run_settings.tsv can replace them.
TRACKS = [
    dict(key="mmseqs2", file="mmseqs2.m8", family="mmseqs2", second=False, name="MMseqs2",
         setting="sensitivity 7 (-s 7), one round", database="Swiss-Prot",
         score="evalue", kind="m8"),
    dict(key="mmseqs2_3rounds", file="mmseqs2_3rounds.m8", family="mmseqs2", second=True,
         name="MMseqs2, three rounds", setting="sensitivity 7 (-s 7), three rounds",
         database="Swiss-Prot", score="evalue", kind="m8"),
    dict(key="kmerseek", file="kmerseek.tsv", family="kmerseek", second=False, name="kmerseek",
         setting="", database="Swiss-Prot", score="region_evalue", kind="kmerseek"),
    dict(key="foldseek", file="foldseek.m8", family="foldseek", second=False, name="Foldseek",
         setting="3Di shape letters and amino acids together",
         database="AlphaFold models of Swiss-Prot", score="evalue", kind="m8"),
    dict(key="prostt5", file="prostt5_foldseek.m8", family="foldseek", second=True,
         name="ProstT5 + Foldseek", setting="3Di letters predicted from sequence by ProstT5",
         database="AlphaFold models of Swiss-Prot", score="evalue", kind="m8"),
    dict(key="phmmer", file="phmmer.domtblout", family="phmmer", second=False, name="phmmer",
         setting="phmmer's default options, E-value cutoff 10", database="Swiss-Prot",
         score="i-Evalue", kind="domtbl_seq"),
    dict(key="jackhmmer", file="jackhmmer.domtblout", family="jackhmmer", second=False,
         name="jackhmmer", setting="three rounds of search (-N 3); the last round is drawn",
         database="Swiss-Prot", score="i-Evalue", kind="domtbl_seq"),
    dict(key="hhsearch", file="hhsearch", family="hhsearch", second=False,
         name="HH-suite (hhsearch)", setting="one .hhr file per protein",
         database="Swiss-Prot or PDB70 profiles", score="Prob and E-value", kind="hhr"),
    dict(key="hmmscan", file="hmmscan.domtblout", family="hmmscan", second=False,
         name="hmmscan with Pfam", setting="scans Pfam's models; no other protein searched",
         database="Pfam-A", score="i-Evalue", kind="domtbl_pfam"),
    dict(key="interproscan", file="interproscan.tsv", family="interproscan", second=False,
         name="InterProScan", setting="its member databases' own cutoffs", database="InterPro",
         score="score column", kind="interproscan"),
]
TRACK = {t["key"]: t for t in TRACKS}
FAMILY_ORDER = list(FAMILIES)

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


def m8_hit(line, key):
    cols = M8_COLS + (["alntmscore", "theader"] if key in ("foldseek", "prostt5") else ["theader"])
    d = dict(zip(cols, line.rstrip("\n").split("\t")))
    e = fnum(d["evalue"])
    sc = [("evalue", d["evalue"]), ("fident", d.get("fident", "")), ("bits", d.get("bits", ""))]
    if key in ("foldseek", "prostt5"):
        sc.append(("alntmscore", d.get("alntmscore", "")))
    sc = [(k, v if v not in ("", None) else "not in the file") for k, v in sc]
    return hit(key, d["qstart"], d["qend"], accession(d["target"]), int(d["tstart"]),
               int(d["tend"]), e, sc, e is not None and e <= E_MAX)


def parse_m8(path, key):
    out = {}
    with open(path) as fh:
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
    with open(path) as fh:
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
    return hit(key, f[15], f[16], accession(f[0]), int(f[17]), int(f[18]), e, sc, e is not None and e <= E_MAX)


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
        e = fnum(f[12])
        if f[3] in proteins:
            sc = [("i-Evalue", f[12]), ("score", f[13]), ("Pfam", pf)]
            h = hit("hmmscan", f[19], f[20], f[0], None, None, e, sc, e is not None and e <= E_MAX)
            h["names"] = {pf}
            out.setdefault(f[3], []).append(h)
        elif e is not None and e <= E_MAX:
            target_pfam.setdefault(accession(f[3]), []).append((int(f[19]), int(f[20]), {pf}))
    return out, target_pfam, pfam_names


def parse_interproscan(path, proteins):
    out, target_ann, names = {}, {}, {}
    with open(path) as fh:
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
                h = hit("interproscan", f[6], f[7], desc or sig, None, None, e, sc, passes)
                h["names"] = nm
                h["db"] = analysis
                out.setdefault(prot, []).append(h)
            elif passes and nm:
                target_ann.setdefault(accession(prot), []).append((int(f[6]), int(f[7]), nm))
    return out, target_ann, names


HHR_HIT = re.compile(r"^\s*(\d+)\s+(.*?)\s+(\d+(?:\.\d+)?)\s+(\S+)\s+(\S+)\s+(-?\d+(?:\.\d+)?)\s+"
                     r"(-?\d+(?:\.\d+)?)\s+(\d+)\s+(\d+)-(\d+)\s+(\d+)-(\d+)\s*\(\s*(\d+)\)\s*$")


def parse_hhr(path):
    """One .hhr file: the hit table (query and template ranges, Prob, E-value) and, for each
    hit, the alignment of the query and template residues as hhsearch printed it."""
    query, hits, alns, cur = None, {}, {}, None
    with open(path) as fh:
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
    with open(path) as fh:
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
    with open(path) as fh:
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


AGREE_OTHERS = [[]]  # the other agreeing families of the last protein agreement() looked at


def agreement(hits, target_ann):
    """Returns (families with a hit, status, agreed name id).
    status: 'none', 'one', 'agree', 'overlap_unknown' (overlapping hits, no family known for
    at least one of them, none shared), or 'disagree'.

    Proteins with many hits (a kinase can have thousands) make comparing every pair slow, so
    hits are grouped by the family name they carry, repeated ranges are dropped, and the search
    for an overlapping pair between two tool families stops at the first one found, best
    E-value first. This gives the same answer as comparing every pair."""
    passing = [h for h in hits if h["passes"]]
    fams = sorted({TRACK[h["track"]]["family"] for h in passing}, key=FAMILY_ORDER.index)
    AGREE_OTHERS[0] = []
    if not fams:
        return fams, "none", None
    if len(fams) == 1:
        return fams, "one", None
    by_name, unknown_hits, all_by_fam = {}, {}, {}
    for h in passing:
        fam = TRACK[h["track"]]["family"]
        e = h["evalue"] if h["evalue"] is not None else math.inf
        iv = (h["q0"], h["q1"])
        all_by_fam.setdefault(fam, {})
        all_by_fam[fam][iv] = min(e, all_by_fam[fam].get(iv, math.inf))
        if h["_names"] is None:
            unknown_hits.setdefault(fam, {})[iv] = min(e, unknown_hits.get(fam, {}).get(iv, math.inf))
            continue
        for n in h["_names"]:
            d = by_name.setdefault(n, {}).setdefault(fam, {})
            d[iv] = min(e, d.get(iv, math.inf))

    def best_pair(xa, xb):
        """Lowest E-value of an overlapping pair, one range from each dict, or None."""
        cand = sorted([(e, iv, xb) for iv, e in xa.items()] + [(e, iv, xa) for iv, e in xb.items()])
        for e, (a0, a1), other in cand:
            if any(overlaps(a0, a1, b0, b1) for (b0, b1) in other):
                return e
        return None

    support = {}
    for n, per_fam in by_name.items():
        fs = list(per_fam)
        for i, fa in enumerate(fs):
            for fb in fs[i + 1:]:
                e = best_pair(per_fam[fa], per_fam[fb])
                if e is not None:
                    got, best = support.get(n, (set(), math.inf))
                    support[n] = (got | {fa, fb}, min(best, e))
    if support:
        # The family named by the most tool families. On a tie, the one whose own domain on
        # this protein (an hmmscan or InterProScan hit) has the best E-value, then the best
        # E-value of an agreeing pair.
        own = {}
        for h in passing:
            if h["track"] in ("hmmscan", "interproscan") and h["_names"]:
                for n in h["_names"]:
                    own[n] = min(own.get(n, math.inf), h["evalue"] if h["evalue"] is not None else math.inf)
        order = sorted(support, key=lambda n: (-len(support[n][0]), own.get(n, math.inf), support[n][1], n))
        AGREE_OTHERS[0] = order[1:]
        return fams, "agree", order[0]
    for fa, xa in unknown_hits.items():
        for fb, xb in all_by_fam.items():
            if fb != fa and best_pair(xa, xb) is not None:
                return fams, "overlap_unknown", None
    return fams, "disagree", None


def coverage(hits, length):
    """Number of families with a passing hit over each residue (list of length `length`)."""
    per_fam = {}
    for h in hits:
        if h["passes"]:
            per_fam.setdefault(TRACK[h["track"]]["family"], []).append((h["q0"], h["q1"]))
    diff = [0] * (length + 2)
    for ivs in per_fam.values():
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
    return len(s) * px * 0.56


def lanes_for(hits):
    """Stack boxes that overlap into separate lanes, best E-value first."""
    lanes, out = [], []
    for h in hits:
        for li, ends in enumerate(lanes):
            if all(h["q0"] > b or h["q1"] < a for a, b in ends):
                ends.append((h["q0"], h["q1"]))
                out.append(li)
                break
        else:
            lanes.append([(h["q0"], h["q1"])])
            out.append(len(lanes) - 1)
    return out, max(1, len(lanes))


LANE_H, BOX_H, ROW_PAD = 18, 14, 8


def box_svg(h, idx, x0, x1, y, track):
    fam = FAMILIES[track["family"]]
    col = fam["color"]
    w = max(2.0, x1 - x0)
    if track["second"]:
        rect = ('<rect x="%.1f" y="%d" width="%.1f" height="%d" rx="2" fill="#fff" stroke="%s" '
                'stroke-width="2"/>' % (x0 + 1, y + 1, max(1.0, w - 2), BOX_H - 2, col))
        tcol = "#1a1a1a"
    else:
        rect = '<rect x="%.1f" y="%d" width="%.1f" height="%d" rx="2" fill="%s"/>' % (x0, y, w, BOX_H, col)
        tcol = "#ffffff"
    label = h["label"]
    lab = ""
    if text_w(label) + 6 <= w:
        shown = label
    else:
        n = int((w - 6) / (11 * 0.56)) - 1
        shown = label[:n] + "…" if n >= 3 else ""
    if shown:
        lab = '<text x="%.1f" y="%d" fill="%s" class="bl">%s</text>' % (x0 + 3, y + BOX_H - 3, tcol, esc(shown))
    return ('<g class="hit" tabindex="0" data-i="%d">%s%s</g>' % (idx, rect, lab))


def page_svg(p, hits_by_track, cov, plddt, W, px, run_status):
    """The track area as two SVGs (names on the left, the drawing on the right), with the best
    TOP_N hits of each tool. Each row is a <g class="row"> in its own coordinates, so the
    page's script can redraw a tool's row with all of its hits ("Show all")."""
    L = p["length"]
    rows = []  # (key, height, left svg, right svg)
    r = ['<line x1="0" y1="18" x2="%.1f" y2="18" stroke="%s" stroke-width="2"/>' % (W, PROTEIN_LINE)]
    for pos in range(50, L + 1, 50):
        x = (pos - 0.5) * px
        r.append('<line x1="%.1f" y1="12" x2="%.1f" y2="24" stroke="%s"/>' % (x, x, PROTEIN_LINE))
        if pos % 100 == 0:
            r.append('<text x="%.1f" y="10" class="rl" text-anchor="middle">%d</text>' % (x, pos))
    r.append('<text x="1" y="10" class="rl">1</text>')
    r.append('<text x="%.1f" y="30" class="rl" text-anchor="end">%d</text>' % (W - 1, L))
    rows.append(("ruler", 32, '<text x="4" y="20" class="tn">residue</text>', "".join(r)))
    for t in TRACKS:
        k = t["key"]
        hs = hits_by_track.get(k, [])
        shown = hs[:TOP_N]
        lane_of, n_l = lanes_for(shown)
        h_row = ROW_PAD + n_l * LANE_H + 4
        if run_status[k] == "not run":
            g = ['<text x="4" y="%d" class="nr">not run: no %s in the input folder</text>'
                 % (ROW_PAD + 11, esc(t["file"]))]
        elif not hs:
            g = ['<text x="4" y="%d" class="nh">no hit at the threshold</text>' % (ROW_PAD + 11)]
        else:
            g = ['<line x1="0" y1="%d" x2="%.1f" y2="%d" stroke="%s" stroke-width="1.5"/>'
                 % (ROW_PAD + BOX_H // 2, W, ROW_PAD + BOX_H // 2, PROTEIN_LINE)]
        for h, li in zip(shown, lane_of):
            g.append(box_svg(h, h["i"], (h["q0"] - 1) * px, h["q1"] * px, ROW_PAD + li * LANE_H, t))
        more = len(hs) - len(shown)
        if more > 0:
            h_row = max(h_row, ROW_PAD + 34)  # room for the "+N more" line under the name
        extra = ('<text x="26" y="%d" class="tm more">+%d more</text>' % (ROW_PAD + 26, more)) if more > 0 else ""
        left = ('<g transform="translate(4,%d)">%s</g><text x="26" y="%d" class="tn">%s</text>%s'
                % (ROW_PAD, legend_mark(t), ROW_PAD + 11, esc(t["name"]), extra))
        rows.append((k, h_row, left, "".join(g)))
    if plddt is not None:
        g = []
        for i, v in enumerate(plddt[:L]):
            c = next((col for lo, hi, col, _ in PLDDT if lo <= v < hi), PLDDT[-1][2])
            g.append('<rect x="%.2f" y="%d" width="%.2f" height="%d" fill="%s"/>' % (i * px, ROW_PAD, px + 0.05, BOX_H, c))
        rows.append(("plddt", ROW_PAD + BOX_H + 8, '<text x="4" y="%d" class="tn">pLDDT</text>' % (ROW_PAD + 11), "".join(g)))
    # families with a hit: step line, 0 to 8
    ah, top = 8 * 6, 10
    pts, prev = [], None
    for i, c in enumerate(cov):
        if c != prev:
            if prev is not None:
                pts.append("L%.2f,%.1f" % (i * px, top + ah - prev * 6))
            pts.append(("M" if prev is None else "L") + "%.2f,%.1f" % (i * px, top + ah - c * 6))
            prev = c
    if prev is not None:
        pts.append("L%.2f,%.1f" % (L * px, top + ah - prev * 6))
    grid = "".join('<line x1="0" y1="%.1f" x2="%.1f" y2="%.1f" class="grid"/>' % (top + ah - v * 6, W, top + ah - v * 6) for v in (0, 2, 4, 6, 8))
    right = grid + '<path d="%s" fill="none" stroke="%s" stroke-width="1.6"/>' % ("".join(pts), AGREEMENT)
    left = '<text x="4" y="%d" class="tn">families with a hit</text>' % (top + 12) + "".join(
        '<text x="100%%" dx="-4" y="%.1f" class="rl" text-anchor="end">%d</text>' % (top + ah - v * 6 + 3, v) for v in (0, 4, 8))
    rows.append(("agree", top + ah + 8, left, right))
    y, L_out, R_out = 0, [], []
    for k, hgt, left, right in rows:
        sep_l = '<line x1="0" y1="0" x2="100%" y2="0" class="sep"/>' if y else ""
        sep_r = '<line x1="0" y1="0" x2="%.1f" y2="0" class="sep"/>' % W if y else ""
        L_out.append('<g class="row" data-k="%s" data-h="%d" transform="translate(0,%d)">%s%s</g>' % (k, hgt, y, sep_l, left))
        R_out.append('<g class="row" data-k="%s" data-h="%d" transform="translate(0,%d)">%s%s</g>' % (k, hgt, y, sep_r, right))
        y += hgt
    return y, "".join(L_out), "".join(R_out)


def legend_mark(t, w=16):
    col = FAMILIES[t["family"]]["color"]
    if t["second"]:
        return '<rect x="1" y="1" width="%d" height="%d" rx="2" fill="#fff" stroke="%s" stroke-width="2"/>' % (w - 2, BOX_H - 2, col)
    return '<rect x="0" y="0" width="%d" height="%d" rx="2" fill="%s"/>' % (w, BOX_H, col)


def svg_inline(inner, w, h, cls=""):
    vb = "" if isinstance(w, str) else ' viewBox="0 0 %s %d"' % (w, h)
    return '<svg class="%s" width="%s" height="%d"%s xmlns="http://www.w3.org/2000/svg">%s</svg>' % (
        cls, w, h, vb, inner)


def sparkline(cov, length, w=60, h=16):
    if not cov:
        return svg_inline("", w, h)
    bins = []
    for b in range(w):
        a0 = int(b * length / w)
        a1 = max(a0 + 1, int((b + 1) * length / w))
        bins.append(max(cov[a0:a1]))
    d = []
    for b, c in enumerate(bins):
        yy = h - 1 - c * (h - 2) / 8
        d.append(("M" if b == 0 else "L") + "%d,%.1f L%d,%.1f" % (b, yy, b + 1, yy))
    base = '<line x1="0" y1="%d" x2="%d" y2="%d" stroke="%s"/>' % (h - 0.5, w, h - 0.5, PROTEIN_LINE)
    return svg_inline(base + '<path d="%s" fill="none" stroke="%s" stroke-width="1.2"/>' % (" ".join(d), AGREEMENT), w, h, "spark")


# --------------------------------------------------------------------------- pages ----

CSS_LINKS = (
    '<link rel="preconnect" href="https://fonts.googleapis.com">'
    '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
    '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600'
    '&family=Source+Sans+3:wght@400;600;700&display=swap">'
)


def page_head(title, depth):
    pre = "../" * depth
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            '<title>%s</title>%s<link rel="stylesheet" href="%sassets/site.css"></head><body>'
            % (esc(title), CSS_LINKS, pre))


def legend_html(ctx, has_plddt, has_expr):
    rows = []
    rows.append('<tr><th>Mark</th><th>Tool</th><th>Setting</th><th>Searched</th><th>Score (in the box\'s tooltip)</th><th>Counts as a hit when</th></tr>')
    for t in TRACKS:
        k = t["key"]
        fam = FAMILIES[t["family"]]
        mark = svg_inline(legend_mark(t), 16, BOX_H)
        run = ctx["run_status"][k]
        name = esc(t["name"])
        if t["second"]:
            name += '<div class="sub">open box: second run in the %s family</div>' % esc(fam["name"])
        elif t["family"] in ("mmseqs2", "foldseek"):
            name += '<div class="sub">solid box: first run in the %s family</div>' % esc(fam["name"])
        coords = {"domtbl_seq": "box: hmm from-to, the positions on your protein",
                  "domtbl_pfam": "box: env from-to (the domain's envelope)",
                  "m8": "box: qstart-qend", "kmerseek": "box: region_start-region_end, plus 1 at the start",
                  "hhr": "box: Query HMM range", "interproscan": "box: start-stop"}[t["kind"]]
        name += '<div class="sub">%s</div>' % esc(coords)
        lb = ctx["leaderboard"].get(k)
        if lb:
            name += '<div class="sub">%s</div>' % esc(lb)
        if run == "not run":
            thr = "not run: no %s in the input folder" % t["file"]
            rows.append('<tr class="notrun"><td>%s</td><td>%s</td><td colspan="4">not run: no %s in the input folder</td></tr>'
                        % (mark, name, esc(t["file"])))
            continue
        elif k == "hhsearch":
            thr = "Prob ≥ 50 and E-value ≤ 0.001"
        elif k == "interproscan":
            thr = ("E-value ≤ 0.001 where the score column is an E-value; any match otherwise. "
                   "Coils, disorder, membrane and signal-peptide predictions do not count")
        else:
            thr = "%s ≤ 0.001" % t["score"]
        rows.append('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (
            mark, name, esc(ctx["settings"][k]), esc(ctx["databases"][k]), esc(t["score"]), esc(thr)))
    if has_plddt:
        sw = "".join('<span class="sw" style="background:%s"></span>%s ' % (c, esc(lbl)) for _, _, c, lbl in PLDDT)
        rows.append('<tr><td colspan="2">pLDDT</td><td colspan="4">AlphaFold\'s confidence in each residue of your predicted structure, 0 to 100: %s</td></tr>' % sw)
    line = svg_inline('<line x1="0" y1="7" x2="16" y2="7" stroke="%s" stroke-width="2"/>' % PROTEIN_LINE, 16, 14)
    rows.append('<tr><td>%s</td><td colspan="5">The protein, from residue 1 to its last residue; the ruler has a tick every 50 residues and a number every 100.</td></tr>' % line)
    step = svg_inline('<path d="M0,11 L6,11 L6,4 L12,4 L12,8 L16,8" fill="none" stroke="%s" stroke-width="1.6"/>' % AGREEMENT, 16, 14)
    rows.append('<tr><td>%s</td><td colspan="5">Families with a hit: computed from the tracks above. At each residue, the number of the 8 tool families that have a hit over it at their threshold (0 to 8). Faint grey lines mark 0, 2, 4, 6 and 8; the numbers 0, 4 and 8 are on the left.</td></tr>' % step)
    if has_expr:
        bar = svg_inline('<rect x="0" y="3" width="16" height="8" fill="%s"/>' % EXPRESSION, 16, 14)
        rows.append('<tr><td>%s</td><td colspan="5">Expression of this protein\'s gene in each sample, from expression.tsv.</td></tr>' % bar)
    return '<table class="legend">%s</table>' % "".join(rows)


def header_sentence(n, status, agreed, names):
    if n == 0:
        return "No tool reports a hit."
    s = "%d of 8 tool families report%s a hit." % (n, "s" if n == 1 else "")
    if n == 1:
        return s
    if status == "agree":
        return s + " They agree on %s." % family_label(agreed, names)
    if status == "overlap_unknown":
        return s + " Their hits overlap, but no Pfam family is known for the Swiss-Prot proteins they match."
    return s + " They do not agree on a family."


def family_label(n, names):
    nm = names.get(n)
    return "%s (%s)" % (nm, n) if nm and nm != n else n


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
    fams, status, agreed = agreement(hits, ctx["target_ann"])
    cov = coverage(hits, L)
    best_e = min([h["evalue"] for h in hits if h["passes"] and h["evalue"] is not None], default=None)

    # draw only hits at the threshold; per track, best E-value first
    by_track = {}
    data = []
    score_names = {}
    for t in TRACKS:
        hs = sorted([h for h in hits if h["track"] == t["key"] and h["passes"]],
                    key=lambda h: (h["evalue"] if h["evalue"] is not None else 1e300, h["q0"]))
        for rank, h in enumerate(hs):
            h["i"] = len(data)
            name = ctx["entry_names"].get(h["target"], h["target"]) if t["key"] not in ("hmmscan", "interproscan") else h["target"]
            h["label"] = name
            d = dict(k=t["key"], l=name, q=[h["q0"], h["q1"]], s=[v for _, v in h["scores"]])
            score_names.setdefault(t["key"], [n for n, _ in h["scores"]])
            if t["key"] not in ("interproscan", "hmmscan"):
                d["t"] = h["target"]
                d["r"] = [h["t0"], h["t1"]]
                tseq = ctx["targets"].get(h["target"])
                if t["key"] == "hhsearch":
                    d["a"] = h["aln"]
                elif t["key"] == "kmerseek":
                    d["m"] = h.get("moltype")
                elif tseq is not None and rank < TOP_N:
                    lo, hi = min(h["t0"], h["t1"]), max(h["t0"], h["t1"])
                    d["a"] = recompute(seq[h["q0"] - 1:h["q1"]], h["q0"], tseq[lo - 1:hi], lo)
            data.append(d)
        by_track[t["key"]] = hs

    plddt = ctx["plddt"].get(pid)
    W = int(max(640, min(L * 2.0, 3200)))
    px = W / L
    hgt, left, right = page_svg(p, by_track, cov, plddt, W, px, ctx["run_status"])
    blocks = ['<div class="tracks" id="tracks"><div class="names">%s</div><div class="scroll">%s</div></div>'
              % (svg_inline(left, "100%", hgt), svg_inline(right, W, hgt, "trk"))]
    n_more = sum(max(0, len(v) - TOP_N) for v in by_track.values())

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

    sentence = header_sentence(len(fams), status, agreed, ctx["names"])
    others = AGREE_OTHERS[0]
    other_note = ""
    if others:
        shown_o = [family_label(n, ctx["names"]) for n in others[:8]]
        other_note = ('<p class="note">%d other Pfam famil%s also meet%s the rule: %s%s. The one named has hits '
                      'from the most tool families; on a tie, the best E-value of its own hmmscan or InterProScan '
                      'domain on this protein.</p>' % (len(others), "y" if len(others) == 1 else "ies",
                                                       "s" if len(others) == 1 else "", esc(", ".join(shown_o)),
                                                       " and %d more" % (len(others) - 8) if len(others) > 8 else ""))
    out = [page_head("%s · tool comparison" % pid, 1)]
    out.append('<header><a class="back" href="../index.html">All proteins</a><h1 class="mono">%s</h1>'
               '<p class="len">%d amino acids</p><p class="summary">%s</p>%s</header><main>' % (esc(pid), L, esc(sentence), other_note))
    out.append('<section><h2>Legend</h2>%s</section>' % legend_html(ctx, plddt is not None, bool(expr)))
    toggle = ""
    if n_more:
        toggle = ('<button id="toggle" type="button">Show all %d hits (%d more)</button>'
                  % (sum(len(v) for v in by_track.values()), n_more))
    out.append('<section><h2>What each tool found</h2><p class="note">Each box is one hit at the '
               'threshold: one row of the tool\'s output table, drawn over the residues of this protein it covers. '
               'A Swiss-Prot protein matched in two places gives two boxes. Up to %d hits per tool are drawn, best '
               'E-value first (ties: the leftmost first); "+N more" under a tool\'s name counts the rest. Hover or '
               'tap a box for its numbers; click it to see the residues. Wide proteins scroll sideways inside the '
               'box below, so some boxes may start to the right of what you see.</p>%s%s</section>'
               % (TOP_N, toggle, "".join(blocks)))
    out.append(expr)
    out.append('<div id="tip" role="tooltip" hidden></div><aside id="panel" hidden><button id="close" type="button">Close</button><div id="pbody"></div></aside>')
    out.append('</main>%s' % ctx["footer"])
    page_data = dict(
        q=seq, px=px, w=W, top=TOP_N, hits=data, names=score_names, alpha=ctx["alphabets"],
        tracks={t["key"]: dict(name=t["name"], color=FAMILIES[t["family"]]["color"], open=t["second"]) for t in TRACKS},
        shards=ctx["n_shards"])
    out.append('<script>const PAGE=%s;</script><script src="../assets/protein.js"></script></body></html>'
               % json.dumps(page_data, separators=(",", ":")).replace("</", "<\\/"))
    (ctx["out"] / "proteins" / (safe(pid) + ".html")).write_text("".join(out))
    return dict(id=pid, length=L, n=len(fams), fams=fams, status=status,
                agreed=family_label(agreed, ctx["names"]) if agreed else None,
                best_e=best_e, spark=sparkline(cov, L))


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


def index_page(rows, ctx):
    total = len(rows)
    cats = {"none": 0, "one": 0, "agree": 0, "disagree": 0}
    only = {f: 0 for f in FAMILY_ORDER}
    for r in rows:
        c = r["status"] if r["status"] in ("none", "one", "agree") else "disagree"
        cats[c] += 1
        if r["n"] == 1:
            only[r["fams"][0]] += 1
    assert sum(cats.values()) == total
    out = [page_head("Protein annotation tools, side by side", 0)]
    out.append('<header><h1>Which tools found something in each protein</h1></header><main>')
    nums = [("all", total, "proteins in all"), ("none", cats["none"], "with no hit from any family"),
            ("one", cats["one"], "with a hit from exactly one family"),
            ("agree", cats["agree"], "where two or more families agree on a family"),
            ("disagree", cats["disagree"], "where two or more families have hits but do not agree")]
    out.append('<section class="nums">%s</section>' % "".join(
        '<a class="num" href="#%s"><span class="big">%s</span><span>%s</span></a>' % (k, "{:,}".format(v), esc(t)) for k, v, t in nums))
    fam_rows = "".join('<tr><td>%s %s</td><td class="r">%s</td><td>%s</td></tr>' % (
        svg_inline('<rect width="14" height="12" rx="2" fill="%s"/>' % FAMILIES[f]["color"], 14, 12),
        esc(FAMILIES[f]["name"]), "{:,}".format(only[f]) if ctx["family_run"][f] else "", "" if ctx["family_run"][f] else "not run")
        for f in FAMILY_ORDER)
    out.append('<section><h2>Proteins only one family annotates</h2><table class="small"><tr><th>Tool family</th><th>Proteins where it is the only one with a hit</th><th></th></tr>%s</table></section>' % fam_rows)
    filt = ['<button data-f="all">All</button>', '<button data-f="one">Exactly one family has a hit</button>',
            '<button data-f="none">No family has a hit</button>', '<button data-f="disagree">Families disagree</button>',
            '<button data-f="agree">Families agree</button>']
    for f in FAMILY_ORDER:
        if ctx["family_run"][f]:
            filt.append('<button data-f="only-%s">Only %s has a hit</button>' % (f, esc(FAMILIES[f]["name"])))
    out.append('<section><h2>Proteins</h2><p class="note">One page per protein of your proteome. On each page, every '
               'tool\'s hits are drawn on the same ruler of residue positions. A tool family is one program; its two '
               'runs (MMseqs2 one round and three rounds; Foldseek on predicted structures and ProstT5 + Foldseek) '
               'count once. A hit counts when its E-value is at most 0.001 '
               '(hhsearch: also Prob at least 50). "Agree" means two hits from different families overlap by '
               'at least half of the shorter one and name the same Pfam family or InterPro entry. The small '
               'picture is the "families with a hit" line of each protein page: its height at each point is the '
               'number of families (0 to 8) with a hit over that part of the protein, squeezed to the same 60-pixel width '
               'for every protein whatever its length. "Best E-value, any tool" is the lowest E-value any tool printed '
               '(each tool computes its own, so they are not on one scale); 0 means the tool printed 0.</p>'
               '<div class="filters" role="group" aria-label="Filter the table">%s</div><p id="count" class="note"></p>' % "".join(filt))
    trs = []
    for r in rows:
        trs.append('<tr data-s="%s" data-f="%s"><td class="mono"><a href="proteins/%s.html">%s</a></td><td class="r" data-v="%d">%s</td>'
                   '<td class="r" data-v="%d">%d</td><td>%s</td><td class="r" data-v="%s">%s</td><td>%s</td></tr>' % (
                       "disagree" if r["status"] == "overlap_unknown" else r["status"], " ".join(r["fams"]), safe(r["id"]), esc(r["id"]),
                       r["length"], "{:,}".format(r["length"]), r["n"], r["n"], esc(r["agreed"] or "none"),
                       "%g" % r["best_e"] if r["best_e"] is not None else "1e999", fmt_e(r["best_e"]) if r["best_e"] is not None else "no hit", r["spark"]))
    out.append('<div class="tablewrap"><table id="t"><thead><tr><th data-k="0" data-t="s">Protein</th><th data-k="1">Length (aa)</th>'
               '<th data-k="2">Families with a hit (of 8)</th><th data-k="3" data-t="s">Agreed family</th><th data-k="4">Best E-value, any tool</th>'
               '<th>Families with a hit, along the protein</th></tr></thead><tbody>%s</tbody></table></div></section>' % "".join(trs))
    out.append('</main>%s<script src="assets/index.js"></script></body></html>' % ctx["footer"])
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
    a = ap.parse_args()
    t_start = time.time()
    # Millions of small hit records make Python's cycle collector rescan the heap again and
    # again while the files are read; the records hold no cycles, so it is turned off.
    gc.disable()
    D = a.inputs
    out = a.out
    if out.exists():
        shutil.rmtree(out)
    (out / "proteins").mkdir(parents=True)
    here = Path(__file__).resolve().parent
    shutil.copytree(here / "assets", out / "assets")

    def stamp(name):
        p = D / name
        if not p.exists():
            return (name, None, None)
        if p.is_dir():
            files = sorted(p.glob("*.hhr"))
            size = sum(f.stat().st_size for f in files)
            mt = max((f.stat().st_mtime for f in files), default=p.stat().st_mtime)
            return (name + "/ (%d .hhr files)" % len(files), size, datetime.datetime.fromtimestamp(mt).strftime("%Y-%m-%d %H:%M"))
        st = p.stat()
        return (name, st.st_size, datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M"))

    names_read = ["proteins.fasta", "swissprot.fasta"] + [t["file"] for t in TRACKS] + [
        "plddt.tsv", "expression.tsv", "leaderboard.tsv", "run_settings.tsv", "kmerseek_alphabets.rs"]
    inputs = [stamp(n) for n in names_read]

    prots = {}
    for h, s in read_fasta(D / "proteins.fasta").items():
        s = s.rstrip("*")
        prots[h.split()[0]] = dict(seq=s, length=len(s))
    entry_names, targets = {}, {}
    if (D / "swissprot.fasta").exists():
        for h, s in read_fasta(D / "swissprot.fasta").items():
            first = h.split()[0]
            acc = accession(first)
            parts = first.split("|")
            if len(parts) > 2 and parts[2]:
                entry_names[acc] = parts[2]
            targets[acc] = s
    print("read %d proteins, %d Swiss-Prot sequences (%.0f s)" % (len(prots), len(targets), time.time() - t_start), flush=True)

    hits, run_status, target_ann, names = {}, {}, {}, {}
    kmer_meta = {}
    hmmscan_target = {}
    for t in TRACKS:
        p = D / t["file"]
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
    # target families: InterProScan's output on that target if there is one, else hmmscan's
    for acc, doms in hmmscan_target.items():
        target_ann.setdefault(acc, doms)

    settings = {t["key"]: t["setting"] for t in TRACKS}
    databases = {t["key"]: t["database"] for t in TRACKS}
    if kmer_meta.get("moltype"):
        settings["kmerseek"] = "alphabet %s, k = %s" % (", ".join(sorted(kmer_meta["moltype"])), ", ".join(sorted(kmer_meta["ksize"])))
    if (D / "run_settings.tsv").exists():
        head, rws = read_tsv_table(D / "run_settings.tsv")
        for r in rws:
            d = dict(zip(head, r))
            if d.get("track") in settings:
                if d.get("setting"):
                    settings[d["track"]] = d["setting"]
                if d.get("database"):
                    databases[d["track"]] = d["database"]

    leaderboard = {}
    if (D / "leaderboard.tsv").exists():
        head, rws = read_tsv_table(D / "leaderboard.tsv")
        for r in rws:
            d = dict(zip(head, r))
            k = d.get("track") or d.get("tool")
            key = next((t["key"] for t in TRACKS if k in (t["key"], t["name"])), None)
            if key:
                leaderboard[key] = "On the leaderboard, zebrafish, Pfam truth, read %s: Fmax %s, $%s per 1,000 queries" % (
                    d.get("date_read", "?"), d.get("fmax", "?"), d.get("cost_per_1000_queries", "?"))

    plddt = {}
    if (D / "plddt.tsv").exists():
        with open(D / "plddt.tsv") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[0] in prots:
                    plddt[f[0]] = [float(x) for x in f[1].split(",") if x]
    expression = None
    if (D / "expression.tsv").exists():
        head, rws = read_tsv_table(D / "expression.tsv")
        expression = (head[1:], {r[0]: [fnum(x) for x in r[1:]] for r in rws})

    alphabets = {}
    for cand in (D / "kmerseek_alphabets.rs",):
        if cand.exists():
            alphabets = parse_alphabets(cand)
    used = kmer_meta.get("moltype", set())
    alphabets = {k: v for k, v in alphabets.items() if k in used}

    family_run = {f: any(run_status[t["key"]] == "run" for t in TRACKS if t["family"] == f) for f in FAMILY_ORDER}
    built = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    CTX.update(dict(proteins=prots, hits=hits, target_ann=target_ann, names=names, entry_names=entry_names,
                    targets=targets, run_status=run_status, settings=settings, databases=databases,
                    leaderboard=leaderboard, plddt=plddt, expression=expression, alphabets=alphabets,
                    out=out, footer="", family_run=family_run))
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
    summary = dict(total=len(rows), **cats, only_family=only, run_status=run_status, seconds=round(time.time() - t_start, 1))
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
