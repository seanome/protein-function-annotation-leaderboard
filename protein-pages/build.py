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
COUNT_E = 10.0  # a hit counts (sentence as built, index, summary.json) when its E-value is at most this; --count-at
SHOW_E = 10.0  # the E-value slider on a protein page starts here, so every drawn hit shows at first; --show-at
TARGET_PFAM_E = 1e-3  # a Swiss-Prot protein's Pfam families: its hmmscan or InterProScan domains at or below this
TARGET_PFAM_FILE = "target_pfam_interpro.tsv"  # Pfam domains of the other targets, from InterPro's REST API
HH_PROB_MIN = 50.0  # hhsearch hits also need Prob of at least this
OVERLAP = 0.5  # two intervals overlap when they share at least half of the shorter one
TOP_N = 5  # hits drawn per track before "show all"

# Colour is the kind of program (three hues and a dark grey); inside a kind, the box's shape
# tells the programs apart. Checked in greyscale and with a deuteranopia simulation (Machado
# et al. 2009, severity 1, applied in linear RGB). In normal vision the four differ by at least
# 7 in lightness (CIE L*). Under deuteranopia they differ by at least 24 in colour (CIE76), but
# sequence search and k-mer are then only 1 apart in lightness, so in greyscale those two are
# told apart by hue alone, and by the row each box sits in. None is a pLDDT colour, and none is
# the orange or blue of the residue panel's hydrophobic and polar classes.
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
# Programs drawn only on protein.html, from the whole-proteome tables: never TRACKS, which build.py parses
WHOLE_TRACKS = {"diamond": {"name": "DIAMOND", "cat": "seq", "shape": "half", "color": CATEGORIES["seq"]["color"]},
                "unifire": {"name": "UniFIRE", "cat": "pfam", "shape": "hatched", "color": CATEGORIES["pfam"]["color"]}}
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
    # max_aa: Biomni predicted structures only for proteins of at most 200 aa, so on a longer
    # protein this row says "not run" whatever the input folder holds
    dict(key="foldseek", file="foldseek.m8", program="foldseek", second=False, name="Foldseek, 3-D structures",
         setting="Foldseek on 3-D structures predicted by ESMFold2 (fast model), for proteins of 200 aa or less",
         database="AlphaFold models of Swiss-Prot", score="evalue", kind="m8", max_aa=200),
    dict(key="prostt5", file="prostt5_foldseek.m8", program="foldseek", second=True,
         name="ProstT5 + Foldseek", setting="3Di letters predicted from sequence by ProstT5",
         database="AlphaFold models of Swiss-Prot", score="evalue", kind="m8"),
    dict(key="phmmer", file="phmmer.domtblout", program="phmmer", second=False, name="phmmer",
         setting="phmmer's default options, E-value cutoff 10", database="Swiss-Prot",
         score="i-Evalue", kind="domtbl_seq"),
    dict(key="jackhmmer", file="jackhmmer.domtblout", program="jackhmmer", second=False,
         name="jackhmmer", setting="three rounds of search (-N 3); the last round is shown",
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
# The InterProScan feature predictions drawn on its row, as dashed boxes that never count:
# analysis -> {signature or "*": the word drawn}. Phobius's topology pieces (cytoplasmic,
# non-cytoplasmic, the n-, h- and c-regions of a signal peptide) are left out.
IPS_DRAWN = {"MobiDBLite": {"*": "disorder"}, "Coils": {"*": "coiled coil"},
             "Phobius": {"TRANSMEMBRANE": "membrane span", "SIGNAL_PEPTIDE": "signal peptide"},
             "TMHMM": {"*": "membrane span"}, "DeepTMHMM": {"*": "membrane span"},
             "SignalP_EUK": {"*": "signal peptide"}, "SignalP-noTM": {"*": "signal peptide"},
             "SignalP-TM": {"*": "signal peptide"}, "SignalP": {"*": "signal peptide"}}


def feature_word(analysis, sig):
    """The word drawn for an InterProScan feature prediction, or None if it is not drawn."""
    m = IPS_DRAWN.get(analysis, {})
    return m.get(sig) or m.get("*")


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


AF_NAME = re.compile(r"AF-([A-Z0-9]+(?:-\d+)?)-F\d+")


def af_accession(t):
    """UniProt accession from an AlphaFold model name, keeping the isoform number:
    AF-P12345-F1-model_v6 gives P12345, AF-C3VPR6-2-F1-model_v6 gives C3VPR6-2."""
    m = AF_NAME.match(t)
    return m.group(1) if m else t.split("-")[1]


def accession(t):
    """Swiss-Prot accession from a target name: sp|P12345|NAME, AF-P12345-F1-model_v4, P12345."""
    if t.startswith("AF-"):
        return af_accession(t)
    if t.startswith(("sp|", "tr|")):
        parts = t.split("|")
        return parts[1] if len(parts) > 1 else t
    return t


def counts(h):
    """A drawn hit counts when its E-value is at most COUNT_E. A hit with no E-value (an
    InterProScan match from a member database whose score is not an E-value) counts if drawn.
    An InterProScan feature prediction (disorder, a coil, a membrane span) is drawn, never counted."""
    return h["passes"] and not h.get("feature") and (h["evalue"] is None or h["evalue"] <= COUNT_E)


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
    sequences, from --format-output ...,qaln,taln), give the tool's own alignment. Foldseek on
    structures also has lddt after alntmscore."""
    extra = {"foldseek": ["alntmscore", "lddt"], "prostt5": ["alntmscore"]}.get(key, [])
    cols = M8_COLS + extra + ["theader", "qaln", "taln"]
    d = dict(zip(cols, line.rstrip("\n").split("\t")))
    e = fnum(d["evalue"])
    sc = [("evalue", d["evalue"]), ("fident", d.get("fident", "")), ("bits", d.get("bits", ""))]
    sc += [(k, d.get(k, "")) for k in extra]
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
            word = feature_word(analysis, sig)
            if word:
                passes = True
            elif analysis in IPS_FEATURES:
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
                sc = [("member database", analysis), ("signature", sig), ("its score or E-value", score),
                      ("InterPro entry", (ipr + " " + ipr_desc) if ipr else "none")]
                # the score is kept as an E-value only where it is one
                h = hit("interproscan", f[6], f[7], desc or sig, None, None,
                        e if analysis in IPS_EVALUE else None, sc, passes)
                h["names"] = nm
                h["db"] = analysis
                if word:
                    h.update(feature=True, target=word, names=set())
                out.setdefault(prot, []).append(h)
            elif passes and not word and nm and (e is None or e <= TARGET_PFAM_E):
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
                                t0=int(m.group(11)), t1=int(m.group(12)), tlen=int(m.group(13)))
            continue
        if line.startswith("No "):
            cur = int(line.split()[1])
            alns[cur] = {"q": [], "t": [], "qs": None, "ts": None, "desc": ""}
        elif cur is not None and line.startswith(">"):
            # the template's full description: a PDB70 line is ">6XYZ_A Protein name; ..."
            parts = line[1:].split(None, 1)
            alns[cur]["desc"] = parts[1].strip() if len(parts) > 1 else ""
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
        # Not every template is in swissprot.fasta (PDB70 chains never are): keep its length and
        # description from the .hhr so the page can still name it and place the match on it
        out[-1]["tlen"] = h["tlen"]
        out[-1]["desc"] = (a or {}).get("desc", "")
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
            h["ksize"] = int(d["ksize"]) if (d.get("ksize") or "").isdigit() else None
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
    # dayhoff6's table is in sourmash, not in alphabets.rs: sourmash src/core/src/encodings.rs,
    # Dayhoff's six groups C, AGPST, DENQ, HKR, ILMV, FWY, written as a-f.
    out.setdefault("dayhoff6", {r: c for c, g in zip("abcdef", ["C", "AGPST", "DENQ", "HKR", "ILMV", "FWY"]) for r in g})
    return out


def seed_runs(qseg, tseg, table, k):
    """Where kmerseek's search could have started on an ungapped hit: every window of k columns
    in which the two segments are in the same class of the hit's alphabet at all k positions.
    Returns (runs, n): runs are the 0-based, end-inclusive column ranges covered by such windows
    (overlapping windows merged), and n is the number of windows. A residue outside the
    alphabet (X, U) matches nothing."""
    n_col = min(len(qseg), len(tseg))
    cq = [table.get(c.upper()) for c in qseg[:n_col]]
    ct = [table.get(c.upper()) for c in tseg[:n_col]]
    same = [a is not None and a == b for a, b in zip(cq, ct)]
    runs, n, i = [], 0, 0
    while i < n_col:
        if not same[i]:
            i += 1
            continue
        j = i
        while j + 1 < n_col and same[j + 1]:
            j += 1
        if j - i + 1 >= k:
            runs.append((i, j))
            n += j - i + 2 - k
        i = j + 1
    return runs, n


def read_target_pfam(path):
    """target_pfam_interpro.tsv (made by test_inputs/fetch_target_pfam.py): Pfam domains of
    Swiss-Prot targets from the InterPro REST API. Returns ({accession: [(start, end, {Pfam id}),
    ...]}, {Pfam id: name}). An accession listed with Pfam id '-' gets an empty list: InterPro
    was asked and lists no Pfam match, which is not the same as not knowing."""
    doms, names = {}, {}
    with opentext(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    col = {c: i for i, c in enumerate(rows[0])}
    for f in rows[1:]:
        acc, pf = f[col["accession"]], f[col["pfam"]]
        lst = doms.setdefault(acc, [])
        if pf == "-":
            continue
        lst.append((int(f[col["start"]]), int(f[col["end"]]), {pf}))
        if f[col["pfam_name"]]:
            names.setdefault(pf, f[col["pfam_name"]])
    return doms, names


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
# For an 'agree' protein: one dict per family that meets the rule, the named one first:
# acc, programs (keys), n (programs), hits (distinct hits carrying it), a, b (residues the hits
# span on your protein), e (best E-value of an agreeing pair), lists (program -> its hits, as a set).
AGREE_INFO = [[]]
AGREE_DECIDER = [None]  # what picked the named family over the next: "programs", "hits", "evalue" or None (no rival)
PROG_TRACKS = [{}]  # program -> the names of its rows with a counted hit ("ProstT5 + Foldseek")


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
    AGREE_INFO[0] = []
    AGREE_DECIDER[0] = None
    PROG_TRACKS[0] = {}
    for t in TRACKS:
        if any(h["track"] == t["key"] for h in passing):
            PROG_TRACKS[0].setdefault(t["program"], []).append(t["name"])
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
        # The Pfam family carried by the most programs; on a tie, the one with the most hits (a hit
        # = a Swiss-Prot protein and a range on your protein, counted once even when two programs
        # report it); then the best E-value of an agreeing pair; then the accession.
        info = {}
        for n, (ps, e) in support.items():
            hs = [h for h in passing if h["_names"] and n in h["_names"] and TRACK[h["track"]]["program"] in ps]
            lists = {}
            for h in hs:
                lists.setdefault(TRACK[h["track"]]["program"], set()).add((h.get("target"), h["q0"], h["q1"]))
            info[n] = dict(acc=n, programs=sorted(ps, key=PROGRAM_ORDER.index), n=len(ps),
                           hits=len(hs),  # every program's hits, as each row counts them
                           a=min(h["q0"] for h in hs), b=max(h["q1"] for h in hs), e=e, lists=lists)
        order = sorted(support, key=lambda n: (-info[n]["n"], -info[n]["hits"], info[n]["e"], n))
        if len(order) > 1:
            a, b = info[order[0]], info[order[1]]
            AGREE_DECIDER[0] = "programs" if a["n"] != b["n"] else "hits" if a["hits"] != b["hits"] else "evalue" if a["e"] != b["e"] else "accession"
        AGREE_OTHERS[0] = order[1:]
        AGREE_INFO[0] = [info[n] for n in order]
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


# A Pfam (PF00069, PF00069.32) or InterPro (IPR000719) accession inside a piece of text.
ACC_RE = re.compile(r"\b(PF\d{5}|IPR\d{6})(\.\d+)?\b")


def acc_url(acc):
    """InterPro's page for a Pfam or InterPro accession (without a version number)."""
    return "https://www.ebi.ac.uk/interpro/entry/%s/%s/" % ("pfam" if acc.startswith("PF") else "InterPro", acc)


def acc_html(s):
    """s escaped for HTML, with every Pfam and InterPro accession in it made a link. The
    accessions are letters and digits only, so they read the same before and after escaping."""
    return ACC_RE.sub(lambda m: '<a class="mono" href="%s" target="_blank" rel="noopener">%s</a>' % (
        acc_url(m.group(1)), m.group(0)), esc(s))


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
    if shape == "half":  # an outline with its lower half filled
        return ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="2" fill="#fff" stroke="%s"%s/><rect x="%.1f" y="%.1f" width="%.1f" '
                'height="%.1f" fill="%s"%s/>' % (x + 0.5, y + 0.5, max(1.0, w - 1), h - 1, col, extra, x, y + h / 2, w, h / 2, col, extra))
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


def feature_svg(x, y, w, h, col):
    """An InterProScan feature prediction: a dashed box in InterProScan's colour, no fill (transparent, so the whole box takes a click, not only its dashed edge)."""
    return ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="2" fill="transparent" stroke="%s" stroke-width="1.5" '
            'stroke-dasharray="3 2"/>' % (x + 1, y + 1, max(1.0, w - 2), h - 2, col))


def box_svg(h, idx, x0, x1, y, track, side="in", lx=None):
    prog = PROGRAMS[track["program"]]
    w = max(2.0, x1 - x0)
    op = e_opacity(h["evalue"])
    if h.get("feature"):
        op, rect = 1.0, feature_svg(x0, y, w, BOX_H, prog["color"])
    else:
        rect = shape_svg(x0, y, w, BOX_H, prog, ' opacity="%.2f"' % op if op < 1 else "")
    ty = y + BOX_H - 3
    if side == "in":
        solid = prog["shape"] in ("filled", "rounded") and op >= 0.6 and not h.get("feature")  # white text needs a dark box
        tcol, cls = ("#ffffff", "bl") if solid else ("#1a1a1a", "bl halo")
        lab = '<text x="%.1f" y="%d" fill="%s" class="%s">%s</text>' % (x0 + 3, ty, tcol, cls, esc(h["label"]))
    else:  # beside the box, with a leader line in the box's colour
        lw = text_w(h["label"])
        l0, l1 = (x0 + w, lx - 1) if side == "right" else (lx + lw + 1, x0)
        lead = ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>'
                % (l0, y + BOX_H / 2, l1, y + BOX_H / 2, prog["color"])) if l1 > l0 else ""
        lab = lead + '<text x="%.1f" y="%d" fill="#1a1a1a" class="bl out">%s</text>' % (lx, ty, esc(h["label"]))
    return '<g class="hit" tabindex="0" data-i="%d">%s%s</g>' % (idx, rect, lab)


# Programs from the Biomni whole-proteome run that give no position on the protein
# (test_inputs/make_biomni_best.py writes biomni_best.json). Each gets a row in the drawing that
# names its best hit, and counts in the sentence at the top; their cutoffs are their own, so the
# E-value slider does not change them. (key, name, how it was run, what counts as a hit); the name
# and how-line are short to fit the names column.
WHOLE = [
    ("diamond", "DIAMOND", "against Swiss-Prot", "E-value 0.001 or below"),
    ("jackhmmer", "jackhmmer, immune", "from 24 immune seeds", "E-value 0.001 or below"),
    ("fantasia", "FANTASIA", "protein language model", "distance below 0.15"),
    ("unifire", "UniFIRE", "UniProt's rules", "any prediction"),
]
# FANTASIA always names its most similar protein. Below this distance its pick shares a GO term
# with DIAMOND's best hit for about 70% of proteins with both, above it for 30-40%
# (figures/fantasia_vs_diamond.png, full_inputs/plot_fantasia_vs_diamond.py); only then it counts.
FANTASIA_MAX_DIST = 0.15


def svg_trunc(s, width, px=12):
    """s cut to fit width pixels of 12 px text (0.55 of the font size per character), with "…"."""
    n = int(width / (px * 0.55))
    return s if len(s) <= n else s[:max(1, n - 1)] + "…"


def whole_rows(b, W):
    """The rows of the programs in WHOLE for one protein: [(key, name, how, found, run, right svg)]."""
    if not b:
        return []
    out = []
    uni = lambda acc: ('<a href="https://www.uniprot.org/uniprotkb/%s/entry" target="_blank" rel="noopener">'
                       '<tspan class="lk">%s</tspan></a>' % (esc(acc), esc(acc)))
    def hit_text(acc, rest, y):
        return '<text x="4" y="%d" class="wt">%s %s</text>' % (y, uni(acc), esc(svg_trunc(rest, W - 12 * 0.55 * (len(acc) + 2))))
    for key, name, how, cutoff in WHOLE:
        v = b.get(key)
        y1, y2 = ROW_PAD + 11, ROW_PAD + 27
        if v is None:
            continue
        found, run, g = bool(v), v != "not run", []
        if v == "not run":
            found = False
            g.append('<text x="4" y="%d" class="nr">%s</text>' % (y1, esc(
                "not run: FANTASIA was run only on proteins with no Swiss-Prot match or in the immune-gene table"
                if key == "fantasia" else "not run: proteins over 5,000 amino acids were skipped")))
        elif not v:
            g.append('<text x="4" y="%d" class="wn">no hit (%s)</text>' % (y1, cutoff))
        elif key == "diamond":
            acc, nm, gene, org, pid_, qcov, e = v
            g.append(hit_text(acc, "%s%s, %s" % (nm, " (%s)" % gene if gene else "", org), y1))
            g.append('<text x="4" y="%d" class="wn">E-value %s%s · whole protein, no position</text>' % (
                y2, esc(e), "" if pid_ is None else ", %.0f%% identity over %.0f%% of your protein" % (pid_, qcov or 0)))
        elif key == "jackhmmer":
            fam, seed, rnd, e = v
            g.append('<text x="4" y="%d" class="wt">%s</text>' % (y1, esc(svg_trunc("found from the seed %s (%s)" % (seed, fam), W))))
            g.append('<text x="4" y="%d" class="wn">E-value %s, round %s · whole protein, no position</text>' % (y2, esc(e), rnd))
        elif key == "fantasia":
            acc, nm, gene, org, dist = v
            found = dist < FANTASIA_MAX_DIST
            g.append(hit_text(acc, "%s%s, %s" % (nm, " (%s)" % gene if gene else "", org), y1))
            g.append('<text x="4" y="%d" class="wn">%s · whole protein, no position</text>' % (y2, esc(
                "distance %.2f, close enough to count (below %.2f)" % (dist, FANTASIA_MAX_DIST) if found else
                "distance %.2f: not close enough to count (below %.2f counts); FANTASIA always names one protein" % (dist, FANTASIA_MAX_DIST))))
        elif key == "unifire":
            what, n = v
            g.append('<text x="4" y="%d" class="wt">%s</text>' % (y1, esc(svg_trunc(what, W))))
            g.append('<text x="4" y="%d" class="wn">%d prediction%s, listed under the rows · whole protein, no position</text>' % (y2, n, "" if n == 1 else "s"))
        out.append((key, name, how, found, run, "".join(g)))
    return out


def read_structure(path):
    """A predicted structure (PDB format): its text, the per-residue confidence from the B-factor
    column of each CA atom on a 0 to 100 scale (a file with values of at most 1, as ESMFold writes,
    is multiplied by 100), and the median distance between neighbouring CA atoms, which is about
    3.8 angstrom in any real protein chain. A first `REMARK   1` line names the predictor; it is
    shown on the page."""
    text = Path(path).read_text()
    made = next((l[11:].strip() for l in text.splitlines() if l.startswith("REMARK   1 ")), "")
    ca = {}
    for line in text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            ca[int(line[22:26])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[60:66]))
    res = sorted(ca)
    bs = [ca[r][3] for r in res]
    scale = 100.0 if bs and max(bs) <= 1.0 else 1.0
    gaps = sorted(math.dist(ca[a][:3], ca[b][:3]) for a, b in zip(res, res[1:]) if b == a + 1)
    return dict(text=text, file=Path(path).name, made=made, plddt=[b * scale for b in bs], scale=scale, n=len(res),
                ca_ca=gaps[len(gaps) // 2] if gaps else None)


HIGHLIGHT = "#E6007E"  # residues pointed at or dragged across: drawing, pLDDT strip and structure


def structure_html(pid, st):
    """The 3-D structure section of a drawn page (assets/structure.js draws it in Mol*, PDBe's build, with Mol*'s own controls): a key,
    the pLDDT strip along the protein right above the structure, the structure, then the text.
    Pointing at a hit, or dragging across the strip, colours those residues on the structure."""
    mean = sum(st["plddt"]) / max(1, len(st["plddt"]))
    n = len(st["plddt"])
    warn = ""
    if st["ca_ca"] is not None and not 3.6 <= st["ca_ca"] <= 4.0:
        warn = ('<p class="warn">These coordinates are not a physical protein chain: neighbouring residues\' '
                'central carbon atoms (CA) are a median %.1f &#8491; apart, where every real protein has 3.8 &#8491;. '
                'The ribbons follow the file, so their shape cannot be trusted.</p>' % st["ca_ca"])
    made = (' The file says how it was made: <span class="mono">%s</span>. pTM, where given, is the predictor\'s '
            'confidence in the overall fold, 0 to 1; below 0.5 the arrangement of the parts is a guess even where each '
            'part is confident.' % esc(st["made"])) if st.get("made") else ""
    sw = "".join('<span class="sw" style="background:%s"></span>%s ' % (c, esc(lbl)) for _, _, c, lbl in PLDDT)
    key = ('<p class="note skey"><b>pLDDT</b>, the predictor\'s confidence in each residue: %s'
           '<span class="sw" style="background:%s"></span>residues you point at or select. The strip runs from residue 1 '
           '(left) to %d (right); drag across it, or click a box in the rows, to select residues. '
           '<b id="molsel" aria-live="polite"></b></p>' % (sw, HIGHLIGHT, n))
    rects = []
    for i, v in enumerate(st["plddt"]):
        c = next((col for lo, hi, col, _ in PLDDT if lo <= v < hi), PLDDT[-1][2])
        rects.append('<rect x="%d" y="0" width="1.05" height="10" fill="%s"/>' % (i, c))
    strip = ('<svg id="pstrip" viewBox="0 0 %d 10" preserveAspectRatio="none" role="img" '
             'aria-label="pLDDT of each residue, residue 1 to %d">%s<rect id="pmask1" x="0" y="0" width="0" height="10" '
             'fill="#fff" opacity="0.7"/><rect id="pmask2" x="0" y="0" width="0" height="10" fill="#fff" opacity="0.7"/>'
             '<rect id="pbrush" x="0" y="0" width="0" height="10" fill="none" stroke="%s" stroke-width="2.5" '
             'vector-effect="non-scaling-stroke"/></svg>' % (n, n, "".join(rects), HIGHLIGHT))
    return ('<section id="structure"><h2>3-D structure (predicted)</h2>%s%s<div id="mol" data-scale="%g" data-n="%d" '
            'data-hl="%s"></div>'
            '<p class="note">Ribbons made from <span class="mono">%s</span> (%d residues), coloured by pLDDT as in the '
            'strip above them (mean %.0f here). Point at a box in the rows or a row of the table, or drag across the '
            'strip, to colour those residues %s on the structure and the strip; the rest turns grey. Clicking a box or '
            'a row, or dragging across the strip, keeps the selection until you click the strip again. Point at the structure to see a residue\'s number and mark it on the strip. '
            'The viewer is Mol* (PDBe\'s build, the one the PDB and AlphaFold DB use): drag to turn, scroll to zoom, '
            'and the sliders button on its right edge opens the sequence and the Structure Tools panels (display style, colours, selections, measurements).%s</p>%s</section>'
            % (key, strip, st["scale"], n, HIGHLIGHT, esc(st["file"]), st["n"], mean, "magenta", made, warn))


def page_svg(p, hits_by_track, cov, plddt, W, px, run_status, m=8, whole=()):
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
        real = [h for h in hs if not h.get("feature")]
        shown = real[:TOP_N] + [h for h in hs if h.get("feature")]
        place, n_l = place_labels(shown, px, W)
        h_row = ROW_PAD + n_l * LANE_H + 4
        if run_status[k] == "too long":
            g = ['<text x="4" y="%d" class="nr">not run: structures were predicted only for proteins of %d aa or less</text>'
                 % (ROW_PAD + 11, t["max_aa"])]
        elif run_status[k] == "not run" and k == "foldseek" and run_status.get("prostt5") != "not run":
            # Readers looked here for Foldseek and missed the next row
            g = ['<text x="4" y="%d" class="nr">not run on 3-D structures. Foldseek on shape letters '
                 'predicted from the sequence is the next row, ProstT5 + Foldseek.</text>' % (ROW_PAD + 11)]
        elif run_status[k] == "not run":
            g = ['<text x="4" y="%d" class="nr">not run: no %s in the input folder</text>'
                 % (ROW_PAD + 11, esc(t["file"]))]
        elif not hs:
            g = ['<text x="4" y="%d" class="nh">no hit with E-value at most %s</text>' % (ROW_PAD + 11, fmt_thr(COUNT_E))]
        elif not real:  # only feature predictions: the row still says it has no hit
            g = ['<line x1="0" y1="%d" x2="%.1f" y2="%d" stroke="%s" stroke-width="1.5"/>'
                 % (ROW_PAD + BOX_H // 2, W, ROW_PAD + BOX_H // 2, PROTEIN_LINE)]
            n_l += 1
            h_row = ROW_PAD + n_l * LANE_H + 4
            g.append('<text x="4" y="%d" class="nh">no hit with E-value at most %s; dashed: predictions, not counted</text>'
                     % (ROW_PAD + (n_l - 1) * LANE_H + 11, fmt_thr(COUNT_E)))
        else:
            g = ['<line x1="0" y1="%d" x2="%.1f" y2="%d" stroke="%s" stroke-width="1.5"/>'
                 % (ROW_PAD + BOX_H // 2, W, ROW_PAD + BOX_H // 2, PROTEIN_LINE)]
        for h, (li, side, lx) in zip(shown, place):
            g.append(box_svg(h, h["i"], (h["q0"] - 1) * px, h["q1"] * px, ROW_PAD + li * LANE_H, t, side, lx))
        more = len(real) - min(len(real), TOP_N)
        if more > 0:
            h_row = max(h_row, ROW_PAD + 34)  # room for the "+N more" line under the name
        extra = ('<text x="26" y="%d" class="tm more">+%d more hit%s</text>' % (ROW_PAD + 26, more, "" if more == 1 else "s")) if more > 0 else ""
        left = ('<g transform="translate(4,%d)">%s</g><text x="26" y="%d" class="tn">%s</text>%s'
                % (ROW_PAD, legend_mark(t), ROW_PAD + 11, esc(t["name"]), extra))
        rows.append((k, h_row, left, "".join(g)))
    for key, name, how, found, run, right in whole:
        rows.append(("w_" + key, ROW_PAD + 34, '<text x="26" y="%d" class="tn">%s</text><text x="26" y="%d" class="tm">%s</text>'
                     % (ROW_PAD + 11, esc(name), ROW_PAD + 27, esc(how)), right))
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
    rows.insert(1, ("agree", top + ah + 8, left, right))  # under the ruler, above the programs' rows
    y, L_out, R_out = 0, [], []
    for k, hgt, left, right in rows:
        sep_l = '<line x1="0" y1="0" x2="100%" y2="0" class="sep"/>' if y else ""
        sep_r = '<line x1="0" y1="0" x2="%.1f" y2="0" class="sep"/>' % W if y else ""
        L_out.append('<g class="row" data-k="%s" data-h="%d" transform="translate(0,%d)">%s%s</g>' % (k, hgt, y, sep_l, left))
        R_out.append('<g class="row" data-k="%s" data-h="%d" transform="translate(0,%d)">%s%s</g>' % (k, hgt, y, sep_r, right))
        y += hgt
    return y, "".join(L_out), "".join(R_out)


def legend_mark(t, w=16, h=BOX_H):
    return shape_svg(0, 0, w, h, PROGRAMS.get(t["program"]) or WHOLE_TRACKS[t["program"]])


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
    # media="print" keeps this request from holding up the first paint; the page shows in system fonts until it arrives
    '<link rel="stylesheet" media="print" onload="this.media=\'all\'" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600'
    '&family=Source+Sans+3:wght@400;600;700&display=swap">'
)


ASSET_VER = [""]  # "?v=<hash of assets/>": a rebuilt site makes browsers fetch the new scripts


def page_head(title, depth):
    pre = "../" * depth
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            '<title>%s</title>%s<link rel="stylesheet" href="%sassets/site.css%s"></head><body>%s'
            % (esc(title), CSS_LINKS, pre, ASSET_VER[0], hatch_defs()))


# Which Swiss-Prot each search used: not the same release for every program, so each drawn page
# shows this table under the drawing. Releases and protein counts as test_inputs/PROVENANCE.md
# says how each was checked (2025_01 by comparing accessions, 2026-11-03).
SWISSPROT_USED = [
    (("kmerseek",), "Swiss-Prot release 2026_03", "all of it", "575,748"),
    (("mmseqs2_3rounds", "phmmer", "jackhmmer"), "Swiss-Prot release 2025_01 (5 February 2025)",
     "without its 120 sea-squirt (Ascidiacea) proteins", "572,850"),
    (("mmseqs2",), "Swiss-Prot release 2026_03 (Biomni's run)", "all of it", "575,748"),
    (("diamond",), "Swiss-Prot release 2026_03 (Biomni's run)", "all of it", "575,748"),
    (("foldseek", "prostt5"), "AlphaFold's models of Swiss-Prot proteins, Foldseek's prebuilt database "
     "(not a sequence release; it keeps models of retired accessions)", "sea-squirt proteins kept", "590,183"),
]


def swissprot_html(ctx, whole):
    """The table of which Swiss-Prot each program run on this page searched."""
    names = {t["key"]: t["name"] for t in TRACKS if ctx["run_status"][t["key"]] not in ("not run", "too long")}
    names.update({w[0]: w[1] for w in whole if w[4] and w[0] == "diamond"})
    rows = []
    for keys, what, part, n in SWISSPROT_USED:
        progs = [names[k] for k in keys if k in names]
        if progs:
            rows.append("<tr><td>%s</td><td>%s</td><td>%s</td><td class=\"num\">%s</td></tr>"
                        % (esc(", ".join(progs)), esc(what), esc(part or "all of it"), esc(n)))
    if not rows:
        return ""
    return ('<section id="swissprot"><h2>Which Swiss-Prot each program searched</h2><p class="note">The programs did '
            'not all search the same copy of Swiss-Prot, so a protein can be missing from one program\'s hits only '
            'because it was not in that copy. Protein names, genes and species in the boxes and the table come from '
            'Swiss-Prot 2026_03, or from UniProt\'s website on 2026-11-01 for proteins not in it, whatever copy was searched.</p><div class="tablewrap"><table><thead><tr>'
            '<th>Program</th><th>Database searched</th><th>Part searched</th><th>Proteins</th></tr></thead><tbody>%s'
            '</tbody></table></div></section>' % "".join(rows))


def legend_html(ctx, has_plddt, has_expr, has_struct=False):
    """The legend, in a narrow column left of the drawing (above it on a phone), so the reader
    meets the key before the marks: one line per tool, its mark then its name. The rule for a
    hit, the setting, the database and the columns that give the box open on a click."""
    items = []
    for t in TRACKS:
        k = t["key"]
        prog = PROGRAMS[t["program"]]
        mark = svg_inline(legend_mark(t, 12, 10), 12, 10)
        if ctx["run_status"][k] in ("not run", "too long"):
            why = " on proteins over %d aa" % t["max_aa"] if ctx["run_status"][k] == "too long" else ""
            items.append('<li class="off">%s<span>%s: not run%s</span></li>' % (mark, esc(t["name"]), why))
            continue
        c, d = fmt_thr(COUNT_E), fmt_thr(E_MAX)
        if k == "hhsearch":
            thr = "counts at Prob ≥ 50 and E-value at most %s; shown up to E-value %s" % (c, d)
        elif k == "interproscan":
            thr = ("counts at E-value at most %s where there is one, shown up to %s; disorder, coiled-coil, "
                   "membrane-span and signal-peptide predictions are shown as dashed boxes and do not count" % (c, d))
        else:
            thr = "counts at %s at most %s; shown up to %s" % (t["score"], c, d)
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
    if ctx["run_status"].get("interproscan") != "not run":
        dash = svg_inline(feature_svg(0, 0, 12, 10, PROGRAMS["interproscan"]["color"]), 12, 10)
        items.append('<li>%s<span>on the InterProScan row: a disorder, coiled-coil, membrane-span or signal-peptide '
                     'prediction; shown, never counted as a hit</span></li>' % dash)
    if has_plddt:
        sw = "".join('<span class="sw" style="background:%s"></span>%s ' % (c, esc(lbl)) for _, _, c, lbl in PLDDT)
        items.append('<li><span></span><span><b>pLDDT</b>, the structure predictor\'s confidence in each residue, 0 to 100: %s</span></li>' % sw)
    if has_struct:
        hl = svg_inline('<rect x="0.5" y="1.5" width="11" height="7" fill="%s"/>' % HIGHLIGHT, 12, 10)
        items.append('<li>%s<span>on the pLDDT strip and the 3-D structure: the residues of the hit you point at or click, '
                     'or the ones you drag across on the strip</span></li>' % hl)
    step = svg_inline('<path d="M0,9 L4,9 L4,2 L8,2 L8,6 L12,6" fill="none" stroke="%s" stroke-width="1.4"/>' % AGREEMENT, 12, 10)
    items.append('<li>%s<span>programs with a hit at each residue, 0 to %d (the programs run), each program counted once</span></li>'
                 % (step, ctx["n_run"]))
    if has_expr:
        bar = svg_inline('<rect x="0.5" y="2.5" width="11" height="5" fill="none" stroke="%s"/>' % EXPRESSION, 12, 10)
        items.append('<li>%s<span>expression of this protein\'s gene in each sample, from expression.tsv</span></li>' % bar)
    return ('<aside class="side legend"><h2>Legend</h2>%s<p class="note">Click a tool for its setting, database '
            'and what makes a hit. The programs searched different copies of Swiss-Prot: '
            '<a href="#swissprot">which one each searched</a>.</p><ul>%s</ul></aside>' % (kind_key(), "".join(items)))


def of_programs(m, full=8):
    """'of 8 programs' when every program was run, else 'of the 5 programs run'. full: how many
    programs there are (8 with positions, plus the WHOLE programs when biomni_best.json is read)."""
    return "of %d programs" % full if m == full else "of the %d programs run" % m


def header_sentence(n, status, agreed, names, m=8, full=8):
    """The sentence at the top of a protein page. n: programs with a hit, counting the WHOLE
    programs; status and agreed: from the hits along the protein only. m: programs with at least
    one input file."""
    if n == 0:
        return "No program reports a hit." if m == full else "None of the %d programs run reports a hit." % m
    s = "%d %s report%s a hit." % (n, of_programs(m, full), "s" if n == 1 else "")
    if status in ("none", "one"):
        return s
    if status == "agree":
        return s + " Their hits along the protein agree on %s." % family_label(agreed, names)
    if status == "overlap_unknown":
        return s + (" Their hits overlap, but for one of them no Pfam family is known on the matched part of the "
                    "Swiss-Prot protein.")
    if status == "disagree":
        fams = DISAGREE_NAMES[0]
        # "Pkinase PF00069 (phmmer)": name and accession without a second pair of brackets
        listed = ", ".join("%s (%s)" % (" ".join(dict.fromkeys([names.get(n) or n, n])),
                                        ", ".join(prog_label(p) for p in ps))

                           for n, ps in fams[:6])
        more = " and %d more" % (len(fams) - 6) if len(fams) > 6 else ""
        return s + " Their hits overlap but name different Pfam families: %s%s." % (listed, more)
    return s + (" Their hits are on different parts of the protein: no two programs' hits share at least half "
                "of the shorter one.")



def prog_label(p):
    """A program as its rows with a counted hit are named on the page ("ProstT5 + Foldseek", not
    "Foldseek"); the program's own name when none is known."""
    names = PROG_TRACKS[0].get(p) or [PROGRAMS[p]["name"]]
    if len(names) > 1:
        return "%s (both runs)" % PROGRAMS[p]["name"]
    return names[0].replace(", ", " (", 1) + ")" if ", " in names[0] else names[0]


def and_list(xs):
    xs = list(xs)
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1] if xs else ""


def agree_note(names, progs):
    """The note under an 'agree' sentence, from AGREE_INFO: what agree means, which programs carry
    the named family and whether their hit lists are the same, the programs that name other
    families, the other families that meet the rule, and which rule picked the named one."""
    info = AGREE_INFO[0]
    if not info:
        return ""
    top = info[0]
    ps = top["programs"]
    s = ('"Agree" means: the matched proteins of hits from at least two programs carry the same Pfam family, and '
         'the hits overlap on your protein. Here that is %s' % and_list([prog_label(p) for p in ps]))
    lists = [top["lists"].get(p, set()) for p in ps]
    if len(ps) == 2 and lists[0] == lists[1]:
        s += (", and the %d %s hit%s on proteins carrying %s %s the same as %s's, so one hit list stands behind the headline"
              % (len(lists[1]), prog_label(ps[1]), "" if len(lists[1]) == 1 else "s", names.get(top["acc"]) or top["acc"],
                 "is" if len(lists[1]) == 1 else "are", prog_label(ps[0])))
    elif len(ps) > 2 and all(x == lists[0] for x in lists[1:]):
        s += ", and their hit lists are the same, so one hit list stands behind the headline"
    rest = [p for p in progs if p not in ps]
    s += ("; %s name%s other families." % (and_list([prog_label(p) for p in rest]), "s" if len(rest) == 1 else "")) if rest else "."
    others = info[1:]
    if others:
        same = all(o["programs"] == ps for o in others)
        s += " %d other famil%s also meet%s the rule%s: %s%s." % (
            len(others), "y" if len(others) == 1 else "ies", "s" if len(others) == 1 else "",
            ", from the same %s" % ("two programs" if len(ps) == 2 else "programs") if same else "",
            ", ".join(family_label(o["acc"], names) for o in others[:8]),
            " and %d more" % (len(others) - 8) if len(others) > 8 else "")
        nm = names.get(top["acc"]) or top["acc"]
        d = AGREE_DECIDER[0]
        rivals = [o for o in others if o["n"] == top["n"]] if d != "programs" else others[:1]
        if d == "programs":
            why = "it is carried by the most programs (%d, against %d)" % (top["n"], others[0]["n"])
        elif d == "hits":
            why = "it has the most hits (%d, against %s for %s, also carried by %d programs)" % (
                top["hits"], and_list([str(o["hits"]) for o in rivals[:3]]),
                and_list([names.get(o["acc"]) or o["acc"] for o in rivals[:3]]) + (" and %d more" % (len(rivals) - 3) if len(rivals) > 3 else ""), top["n"])
        elif d == "evalue":
            why = "the programs and hits tie, and it has the best E-value of an agreeing pair (%s)" % fmt_e(top["e"])
        else:
            why = "the programs, hits and E-values tie, so the first accession is named"
        s += " %s is named because %s. The rule: most programs, then most hits, then best E-value." % (nm, why)
    return '<p class="note">%s</p>' % acc_html(s)


def sentence_steps(hits, ctx, lo, n_whole=0, m_whole=None):
    """The sentence at the top and its note at each E-value slider step from 10**lo to 10, so the
    page can change them with the slider. Both are HTML, with Pfam and InterPro accessions as
    links. The hits that count change only at a step equal to some hit's E-value rounded up to a
    power of 10, so only those steps are worked out; the page uses the last entry at or below the
    slider's step. Returns [[step, sentence, note], ...]."""
    global COUNT_E
    keep = COUNT_E
    cuts = {lo} | {min(1, max(lo, math.ceil(math.log10(h["evalue"]) - 1e-9)))
                   for h in hits if h["passes"] and h["evalue"]}
    out, prev = [], None
    try:
        for pos in sorted(cuts):
            COUNT_E = 10.0 ** pos
            progs, status, agreed = agreement(hits, ctx["target_ann"])
            cur = (acc_html(header_sentence(len(progs) + n_whole, status, agreed, ctx["names"], ctx["n_run"] + (ctx.get("n_whole", 0) if m_whole is None else m_whole),
                                            8 + ctx.get("n_whole", 0))),
                   agree_note(ctx["names"], progs) if status == "agree" else "")
            if cur != prev:
                out.append([pos, cur[0], cur[1]])
                prev = cur
    finally:
        COUNT_E = keep
    return out


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


# ------------------------------------------------------------------- dense page ----
# Every drawn protein page in the design of the BHF page ("BHF: what every annotation program
# found", claude.ai artifact JgFXM7JTxrdpgiw9B5dEXF, 2026-11-03; the dense-results-page skill):
# a title block that scrolls away, a sticky header (E-value ladder, search, buttons, three-line
# key, column names), one grid of program rows that unfold into hit rows, the 3-D model sticky at
# the right, then the Swiss-Prot table and the input files. assets/dense-bhf.js draws the grid
# from PAGE and the globals WHOLE, NOTRUN and HOW written here.

WHOLE_DENSE = {  # key -> (name on the row, how it was run)
    "diamond": ("DIAMOND against Swiss-Prot", "run by Biomni; E-value 0.001 or below"),
    "jackhmmer": ("jackhmmer with immune-gene seeds", "run by Biomni in the reverse direction: immune-gene seed proteins "
                  "(24 in Biomni's run; 20 have a hit in its table) as queries against the Botryllus proteome; E-value 0.001 or below"),
    "fantasia": ("FANTASIA (language model)", "run by Biomni; always names its most similar protein"),
    "unifire": ("UniFIRE (UniProt's rules)", "run by Biomni on InterProScan's matches"),
}


def uni_link(acc):
    return ('<a class="mono" href="https://www.uniprot.org/uniprotkb/%s/entry" target="_blank" rel="noopener">%s</a>'
            % (esc(acc), esc(acc)))


def dense_whole(b):
    """window.WHOLE: [key, name, how, text, flag, second line] per whole-protein program run on this
    protein; flag 1 = a hit that counts, 2 = named but not counted (FANTASIA above its cutoff), 0 = nothing."""
    out = []
    for key, name, how, cutoff in WHOLE:
        v = (b or {}).get(key)
        if v is None:
            continue
        nm, hw = WHOLE_DENSE[key]
        if v == "not run":
            continue  # not on the no-hit line; twelve_text names the programs run
        elif not v:
            out.append(["w_" + key, nm, hw, "no prediction" if key == "unifire" else "no hit (%s)" % cutoff, 0])
        elif key == "diamond":
            acc, pn, gene, org, pid_, qcov, e = v
            out.append(["w_" + key, nm, hw, "%s %s%s, <i>%s</i>" % (uni_link(acc), esc(pn), " (%s)" % esc(gene) if gene else "", esc(org)), 1,
                        "E-value %s%s · whole protein, no position" % (esc(e), "" if pid_ is None else
                                                                        ", %.0f%% identity over %.0f%% of your protein" % (pid_, qcov or 0))])
        elif key == "jackhmmer":
            fam, seed, rnd, e = v
            out.append(["w_" + key, nm, hw, "found from the seed %s (%s)" % (esc(seed), esc(fam)), 1,
                        "E-value %s, round %s · whole protein, no position" % (esc(e), rnd)])
        elif key == "fantasia":
            acc, pn, gene, org, dist = v
            counts_ = dist < FANTASIA_MAX_DIST
            out.append(["w_" + key, nm, hw, "%s %s%s, <i>%s</i> · its most similar protein" % (
                uni_link(acc), esc(pn), " (%s)" % esc(gene) if gene else "", esc(org)), 1 if counts_ else 2,
                "embedding distance %.2f; this site counts FANTASIA's pick only below %.2f · whole protein, no position" % (dist, FANTASIA_MAX_DIST)])
        elif key == "unifire":
            what, n = v
            out.append(["w_" + key, nm, hw, esc(what), 1, "%d prediction%s · whole protein, no position" % (n, "" if n == 1 else "s")])
    return out


def dense_notrun(run_status, has_struct=False):
    """window.NOTRUN: track -> why it was not run on this protein."""
    out = {}
    for t in TRACKS:
        s = run_status[t["key"]]
        if s == "too long":
            out[t["key"]] = "not run: Biomni predicted structures only for proteins of %d aa or less%s" % (
                t["max_aa"], "; the structure at the right is a separate prediction" if has_struct else "")
        elif s == "not run":
            out[t["key"]] = ("not run on 3-D structures; Foldseek on shape letters predicted from the sequence is the next row"
                             if t["key"] == "foldseek" and run_status.get("prostt5") != "not run" else "not run: no %s in the input folder" % t["file"])
    return out


def dense_key(ctx, has_plddt, has_struct, extra=()):
    """The three-line key: colour = kind, mark = program, only the programs run on this protein; then
    identity, predictions, the seed, the line under a hit, and pLDDT when there is a structure."""
    kinds = {}
    for t in TRACKS:
        if ctx["run_status"][t["key"]] != "run":
            continue
        prog = PROGRAMS[t["program"]]
        names = kinds.setdefault(prog["cat"], {})
        names.setdefault(t["program"], []).append(t)
    for k, name in extra:  # protein.html's programs drawn from the whole-proteome tables
        kinds.setdefault(WHOLE_TRACKS[k]["cat"], {})[k] = [dict(program=k, name=name)]
    parts = ['<span class="lead">colour = kind, mark = program:</span>']
    for cat in ("seq", "struct", "kmer", "pfam"):
        if cat not in kinds:
            continue
        items = []
        for prog, ts in kinds[cat].items():
            nm = (PROGRAMS.get(prog) or WHOLE_TRACKS[prog])["name"]
            if len(ts) > 1:
                nm = "%s (%s)" % (nm, " and ".join(t["name"].replace(PROGRAMS[prog]["name"], "").strip(", ") or "one round" for t in ts)) \
                    if prog == "mmseqs2" else "Foldseek (on 3-D structures, and ProstT5 + Foldseek)"
            else:
                nm = ts[0]["name"]
            items.append("%s%s" % (svg_inline(legend_mark(ts[0], 12, 10), 12, 10), esc(nm)))
        parts.append('<span class="kind"><b>%s:</b> %s</span>' % (esc({"seq": "sequence search", "struct": "structure", "kmer": "k-mer",
                                                                       "pfam": "domain scan"}[cat]), " ".join(items)))
    parts.append('<span>fainter = fewer identical residues: <i class="ramp"><i style="opacity:.3"></i><i style="opacity:.5"></i>'
                 '<i style="opacity:.75"></i><i style="opacity:1"></i></i> &lt;20% · 20–40% · 40–60% · 60%+</span>')
    if ctx["run_status"].get("interproscan") == "run":
        parts.append('<span><svg width="12" height="10"><rect x="1" y="1" width="10" height="8" rx="2" fill="none" stroke="#3D3D3D" '
                     'stroke-width="1.5" stroke-dasharray="3 2"/></svg> prediction (disorder, coil, membrane, signal): shown, not counted</span>')
    if ctx["run_status"].get("kmerseek") == "run":
        parts.append('<span><svg width="18" height="10"><rect width="18" height="10" rx="2" fill="#C0392B" opacity=".3"/>'
                     '<rect x="4" y="3" width="10" height="4" fill="#b8860b"/></svg> seed inside a kmerseek box: where its search started</span>')
    parts.append('<span><svg width="26" height="10"><line x1="0" x2="26" y1="5" y2="5" stroke="#d6d5cf" stroke-width="2"/><rect x="2" y="2" '
                 'width="7" height="6" fill="none" stroke="#8a8a8a" stroke-width=".8"/><line x1="14" x2="22" y1="5" y2="5" stroke="#52514e" '
                 'stroke-width="3"/></svg> under a hit: the other protein, dark = matched part, boxes = its Pfam domains</span>')
    if has_plddt:
        sw = " ".join('<i class="pl" style="background:%s"></i>%s' % (c, "%d+" % lo if hi > 100 else "under %d" % hi if lo <= 0 else "%d–%d" % (lo, hi))
                      for lo, hi, c, _ in PLDDT)
        parts.append('<span class="kind" title="the predicted structure\'s confidence per residue, in the strip under the ruler%s"><b>pLDDT:</b> %s%s</span>'
                     % (" and on the 3-D model" if has_struct else "", sw,
                        ' <i class="pl" style="background:%s"></i>pointed-at hit on the model' % HIGHLIGHT if has_struct else ""))
    parts.append('<span class="kind" id="sskey"></span>')  # filled by dense-bhf.js when the page has s4pred's prediction
    parts.append('<span class="kind" id="alphakey"></span>')
    return '<div class="key">%s</div>' % "".join(parts)


def viewer_html(pid, st):
    """The 3-D model, sticky at the right of the rows (assets/structure.js draws it in Mol*)."""
    mean = sum(st["plddt"]) / max(1, len(st["plddt"]))
    n = len(st["plddt"])
    warn = ""
    if st["ca_ca"] is not None and not 3.6 <= st["ca_ca"] <= 4.0:
        warn = ('<p class="warn">These coordinates are not a physical protein chain: neighbouring residues\' central carbon atoms '
                '(CA) are a median %.1f &#8491; apart, where every real protein has 3.8 &#8491;. The ribbons follow the file, so their '
                'shape cannot be trusted.</p>' % st["ca_ca"])
    rects = "".join('<rect x="%d" y="0" width="1.05" height="10" fill="%s"/>' % (i, next((c for lo, hi, c, _ in PLDDT if lo <= v < hi), PLDDT[-1][2]))
                    for i, v in enumerate(st["plddt"]))
    strip = ('<svg id="pstrip" viewBox="0 0 %d 10" preserveAspectRatio="none" role="img" aria-label="pLDDT of each residue, residue 1 to %d">'
             '%s<rect id="pmask1" x="0" y="0" width="0" height="10" fill="#fff" opacity="0.7"/><rect id="pmask2" x="0" y="0" width="0" '
             'height="10" fill="#fff" opacity="0.7"/><rect id="pbrush" x="0" y="0" width="0" height="10" fill="none" stroke="%s" '
             'stroke-width="2.5" vector-effect="non-scaling-stroke"/></svg>' % (n, n, rects, HIGHLIGHT))
    made = (' The file says how it was made: <span class="mono">%s</span>. pTM, where given, is the predictor\'s confidence in the '
            'overall fold, 0 to 1; below 0.5 the arrangement of the parts is a guess even where each part is confident.' % esc(st["made"])) if st.get("made") else ""
    return ('<aside class="viewer" id="structure"><h2 id="struct-h">3-D structure (predicted); the hit you point at turns magenta</h2>'
            '<p class="note skey">Colours as in the key\'s pLDDT entry. Drag across the strip, or click a hit, to colour residues magenta on '
            'the model; the icons at the model\'s top right open Mol*\'s own panels. <b id="molsel" aria-live="polite"></b></p>%s%s'
            '<div id="mol" data-scale="1" data-n="%d" data-hl="%s"></div>'
            '<details class="more"><summary>How the structure was made and how to use the viewer</summary><p class="note">Ribbons made '
            'from <span class="mono">%s</span> (%d residues), coloured by pLDDT as in the strip above them (mean %.0f here). Point at a hit '
            'row, or drag across the strip, to colour those residues magenta on the structure and the strip; the rest turns grey. Clicking '
            'a row, or dragging across the strip, keeps the selection until you click the strip again. Point at the structure to see a '
            'residue\'s number and mark it on the strip. The viewer is Mol* (PDBe\'s build, the one the PDB and AlphaFold DB use): drag to '
            'turn, scroll to zoom, and its panels (the icons at its top right) change the display style, colours and selections.%s</p>'
            '</details></aside>' % (warn, strip, n, HIGHLIGHT, esc(st["file"]), n, mean, made))


# Which Swiss-Prot each program searched, from Biomni's note of 2026-11-03 (full_inputs/biomni_databases_2026-11-03.md)
# and the check of 2025_01 by accession on Sherlock (test_inputs/PROVENANCE.md). (track keys, program, database, part, proteins)
SWISSPROT_DENSE = [
    (("kmerseek",), "kmerseek", "Swiss-Prot release 2026_03 (2 September 2026)", "all of it", "575,748"),
    (("mmseqs2_3rounds", "phmmer", "jackhmmer"), "MMseqs2, three rounds; phmmer; jackhmmer",
     "Swiss-Prot release 2025_01 (5 February 2025), 572,970 entries; every one of the 572,850 accessions is in 2025_01, which has "
     "exactly 120 more (checked on Sherlock, 2026-11-03)", "without its 120 sea-squirt (Ascidiacea) proteins", "572,850"),
    (("mmseqs2",), "MMseqs2 (Biomni's run)", "Swiss-Prot release 2026_03, downloaded 2026-10-31 from UniProt's current release; Biomni "
     "says three rounds, but the settings were not saved with the table", "all of it, sea-squirt proteins included", "575,748"),
    (("w_diamond",), "DIAMOND (Biomni's run)", 'the same 2026_03 file; <span class="mono">--sensitive</span>, E ≤ 0.001, at most 25 hits per protein',
     "all of it", "575,748 in the database; 533,805 distinct proteins appear in its hit table"),
    (("prostt5", "foldseek"), "ProstT5 + Foldseek; Foldseek on structure (Biomni's runs)",
     "Foldseek's prebuilt AlphaFold/Swiss-Prot structure database, downloaded 2026-10-31; more models than Swiss-Prot entries, "
     "because it keeps models of merged and retired accessions", "sea-squirt proteins kept", "590,183 (Biomni's corrected count; an earlier version said 570,421)"),
    (("w_jackhmmer",), "jackhmmer with immune-gene seeds (Biomni's run)", "the Botryllus proteome itself: the immune-gene seed proteins "
     "were the queries and Botryllus the database, the reverse of every other row", "all 45,339 Botryllus proteins", "45,339"),
]


def dense_swissprot(ctx, whole):
    run = {t["key"] for t in TRACKS if ctx["run_status"][t["key"]] == "run"} | {w[0] for w in whole}
    rows = "".join('<tr><td>%s</td><td>%s</td><td>%s</td><td class="num">%s</td></tr>' % (esc(n), d, esc(p), esc(c))
                   for keys, n, d, p, c in SWISSPROT_DENSE if run & set(keys))
    return ('<section class="below" id="swissprot"><h2>Which Swiss-Prot each program searched</h2><p class="note">The programs did not '
            'all search the same copy of Swiss-Prot, so a protein can be missing from one program\'s hits only because it was not in that '
            'copy. Names, genes and species of the matched proteins come from Swiss-Prot 2026_03, or from UniProt\'s website on 2026-11-01 '
            'for proteins not in it.</p><div class="tablewrap"><table><thead><tr><th>Program</th><th>Database searched</th><th>Part '
            'searched</th><th class="num">Proteins</th></tr></thead><tbody>%s</tbody></table></div><p class="note">Source for the Biomni rows: '
            'Biomni\'s report to Olga, 2026-11-03, checked by it against UniProt\'s release notes; the 2025_01 row was also checked here by '
            'comparing accessions. Biomni is the agent that ran the whole-proteome programs (MMseqs2, DIAMOND, ProstT5 + Foldseek, Foldseek '
            'on structure, jackhmmer with immune seeds, FANTASIA, UniFIRE).</p></section>' % rows)


def dense_help(ctx):
    """The ? Help panel: what a row is, the controls, the keys, the lanes, and each program's rule,
    setting, database and box."""
    progs = []
    for t in TRACKS:
        k = t["key"]
        if ctx["run_status"][k] != "run":
            continue
        c, d = fmt_thr(COUNT_E), fmt_thr(E_MAX)
        thr = ("counts at Prob ≥ 50 and E-value at most %s; shown up to E-value %s" % (c, d) if k == "hhsearch" else
               "counts at E-value at most %s where there is one, shown up to %s; disorder, coiled-coil, membrane-span and signal-peptide "
               "predictions are shown dashed and do not count" % (c, d) if k == "interproscan" else
               "counts at %s at most %s; shown up to %s" % (t["score"], c, d))
        coords = {"domtbl_seq": "hmm from-to (positions on your protein)",
                  "domtbl_pfam": "env from-to, the domain's envelope: the stretch hmmscan assigns to the domain, a little wider than the aligned part",
                  "m8": "qstart-qend", "kmerseek": "region_start-region_end, plus 1 at the start", "hhr": "Query HMM range",
                  "interproscan": "start-stop"}[t["kind"]]
        progs.append("<dt>%s</dt><dd>%s</dd>" % (esc(t["name"]), "<br>".join(esc(x) for x in (thr, ctx["settings"][k], ctx["databases"][k], "box: " + coords))))
    return ('<aside class="help" id="help"><h2>What one row is</h2><p>A <b>program row</b> is folded: its track is a histogram of how many '
            'of its hits cover each residue, with its best hit shown as a box on top. The columns give the count of hits at the E-value '
            'chosen, the best E-value, that hit\'s score (bits, or hhsearch\'s Prob) and identity where known, and the best hit\'s name. '
            'Click the row to unfold one row per hit, best E-value first.</p><p>A <b>hit row</b> gives the residues of your protein it '
            'covers, its E-value and score, the identity, and the matched part of the other protein. Click it to see the residues: a map of '
            'every hit between the two proteins, and the alignment with the same residues in bold. For a kmerseek hit the alignment also '
            'shows each residue\'s class in the alphabet kmerseek used and the seed where its search started.</p><p>The jackhmmer '
            'with immune seeds and FANTASIA rows come from Biomni\'s run on the whole proteome (Biomni is the agent that ran the '
            'whole-proteome programs except kmerseek). That run kept each protein\'s best hit and not where it lies, so they name the hit with '
            'no box, and each counts in the sentence by its own cutoff, which the E-value chips do not change. DIAMOND\'s table kept no positions '
            'either, so its best hit is aligned again and shown, marked realigned. UniFIRE\'s rules work on signatures, so each prediction is '
            'placed where the rule\'s signatures lie on this protein; a site uses the rule\'s own position.</p>'
            '<h2>Controls</h2><p>The <b>E-value cutoff</b> chips keep the hits at or below that value; the sentence at the top and every '
            'count follow it. <b>Search</b> looks through the matched protein\'s name, gene, accession, species, Pfam families and GO terms; '
            'every word typed must appear. The 3-D structure, where there is one, colours the residues of the hit you point at; a click '
            'keeps them until you click the pLDDT strip.</p><h2>Keyboard</h2><dl><dt><kbd>j</kbd> <kbd>k</kbd></dt><dd>next and previous '
            'row</dd><dt><kbd>Enter</kbd> <kbd>Space</kbd></dt><dd>unfold a program, or open a hit\'s residues</dd><dt><kbd>Esc</kbd></dt>'
            '<dd>close every open hit and this help</dd><dt><kbd>s</kbd> <kbd>t</kbd></dt><dd>scroll the 3-D model into view; back to the '
            'top</dd><dt><kbd>?</kbd></dt><dd>this help</dd></dl><h2>Lanes under the ruler</h2><p>The row "helix and strand, predicted", where s4pred was run: its prediction of helix (tall dark '
            'bar), strand (lower grey bar) and coil (the thin line) for each residue, made from this protein\'s sequence alone; a prediction, '
            'not a measured structure. In the ruler\'s own row: pLDDT, where there is a structure: the '
            'predictor\'s confidence per residue. Then your protein rewritten in a kmerseek alphabet, one colour per class (the alphabet of '
            'the kmerseek hit under the pointer; by default its best hit\'s): a long run of one colour is a stretch kmerseek matches on '
            'composition alone. Inside a kmerseek hit\'s box, the dark bar is its seed: a run of at least k residues whose classes are the '
            'same in both proteins at every position, where the search started.</p><p>Pointing at any hit fades every hit that does not '
            'overlap it on your protein. A program whose only output is a prediction (InterProScan\'s disorder call) or a named protein that '
            'falls short of its cutoff (FANTASIA) keeps its row; programs that returned nothing share one line at the bottom of the grid.</p>'
            '<h2>Each program: what counts as a hit, its setting, its database, what the box spans</h2><dl>%s</dl></aside>' % "".join(progs))


DENSE_FILTERS = (
    '<div class="filters"><input id="emax" type="hidden" value="%d" data-start="%d"><span class="lead2">E-value cutoff:</span>'
    '<span class="ladder" id="ladder" title="The E-value setting: how the sentence at the top changes as the cutoff tightens. Click one."></span>'
    '<span class="search"><svg width="13" height="13" viewBox="0 0 16 16" aria-hidden="true"><circle cx="7" cy="7" r="5" fill="none" '
    'stroke="#898781" stroke-width="1.6"/><line x1="11" y1="11" x2="15" y2="15" stroke="#898781" stroke-width="1.6"/></svg>'
    '<input id="hq" type="search" placeholder="protein, gene, species, Pfam family or GO term" autocomplete="off" spellcheck="false" '
    'aria-label="Search the hits: every word must match"></span>'
    '<button class="btn" id="expand-all" title="Show every hit of every program">&#9662; Expand all</button>'
    '<button class="btn" id="collapse-all" title="Fold every program to one row">&#9656; Collapse</button>'
    '<button class="btn" id="dl" title="The hits shown, one line each, as a tab-separated file">&#8595; TSV of the hits shown</button>'
    '<button class="btn" id="help-btn" aria-expanded="false" title="What the columns and marks mean; keys">? Help</button>'
    '<span class="count" id="count"></span><p class="note" id="overlap" aria-live="polite"></p>'
    '<button class="btn" id="clear" title="E-value back to 10, search and domain filter cleared, sort by E-value">&#215; clear all</button></div>')

DENSE_COLHEAD = (
    '<div class="g colhead"><div>program, then its hits when unfolded</div><div class="num"><button data-sort="aa" title="program row: hits at '
    'this cutoff, distinct proteins, settings, peak coverage; hit row: residues of your protein, their count, and the two most common residues '
    'with their share. Click to sort hits by start">hits · residues <i></i></button></div><div class="num"><button data-sort="e" '
    'title="program row: best E-value, and one tick per hit on a log axis from 0.001 to 10; hit row: its E-value. Click to sort">E-value '
    '<i></i></button></div><div class="num"><button data-sort="score" title="score in bits; on hhsearch\'s rows its Prob, the chance out of 100 that the hit is a true relative. On a kmerseek hit row: the '
    'alphabet and k it was found with. Click to sort">score (bits or Prob) · setting <i></i></button></div><div class="num"><button data-sort="id" '
    'title="identical residues: the program\'s own figure (ProstT5, MMseqs2) or counted on its alignment. Click to sort">identical <i></i>'
    '</button></div><div><button data-sort="tl" title="program row: its best hit; hit row: the matched part of the other protein. Click to '
    'sort by the other protein\'s length">best hit · matched part <i></i></button></div><div class="track" title="Boxes sit at the residues '
    'of your protein the program reported. Under each hit, the other protein is shown at its own length with the matched part dark and its '
    'Pfam domains as ticks. Click a row for the residues.">where on your protein · under each hit, the other protein to its own length</div></div>')


def dense_page(pid, pname, L, ctx, st, whole, page_data):
    """The page's HTML, in the design described above."""
    title = ("the ID as written in the proteome FASTA: gene ID, underscore, transcript ID" +
             ("; for this protein the two are the same" if pid.split("_")[0] == pid.split("_", 1)[-1].split(".")[0] else ""))
    pre = "../assets/"
    head = page_head("%s%s: what every annotation program found" % (pname + " (" if pname else "", pid + (")" if pname else "")), 1)
    head = head.replace("</head>", '<link rel="stylesheet" href="%sdense-bhf.css%s">%s</head>' % (
        pre, ASSET_VER[0], '<link rel="stylesheet" href="%smolstar/pdbe-molstar-light.css">' % pre if st else ""), 1)
    out = [head,
           '<div class="titleblock"><div class="titlebar"><a class="back" href="../index.html" title="back to the table of all proteins">'
           '<span aria-hidden="true">←</span> All proteins</a><h1 title="%s">%s</h1>%s<span class="len">%s amino acids</span>'
           '</div><div class="sentence"><p id="sent"></p><div id="sentnote"></div></div></div>'
           % (esc(title), esc(pid), '<span class="pn">%s</span>' % esc(pname) if pname else "", "{:,}".format(L)),
           '<header class="top">', DENSE_FILTERS % (round(math.log10(SHOW_E)), round(math.log10(COUNT_E))),
           dense_key(ctx, bool(page_data.get("plddt")), bool(st), (("diamond", "DIAMOND"), ("unifire", "UniFIRE")) if ctx["biomni"] else ()), DENSE_COLHEAD, '</header>',
           '<main><div class="withviewer%s"><div id="grid"></div>%s</div>' % ("" if st else " noviewer", viewer_html(pid, st) if st else ""),
           dense_swissprot(ctx, whole), '</main>', ctx["footer"],
           '<div id="tip" role="tooltip" hidden></div>', dense_help(ctx)]
    glob = "window.WHOLE=%s;window.NOTRUN=%s;window.HOW=%s;" % (
        json.dumps(dense_whole(ctx["biomni"].get(pid))), json.dumps(dense_notrun(ctx["run_status"], bool(st))),
        json.dumps({t["key"]: (ctx["settings"][t["key"]] if ctx["run_status"][t["key"]] == "run" else "not run") for t in TRACKS}))
    out.append('<script>const PAGE=%s;%s</script><script src="%sdense-bhf.js%s"></script>'
               % (json.dumps(page_data, separators=(",", ":")).replace("</", "<\\/"), glob.replace("</", "<\\/"), pre, ASSET_VER[0]))
    if st:
        out.append('<script src="%smolstar/pdbe-molstar-plugin.js"></script><script src="%s.pdb.js"></script><script src="%sstructure.js%s"></script>'
                   % (pre, safe(pid), pre, ASSET_VER[0]))
    out.append("</body></html>")
    return "".join(out)


def identity_count(k, h, seq, tseq, a):
    """[identical, aligned columns] of one hit, for the identity column: kmerseek's gapless pair; the
    alignment the program wrote (a lower-case letter of hmmscan's model consensus counts when it is
    the same residue); otherwise a local realignment of the two ranges (BLOSUM62, gap open -11,
    extend -1), skipped when the two ranges together are too long to align quickly. None if neither."""
    if k == "kmerseek" and tseq and h.get("t0"):
        q = seq[h["q0"] - 1:h["q1"]]
        t = tseq[h["t0"] - 1:h["t1"]]
        return [sum(1 for x, y in zip(q, t) if x == y), min(len(q), len(t))]
    if a is None and tseq and h.get("t0") and k not in ("hmmscan", "interproscan"):
        lo, hi = min(h["t0"], h["t1"]), max(h["t0"], h["t1"])
        if (h["q1"] - h["q0"] + 1) * (hi - lo + 1) <= 4_000_000:
            a = recompute(seq[h["q0"] - 1:h["q1"]], h["q0"], tseq[lo - 1:hi], lo)
    if not a:
        return None
    cols = [(x, y) for x, y in zip(a["q"], a["t"]) if x.isalpha() or y.isalpha()]  # every alignment column, gaps included
    return [sum(1 for x, y in cols if x.isalpha() and x.upper() == y.upper()), len(cols)] if cols else None


def twelve_text(ctx, whole):
    """The hover on the sentence's denominator: which programs it counts."""
    notes = []
    for prog in PROGRAM_ORDER:
        ts = [t for t in TRACKS if t["program"] == prog]
        if len(ts) < 2:
            continue
        run = [t for t in ts if ctx["run_status"][t["key"]] == "run"]
        if len(run) == len(ts):
            notes.append("%s's two runs counted once" % PROGRAMS[prog]["name"])
        for t in ts:
            if t not in run:
                notes.append("%s was not run%s" % (t["name"], " on this protein, which is over %d aa" % t["max_aa"]
                                                    if ctx["run_status"][t["key"]] == "too long" else ""))
    ran = whole  # every whole-protein program run here counts, hit or not (dense_whole leaves out the ones not run)
    s = "the %d programs that place hits along the protein%s" % (ctx["n_run"], " (%s)" % "; ".join(notes) if notes else "")
    if ran:
        s += " plus the %d whole-protein program%s from Biomni's run: %s" % (len(ran), "" if len(ran) == 1 else "s",
                                                                             ", ".join(WHOLE_DENSE[w[0][2:]][0] for w in ran))
    return s


def build_protein(task):
    """Runs in a worker: draws one protein page, returns its row for the index."""
    pid, ctx_key = task
    ctx = CTX
    p = ctx["proteins"][pid]
    hits = [materialize(x) for x in ctx["hits"].get(pid, [])]
    L = p["length"]
    long_for = {t["key"]: "too long" for t in TRACKS
                if t.get("max_aa") and L > t["max_aa"] and ctx["run_status"][t["key"]] == "run"}
    if long_for:  # tracks run only on shorter proteins: "not run" on this page
        ctx = dict(ctx, run_status={**ctx["run_status"], **long_for})
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
    tdom = {}  # target -> its Pfam domains, once per page, not once per hit (TTN has 10,000+ hits)
    tsrc = {}  # target -> where its Pfam domains come from (interproscan, hmmscan or interpro)
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
            if h.get("feature"):
                d["x"] = 1  # a feature prediction: drawn as a dashed box, never counted
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
                elif h.get("desc"):
                    d["n"] = h["desc"]
            d["f"] = None if h["_names"] is None else [family_full(n, ctx["names"]) for n in sorted(h["_names"])]
            bit = 1 << PROGRAM_ORDER.index(t["program"])
            for txt in [name, h["target"], d.get("n"), d.get("g"), d.get("o")] + (d["f"] or []):
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
                elif h.get("tlen"):
                    d["tl"] = h["tlen"]
                doms = ctx["target_ann"].get(h["target"])
                if doms and h["target"] not in tdom:
                    tdom[h["target"]] = non_overlapping(doms, lambda nm: family_label(sorted(nm)[0], ctx["names"]) if nm else "")
                if doms is not None:  # where the target's Pfam domains come from, also once per page
                    tsrc[h["target"]] = ctx["target_src"].get(h["target"])
                if t["key"] == "hhsearch":
                    a = dict(h["aln"], src="the alignment hhsearch wrote (the .hhr file)") if h.get("aln") else None
                elif t["key"] == "kmerseek":
                    d["m"] = h.get("moltype")
                    d["kk"] = h.get("ksize")
                    tab = ctx["alphabets"].get(d["m"])
                    if tab and tseq and d["kk"]:
                        # the seed: k columns in the same class of the hit's alphabet in both
                        # (kmerseek's hits are ungapped, so column i is residue q0 + i and t0 + i)
                        runs, n = seed_runs(seq[h["q0"] - 1:h["q1"]], tseq[h["t0"] - 1:h["t1"]], tab, d["kk"])
                        d["sd"] = [[h["q0"] + i, h["q0"] + j] for i, j in runs]
                        d["sn"] = n
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
            if not d.get("x"):
                ic = identity_count(t["key"], h, seq, ctx["targets"].get(h.get("target")) if t["key"] not in ("hmmscan", "interproscan") else None, a)
                if ic:
                    d["id"] = ic
            data.append(d)
        by_track[t["key"]] = hs

    plddt = ctx["plddt"].get(pid)
    st = ctx["structures"].get(pid)
    W = int(max(640, min(L * 2.0, 3200)))
    px = W / L
    whole = whole_rows(ctx["biomni"].get(pid), W)
    n_whole = sum(1 for w in whole if w[3])
    m_whole = sum(1 for w in whole if w[4])  # WHOLE programs run on this protein
    pos_e = [h["evalue"] for v in by_track.values() for h in v if h["evalue"]]
    c_pos = int(round(math.log10(COUNT_E)))
    s_pos = int(round(math.log10(SHOW_E)))
    lo = min(c_pos, s_pos, max(-300, math.floor(math.log10(min(pos_e))))) if pos_e else min(c_pos, s_pos)
    steps = sentence_steps(hits, ctx, lo, n_whole, m_whole)
    agreement(hits, ctx["target_ann"])  # AGREE_INFO at the counting cutoff, for the brackets under the ruler
    agree = [dict(acc=a["acc"], name=ctx["names"].get(a["acc"]) or a["acc"], a=a["a"], b=a["b"],
                  programs=sorted(PROGRAMS[p]["name"] for p in a["programs"]), n=a["n"], hits=a["hits"]) for a in AGREE_INFO[0][:3]]  # the bracket lane: the named family and the next two
    ptm = re.search(r"pTM ([0-9.]+)", st["made"]) if st and st.get("made") else None
    page_data = dict(
        id=pid, pn=pname or "", q=seq, px=px, w=W, top=TOP_N, m=ctx["n_run"], sent=steps, hits=data, td=tdom, ds=tsrc, names=score_names,
        alpha={m: ctx["alphabets"][m] for m in sorted({d["m"] for d in data if d.get("m")}) if m in ctx["alphabets"]},
        tracks={t["key"]: dict(name=t["name"], color=PROGRAMS[t["program"]]["color"], shape=PROGRAMS[t["program"]]["shape"],
                               cat=PROGRAMS[t["program"]]["cat"], p=t["program"]) for t in TRACKS},
        pfam_color=CATEGORIES["pfam"]["color"], tseq_files=ctx["n_shards"],
        agree=agree, twelve=twelve_text(ctx, dense_whole(ctx["biomni"].get(pid))))
    if plddt:
        page_data["plddt"] = [round(v, 2) for v in plddt]
    if ptm:
        page_data["ptm"] = ptm.group(1).rstrip(".")
    if st:  # the structure's text in its own file, read by assets/structure.js after Mol* (PDBe's pdbe-molstar, bundled in assets/molstar)
        (ctx["out"] / "proteins" / (safe(pid) + ".pdb.js")).write_text("window.PDB=%s;" % json.dumps(st["text"]))
    (ctx["out"] / "proteins" / (safe(pid) + ".html")).write_text(dense_page(pid, pname, L, ctx, st, dense_whole(ctx["biomni"].get(pid)), page_data))
    (ctx["out"] / "proteins" / (safe(pid) + ".aln.js")).write_text(
        "window.ALN=%s;" % json.dumps({str(k): v for k, v in alns.items()}, separators=(",", ":")))
    return dict(id=pid, length=L, n=len(progs) + n_whole, progs=progs, status=status,
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
               'tool\'s hits are shown on the same ruler of residue positions. Each program counts once, even when it was '
               'run twice (MMseqs2 one round and three rounds; Foldseek on predicted structures and ProstT5 + Foldseek), '
               'so there are 8 programs for 10 rows; %d of the 8 were run on this proteome. A hit counts when its E-value is at most %s '
               '(hhsearch: also Prob at least 50); the protein pages also show weaker hits, up to E-value %s. "Agree" means two hits from different programs overlap by '
               'at least half of the shorter one and name the same Pfam family or InterPro entry. The small '
               'picture is the "programs with a hit" line of each protein page: its height at each point is the '
               'number of programs (0 to %d) with a hit over that part of the protein, squeezed to the same 90-pixel width '
               'and the same 0-to-%d height for every protein whatever its length; hover over it for the position and count. "Best E-value, any tool" is the lowest E-value any tool printed '
               '(each tool computes its own, so they are not on one scale); 0 means the tool printed 0.</p>'
               '<div class="search"><label for="q">Search</label><input id="q" type="search" autocomplete="off" '
               'placeholder="a gene, protein or Pfam family, e.g. HLA, immune, MHC, RING"></div>'
               '<p class="note">The search looks at every hit shown on the protein pages (E-value at most %s): the name, gene and species of '
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
                       r["length"], "{:,}".format(r["length"]), r["n"], r["n"], acc_html(r["agreed"] or "none"),
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
    ap.add_argument("--show-at", type=float, default=SHOW_E,
                    help="the E-value slider on each protein page starts here (a power of 10, at most 10; default 10)")
    a = ap.parse_args()
    lg = math.log10(a.count_at) if a.count_at > 0 else None
    if lg is None or abs(lg - round(lg)) > 1e-9 or a.count_at > E_MAX:
        print("--count-at must be a power of 10 no larger than %g (the slider moves in steps of 10 times)" % E_MAX)
        sys.exit(2)
    globals()["COUNT_E"] = a.count_at
    ls = math.log10(a.show_at) if a.show_at > 0 else None
    if ls is None or abs(ls - round(ls)) > 1e-9 or a.show_at > E_MAX:
        sys.exit("--show-at must be a power of 10, at most %s" % fmt_thr(E_MAX))
    globals()["SHOW_E"] = a.show_at
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
            files = sorted(p.glob("*.hhr")) + sorted(p.glob("*.hhr.gz"))
            size = sum(f.stat().st_size for f in files)
            mt = max((f.stat().st_mtime for f in files), default=p.stat().st_mtime)
            return (name + "/ (%d .hhr files)" % len(files), size, datetime.datetime.fromtimestamp(mt).strftime("%Y-%m-%d %H:%M"))
        st = p.stat()
        return (p.name, st.st_size, datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M"))

    names_read = ["proteins.fasta", "swissprot.fasta"] + [t["file"] for t in TRACKS] + [
        "phmmer.out", "jackhmmer.out", "hmmscan.out", "protein_names.tsv", "plddt.tsv", "expression.tsv", "leaderboard.tsv", "run_settings.tsv", "kmerseek_alphabets.rs",
        TARGET_PFAM_FILE, "biomni_best.json"]
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
            for f in sorted(p.glob("*.hhr")) + sorted(p.glob("*.hhr.gz")):  # gzipped in the test set
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
    # Pfam families of targets: InterProScan's output on that target if there is one, else hmmscan's,
    # else InterPro's Pfam matches read from its REST API (target_pfam_interpro.tsv)
    target_src = {acc: "interproscan" for acc in target_ann}
    for acc, doms in hmmscan_target.items():
        if acc not in target_ann:
            target_ann[acc] = doms
            target_src[acc] = "hmmscan"
    if inp(TARGET_PFAM_FILE).exists():
        api_doms, api_names = read_target_pfam(inp(TARGET_PFAM_FILE))
        for acc, doms in api_doms.items():
            if acc not in target_ann:
                target_ann[acc] = doms
                target_src[acc] = "interpro"
        for pf, nm in api_names.items():
            names.setdefault(pf, nm)
        print("  %s: Pfam domains of %d more Swiss-Prot proteins" % (TARGET_PFAM_FILE, sum(v == "interpro" for v in target_src.values())), flush=True)

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
    # Predicted structures, structures/<protein ID>.pdb (or the ID's first part, FUN008084.pdb);
    # optional. Their confidence fills the pLDDT strip when plddt.tsv has none for that protein.
    structures = {}
    if (D / "structures").is_dir():
        for pid in prots:
            for name in dict.fromkeys([pid, pid.split("_")[0]]):
                f = D / "structures" / (name + ".pdb")
                if f.exists():
                    structures[pid] = read_structure(f)
                    if pid not in plddt:
                        plddt[pid] = structures[pid]["plddt"]
                    break
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
    # Biomni's results with no position (test_inputs/make_biomni_best.py); optional
    bfile = D / "biomni_best.json"
    biomni = json.loads(bfile.read_text()) if bfile.exists() else {}
    n_whole = len(WHOLE) if biomni else 0
    built = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    CTX.update(dict(proteins=prots, prot_names=prot_names, hits=hits, target_ann=target_ann, target_src=target_src, names=names, entry_names=entry_names,
                    target_info=target_info,
                    targets=targets, run_status=run_status, settings=settings, databases=databases,
                    leaderboard=leaderboard, plddt=plddt, structures=structures, expression=expression, alphabets=alphabets,
                    out=out, footer="", program_run=program_run, n_run=n_run,
                    biomni=biomni, n_whole=n_whole))
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
    summary = dict(count_at=COUNT_E, drawn_up_to=E_MAX, programs_run=n_run, programs_total=8 + n_whole,
                   total=len(rows),
 **cats, only_program=only,
                   run_status=run_status, seconds=round(time.time() - t_start, 1))

    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
