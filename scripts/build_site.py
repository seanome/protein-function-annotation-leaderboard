#!/usr/bin/env python3
"""Build the static site in docs/ from the tables in data/.

Every number on every page comes from data/results.parquet or data/compute.parquet. A
filter combination with no row has no entry in the JSON the pages read, and the pages show
"—" with "not run" for it. Nothing is averaged, filled or estimated here except
cost_usd_est, whose formula and prices are printed beside it.
"""

import argparse
import html
import json
import math
import re
import sys
from pathlib import Path

import markdown
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
ROOT_DIR = Path(__file__).resolve().parents[1]
from tools import KINDS, SPECIES, TOOLS  # noqa: E402

SITE = "Protein annotation leaderboard"
REPO = "seanome/protein-function-annotation-leaderboard"

# Default prices for cost_usd_est. Both from AWS's public on-demand price list for
# US East (N. Virginia), Linux, read 2026-10-15 from the file the price page loads
# (publication date 2026-10-12): c7i.large $0.08925 per hour for 2 vCPU, and g5.xlarge
# $1.006 per hour for 1 NVIDIA A10G GPU.
PRICE_PAGE = "https://aws.amazon.com/ec2/pricing/on-demand/"
PRICE_FILE = (
    "https://b0.p.awsstatic.com/pricing/2.0/meteredUnitMaps/ec2/USD/current/"
    "ec2-ondemand-without-sec-sel/US%20East%20(N.%20Virginia)/Linux/index.json"
)
R_CPU = 0.08925 / 2
R_GPU = 1.006

TRUTH = {"swissprot": "Swiss-Prot features", "pfam": "Pfam", "pfamn": "Pfam-N"}
DEFAULTS = {
    "truth": "pfam", "calls": "raw", "species": "zebrafish", "axis": "all",
    "stratum": "all", "x": "cpu", "y": "fmax", "sort": "fmax", "dir": "desc",
    "view": "best", "r_cpu": R_CPU, "r_gpu": R_GPU,
}
# Divergence times from human, in million years (Mya), read by main() from results.parquet's
# target_mya_to_human (TimeTree 5 pairwise ages). E. coli, for which TimeTree has no age, is
# absent. main() fills in DEFAULT_WORDS from it too.
SPECIES_MYA: dict = {}
DEFAULT_WORDS = "Pfam truth set, raw calls, zebrafish as the relative searched"


def mya_label(v: float) -> int:
    """A divergence time as the site shows it: the nearest whole million years."""
    return int(round(v))


def rel_age(sp: str) -> str:
    """A relative's divergence time from human in words: "435 Mya", or that TimeTree has none."""
    return f"{mya_label(SPECIES_MYA[sp]):,} Mya" if sp in SPECIES_MYA else "no TimeTree age"


def rel_order_words(rel_order: list[str]) -> str:
    """Nearest and farthest relative by divergence time, and any with no age, which sort last."""
    aged = [sp for sp in rel_order if sp in SPECIES_MYA]
    no_age = [SPECIES[sp] for sp in rel_order if sp not in SPECIES_MYA]
    w = (f"from the nearest ({SPECIES[aged[0]]}, {rel_age(aged[0])}) to the farthest "
         f"({SPECIES[aged[-1]]}, {rel_age(aged[-1])})")
    return w + (f", then {' and '.join(no_age)}, which has no TimeTree age" if no_age else "")


# A subset of the proteins (a stratum in the tables) enters a card's comparisons and
# questions only with at least this many true instances; below it one or two calls swing
# Fmax from 0 to 1.
MIN_N = 30

IDENTITY = [("0-20%", "0–20%"), ("20-30%", "20–30%"), ("30-40%", "30–40%"),
            ("40-60%", "40–60%"), ("60-100%", "60–100%"),
            ("no_homolog", "no same-family domain")]
LENGTH_SHORT = [("1", "1 aa"), ("2-15", "2–15 aa"), ("16-30", "16–30 aa"), ("31-60", "31–60 aa")]
LENGTH_DOMAIN = [("61-120", "61–120 aa"), ("121-250", "121–250 aa"), ("251+", "251 aa or more")]
PLDDT = [("0-50", "0–50"), ("50-70", "50–70"), ("70-90", "70–90"), ("90-100", "90–100")]
FEATURES = {
    "ACT_SITE": "active site residue", "BINDING": "binding site residue",
    "COILED": "coiled coil", "DNA_BIND": "DNA-binding region", "DOMAIN": "domain",
    "INTRAMEM": "region inside the membrane that does not cross it", "MOTIF": "short motif",
    "REGION": "region of interest", "REPEAT": "repeat", "SITE": "other single-residue site",
    "TRANSMEM": "membrane-spanning segment", "ZN_FING": "zinc finger",
}
# The same features as short words, for chips and figure labels; the Swiss-Prot key is
# kept for the hover text and the JSON.
FEATURE_WORDS = {
    "ACT_SITE": "active site", "BINDING": "binding site", "COILED": "coiled coil",
    "DNA_BIND": "DNA-binding region", "DOMAIN": "domain", "INTRAMEM": "in-membrane region",
    "MOTIF": "short motif", "REGION": "region of interest", "REPEAT": "repeat",
    "SITE": "other site", "TRANSMEM": "membrane-spanning", "ZN_FING": "zinc finger",
}
SPECIES_GROUPS = [
    # Within a group, nearest first by TimeTree 5 age (mouse and rat tie at 84 Mya, cow and dog at 89).
    ("mammals", ["ptroglodytes", "mouse", "rnorvegicus", "btaurus", "clupus", "mdomestica"]),
    ("other vertebrates", ["chicken", "xtropicalis", "zebrafish"]),
    ("invertebrates", ["ciona", "fly", "worm"]),
    ("fungus, plant, bacterium", ["yeast", "arabidopsis", "ecoli"]),
    ("no relative", ["all"]),
]
# What the setting names mean. kmerseek alphabets are listed in the kmerseek README
# (seanome/kmerseek), with their residue groups.
ALPHABETS = {
    "protein20": "the 20 amino acids, not reduced",
    "gbmr4": "4 letters: ADKERNTSQ, YFLIVMCWH, G, P (Solis and Rackovsky 2000)",
    "wwmj5": "5 letters: CMFILVWY, ATH, GP, DE, SNQRK (Wang and Wang 1999)",
    "hp_thomas_dill2": "2 letters, hydrophobic ACFILMVWY and polar DEGHKNPQRST (Thomas and Dill 1996)",
}
SETTING_WORDS = {
    ("hmmer3_phmmer", "default"): "phmmer's default options, E-value cutoff 10",
    ("hmmer3_jackhmmer", "n3"): "three rounds of search (-N 3)",
    ("mmseqs2_seqseq", "s7"): "sensitivity 7 (-s 7), one round",
    ("mmseqs2_iterative", "s7"): "sensitivity 7 (-s 7), three rounds",
    ("hhblits", "single_seq"): "profiles built from single proteins, no family database",
    ("foldseek", "3di_aa"): "3Di shape letters and amino acids together",
    ("prostt5", "3di_from_seq"): "3Di letters predicted from sequence by ProstT5",
    ("hmmscan", "direct"): "scans Pfam's models; no relative searched",
}


# kmerseek alphabets in size groups, smallest first: 2 letters, 4-5 letters, then the 20
# amino acids unreduced. Every list, legend and panel of kmerseek settings uses this order.
ALPHA_ORDER = ["hp_thomas_dill2", "gbmr4", "wwmj5", "protein20"]
# Short plain names for settings, shown in tables and panels; the code name is the hover text.
SETTING_LABEL = {
    ("hmmer3_phmmer", "default"): "default options",
    ("hmmer3_jackhmmer", "n3"): "3 rounds",
    ("mmseqs2_seqseq", "s7"): "sensitivity 7",
    ("mmseqs2_iterative", "s7"): "sensitivity 7, 3 rounds",
    ("hhblits", "single_seq"): "single-protein profiles",
    ("foldseek", "3di_aa"): "shape letters + amino acids",
    ("prostt5", "3di_from_seq"): "shape letters from sequence",
    ("hmmscan", "direct"): "Pfam models",
}
KMER_SETTING = re.compile(r"(.+)_k(\d+)_lc(False|True)")


def setting_label(tool: str, setting: str) -> str:
    """"gbmr4, k = 21" for kmerseek (plus a note when low-complexity k-mers were removed)."""
    m = KMER_SETTING.fullmatch(setting)
    if tool == "kmerseek" and m:
        return f"{m[1]}, k = {m[2]}" + (", low-complexity removed" if m[3] == "True" else "")
    return SETTING_LABEL.get((tool, setting), setting)


def setting_key(tool: str, setting: str) -> tuple:
    """Sort key: kmerseek by alphabet size group, then k; other tools by setting name."""
    m = KMER_SETTING.fullmatch(setting)
    if tool == "kmerseek" and m:
        a = ALPHA_ORDER.index(m[1]) if m[1] in ALPHA_ORDER else len(ALPHA_ORDER)
        return (a, int(m[2]), m[3], setting)
    return (0, 0, "", setting)


def setting_words(tool: str, setting: str) -> str:
    if (tool, setting) in SETTING_WORDS:
        return SETTING_WORDS[(tool, setting)]
    m = re.fullmatch(r"(.+)_k(\d+)_lc(False|True)", setting)
    if tool == "kmerseek" and m:
        lc = "low-complexity k-mers kept" if m[3] == "False" else "low-complexity k-mers removed"
        return f"alphabet {m[1]} ({ALPHABETS.get(m[1], 'see the kmerseek README')}), k = {m[2]}, {lc}"
    return setting


# Marker shape for each kmerseek alphabet on Figure 3 (cost.html); other tools take theirs from TOOLS.
KMERSEEK_SHAPES = {"hp_thomas_dill2": "diamond", "gbmr4": "square", "wwmj5": "triangle",
                   "protein20": "circle"}

METRICS = [
    {"id": "fmax", "label": "Fmax", "field": "fmax"},
    {"id": "auprc", "label": "AUPRC", "field": "auprc"},
    {"id": "iou", "label": "Median IoU", "field": "iou"},
]
# The columns Figure 3 on cost.html can put on its axes. Scores run 0 to 1 on a linear
# axis, higher is better. Cost, time and memory span orders of magnitude and go on an
# axis where each step is ten times the last, lower is better. "less" finishes the
# sentence "no other tool is both ... and better" for the line through the best trade-offs;
# "short" names the quantity in the list of rows not plotted.
AXES = [
    {"id": "fmax", "label": "Fmax, 0 to 1", "short": "Fmax", "kind": "score", "better": "higher"},
    {"id": "auprc", "label": "AUPRC, 0 to 1", "short": "AUPRC", "kind": "score", "better": "higher"},
    {"id": "iou", "label": "Median IoU of correct calls, 0 to 1", "short": "median IoU", "kind": "score",
     "better": "higher"},
    {"id": "cpu", "label": "Core-hours per 1,000 queries, per relative searched", "short": "core-hours",
     "kind": "cost", "better": "lower", "less": "cheaper in core-hours"},
    {"id": "cost", "label": "Estimated cost, USD per 1,000 queries, per relative searched", "short": "cost",
     "kind": "cost", "better": "lower", "less": "cheaper"},
    {"id": "wall", "label": "Wall-clock hours per 1,000 queries, per relative searched", "short": "wall-clock time",
     "kind": "cost", "better": "lower", "less": "faster"},
    {"id": "ram", "label": "Peak RAM of the search, GB", "short": "peak RAM", "kind": "cost", "better": "lower",
     "less": "lighter on memory"},
    {"id": "db", "label": "Database build, core-hours per database built", "short": "database build time",
     "kind": "cost", "better": "lower", "less": "quicker to build its databases"},
]
# Definitions shown beside a column name or label (docs/assets/terms.js). The wording follows
# the Methods page and Table 1's text, so the site says one thing. Each links to its section of
# explainers.html. Values are HTML.
TERMS = {
    "fmax": ("Fmax", "explainers.html#fmax",
             "The best balance of precision and recall at any one score cutoff, with both worked "
             "out protein by protein and averaged. Precision is the share of calls that are correct; recall is "
             "the share of true annotations found. 0 to 1, higher is better."),
    "auprc": ("AUPRC", "explainers.html#auprc",
              "Area under the precision-recall curve: precision summed over every score cutoff, "
              "weighted by the recall gained there. It pools every call of every protein into one "
              "list. 0 to 1, higher is better."),
    "iou": ("Median IoU", "explainers.html#iou",
            "How well a correct call's ends match the true region. IoU is the residues the call and "
            "the region share divided by the residues either covers: 1 is a perfect match. The "
            "median is the middle IoU over all correct calls."),
    "n": ("True instances", "explainers.html#fmax",
          "The number of known annotations on the human proteins in this slice, whether or not the "
          "relative could supply them. Recall in Fmax is counted out of these."),
    "cpu": ("Core-hours per 1,000 queries", "explainers.html#cost",
            "Cores the search reserved × its wall time, for the search alone, divided by the number "
            "of relatives searched and scaled to 1,000 human proteins. A task that sat partly idle "
            "still counts all its cores. Less is cheaper."),
    "cpu_rep": ("CPU-h reported", "explainers.html#cost",
                "CPU-hours per 1,000 queries as reported by a tool's submitter, on hardware they "
                "name. None have been submitted yet."),
    "ram": ("Peak RAM", "explainers.html#cost",
            "The largest resident memory, in GB, of any of the search's tasks, from the run's task "
            "records. Memory is not priced in the cost."),
    "gpu": ("GPU-hours", "explainers.html#cost",
            "GPU time the search used, per 1,000 queries. No search here used a GPU. ProstT5 used "
            "one to predict 3Di shape letters while building its databases; that time is on its "
            "card."),
    "cost": ("Estimated cost", "explainers.html#cost",
             "Estimated cost in US dollars per 1,000 queries: core-hours × r_cpu + GPU-hours × "
             "r_gpu, where r_cpu and r_gpu are prices per hour. The default prices are AWS on-demand "
             "rates."),
    "calls": ("Raw or dedup", "explainers.html#dedup",
              "Raw keeps every call a tool makes. Dedup merges calls of the same family on the same "
              "protein that overlap with an IoU of 0.5 or more, keeping the highest-scoring one, so "
              "one domain hit by many of the relative's proteins counts once."),
    "wall": ("Wall-clock hours per 1,000 queries", "explainers.html#cost",
             "How long the search took from start to finish on the cores it was given (4 to 16), "
             "per relative searched, scaled to 1,000 human proteins. Less is faster."),
    "db": ("Database build", "explainers.html#cost",
           "Core-hours to build one relative's search database, done once before any search. "
           "The mean over the builds the run's task records hold. Not in the per-query cost."),
    "measured": ("measured", "explainers.html#cost",
                 "Measured on Stanford's Sherlock cluster, from the run's own task records."),
    "reported": ("reported", "explainers.html#cost",
                 "Reported by the tool's submitter, on the hardware they name. None yet."),
}


# The same terms in everyday words, for a reader who has never met them: (a few words shown
# under the column name, one or two plain sentences that open each definition). The fuller,
# technical definition from TERMS follows them, then the link to the worked example.
PLAIN = {
    "fmax": ("right calls, few wrong ones",
             "Does the tool find the real domains and features without adding wrong ones? Fmax "
             "balances the two at the tool's best score cutoff. 1 is perfect; 0 means it found none."),
    "auprc": ("right calls stay on top",
              "Line up every call the tool makes, most confident first. AUPRC asks how well the "
              "correct calls stay near the top of that list, whatever cutoff you pick. 1 is "
              "perfect; 0 means none are correct."),
    "iou": ("how well the ends line up",
            "For the calls that are right, how closely the start and end match the real domain. 1 "
            "means exactly the same residues; 0.5 means half of the residues the two cover "
            "together are shared."),
    "n": ("known annotations to find",
          "How many known domains or features the human proteins have here. A tool can find at "
          "most this many."),
    "cpu": ("computer time used",
            "How much computer time the search took for 1,000 human proteins against one relative. "
            "Lower is cheaper."),
    "cpu_rep": ("the authors' own figure",
                "Computer time the tool's authors report on their own hardware. None sent yet."),
    "ram": ("memory needed",
            "The most memory the search needed at one time. A computer with less memory than this "
            "cannot run it."),
    "gpu": ("graphics-card time",
            "Time on a graphics card (GPU). No search here used one."),
    "cost": ("price for 1,000 proteins",
             "What the search would cost, in US dollars, to annotate 1,000 proteins against one "
             "relative at cloud prices. You can change the prices. Lower is cheaper."),
    "calls": ("count each domain once?",
              "Raw counts every call. Dedup counts a domain once, even when several of the "
              "relative's proteins point to it."),
    "wall": ("time to finish",
             "How long you wait for the search to finish, for 1,000 proteins against one relative."),
    "db": ("set-up time",
           "Time to get the relative's proteins ready to search. It is done once, before any search."),
    "measured": ("timed by us",
                 "We timed it ourselves, on Stanford's Sherlock computer cluster."),
    "reported": ("timed by the authors",
                 "The tool's authors timed it on their own computers."),
}

def term_btn(term: str, where: str, up: str = "") -> str:
    """A "?" button and its definition, for terms.js. `where` makes the id unique on the page."""
    name, href, words = TERMS[term]
    simple = PLAIN[term][1]
    tid = f"tip-{term}-{where}"
    return (f'<button type="button" class="term" data-term="{term}" aria-label="What is {esc(name)}?" '
            f'aria-expanded="false" aria-controls="{tid}" aria-describedby="{tid}">?</button>'
            f'<span class="term-tip" role="tooltip" id="{tid}" hidden><b>{esc(name)}.</b> {simple}'
            f'<span class="term-detail"><b>In detail:</b> {words}</span> '
            f'<a href="{up}{href}">how it works →</a></span>')


STRATUM_WORDS = {
    "identity": dict(IDENTITY),
    "feature_length_bin": dict(LENGTH_SHORT + LENGTH_DOMAIN + [("whole_protein", "whole protein")]),
    "plddt": dict(PLDDT),
    "feature_type": dict(FEATURE_WORDS),
}
AXIS_WORDS = {
    "identity": "identity to the closest same-family domain in the relative",
    "feature_length_bin": "length of the annotated region",
    "plddt": "AlphaFold confidence (mean pLDDT) of the human protein",
    "feature_type": "Swiss-Prot feature type",
}


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def chiprows(species_mya: dict) -> list[dict]:
    def c(v, label, title=None):
        d = {"value": v, "label": label}
        if title:
            d["title"] = title
        return d

    return [
        {"id": "truth", "kind": "plain", "label": "Truth set",
         "help": "the known annotations calls are graded against. Pfam: protein domain "
                 "families. Pfam-N: regions a deep-learning model added where Pfam's own "
                 "models found nothing. Swiss-Prot features: 12 kinds of curated feature, "
                 "such as active sites and membrane-spanning segments",
         "groups": [{"chips": [c("swissprot", "Swiss-Prot features", "curated features from Swiss-Prot entries"),
                               c("pfam", "Pfam", "Pfam-A domain families"),
                               c("pfamn", "Pfam-N", "regions added by a deep-learning model where Pfam-A's models found nothing (Pfam 35.0)")]}]},
        {"id": "calls", "kind": "plain", "label": "Dedup or raw",
         "help": "raw keeps every call; dedup merges calls of one family that overlap on the "
                 "same protein, so one domain hit by many relatives' proteins counts once",
         "groups": [{"chips": [c("raw", "raw"), c("dedup", "dedup")]}]},
        {"id": "species", "kind": "plain", "label": "Nearest reviewed relative: time since split",
         "help": "the species whose annotated proteome is searched and whose annotations are "
                 "copied onto the human proteins; the number is how many million years (Mya) "
                 "ago its line split from ours, from TimeTree 5",
         "groups": [{"label": g, "chips": [
             c(s, f"{SPECIES[s]} {mya_label(species_mya[s])}") if s in species_mya
             else c(s, SPECIES[s], None if s == "all" else
                    "TimeTree 5 gives no divergence time for human and a bacterium")
             for s in members]} for g, members in SPECIES_GROUPS]},
        {"id": "identity", "kind": "stratum", "axis": "identity",
         "label": "Nearest reviewed relative: % identity",
         "help": "% identity of each human domain to the closest domain of the same family "
                 "in that relative",
         "groups": [{"chips": [c("all", "any")] + [c(v, w) for v, w in IDENTITY]}]},
        {"id": "length", "kind": "stratum", "axis": "feature_length_bin",
         "label": "Region type and length",
         "help": "length of the annotated region the call must find",
         "groups": [{"chips": [c("all", "any")]},
                    {"label": "short feature, 60 aa or less",
                     "chips": [c(v, w) for v, w in LENGTH_SHORT]},
                    {"label": "folded domain, 61 aa or more",
                     "chips": [c(v, w) for v, w in LENGTH_DOMAIN]},
                    {"chips": [c("whole_protein", "whole protein")]}]},
        {"id": "plddt", "kind": "stratum", "axis": "plddt", "label": "Disorder",
         "help": "mean AlphaFold confidence (pLDDT, 0–100) of the human protein; "
                 "low values mean disordered or unknown shape",
         "groups": [{"chips": [c("all", "any")] + [c(v, w) for v, w in PLDDT]}]},
        {"id": "ftype", "kind": "stratum", "axis": "feature_type", "label": "Feature type",
         "help": "the kind of Swiss-Prot feature; only the Swiss-Prot truth set has these",
         "groups": [{"chips": [c("all", "any")] + [c(k, FEATURE_WORDS[k], f"{k}: {v}") for k, v in FEATURES.items()]}]},
    ]


def mark_empty(rows: list[dict], res: pl.DataFrame) -> list[dict]:
    """Flag stratum chips that no row of the results table has, for any tool or filter."""
    present = {(a, v) for a, v in res.select("stratum_axis", "stratum").unique().iter_rows()}
    for row in rows:
        if row["kind"] != "stratum":
            continue
        for g in row["groups"]:
            for c in g["chips"]:
                if c["value"] != "all" and (row["axis"], c["value"]) not in present:
                    c["empty"] = True
    return rows


# The site's pages in reading order. The last two sit under "More" in the header.
NAV = [("index.html", "Your organism", "organism"), ("leaderboard.html", "Leaderboard", "leaderboard"), ("tools/index.html", "Tool cards", "tools"),
       ("cost.html", "Compare two columns", "cost"), ("explainers.html", "How the scores work", "explainers"),
       ("pileup.html", "Region pile-up", "pileup"), ("methods.html", "Methods", "methods")]
NAV_MORE = [("claims.html", "Claims", "claims"), ("coverage.html", "Coverage", "coverage")]


def layout(title: str, page: str, body: str, src: dict, depth: int = 0,
           head: str = "", description: str = "") -> str:
    up = "../" * depth

    def link(h, t, p):
        keep = " data-keep-filters" if p in ("leaderboard", "cost") else ""
        cur = ' aria-current="page"' if p == page else ""
        return f'<a href="{up}{h}"{keep}{cur}>{t}</a>'

    main_links = "".join(link(*x) for x in NAV)
    more_links = "".join(link(*x) for x in NAV_MORE)
    more_cur = ' class="cur"' if page in {p for _, _, p in NAV_MORE} else ""
    here = next((t for _, t, p in NAV + NAV_MORE if p == page), "Pages")
    commit = src["commit"]
    src_line = (
        f'Numbers from <a href="https://github.com/{src["repo"]}/tree/{commit}/{src["path"]}">'
        f'{src["repo"]}</a> at commit <a class="mono" href="https://github.com/{src["repo"]}/commit/{commit}">'
        f'{commit[:8]}</a> (<span class="mono">{src["path"]}</span>), copied {src["copied_at"][:10]}.'
    )
    # A short "data from" line right after the page's lede; the full source is in the footer.
    asof = (f'<p class="asof">Data from {fmt_date(src.get("commit_date") or src["copied_at"])}, '
            f'commit <span class="mono">{commit[:8]}</span> (<a href="#source-commit">source</a>)</p>')
    m = re.search(r'<p class="lede">.*?</p>', body, flags=re.S)
    if m:
        body = body[:m.end()] + "\n" + asof + body[m.end():]
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>{esc(title)}</title>
<meta name="description" content="{esc(description)}">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Source+Sans+3:wght@400;600;700&display=swap">
<link rel="stylesheet" href="{up}assets/site.css">
<link rel="icon" href="{up}assets/favicon.svg" type="image/svg+xml">
<script defer src="{up}assets/terms.js"></script>
<script defer src="{up}assets/nav.js"></script>
{head}
</head>
<body data-page="{page}">
<header class="site"><div class="wrap">
<a class="brand" href="{up}index.html">{SITE}</a>
<nav class="site" aria-label="Pages">{main_links}<details class="navmore"><summary{more_cur}>More</summary>
<div class="navpop">{more_links}</div></details></nav>
<details class="navmenu"><summary aria-label="Pages: {esc(here)}">{esc(here)}</summary>
<div class="navpop">{main_links}{more_links}</div></details>
</div></header>
<main class="wrap" data-page="{page}">
{body}
</main>
<footer class="site"><div class="wrap">
<p id="source-commit" data-commit="{commit}">{src_line}</p>
<p>Site source: <a href="https://github.com/{REPO}">{REPO}</a>. Code font: Fantasque Sans Mono 1.8.0 (SIL Open Font License).</p>
</div></footer>
</body>
</html>
"""


def legend_kinds(tools_used: set[str]) -> str:
    items = []
    for kind, k in KINDS.items():
        names = [TOOLS[t]["name"] for t in TOOLS if TOOLS[t]["kind"] == kind and t in tools_used]
        if names:
            items.append(f'<span class="item"><span class="kind k-{kind}" aria-hidden="true"></span>'
                         f'{esc(k["label"])}: {esc(", ".join(names))}</span>')
    return '<div class="legend" aria-label="Colour key">' + "".join(items) + "</div>"


def svg_marker(shape: str, color: str, filled: bool = True) -> str:
    fill = color if filled else "var(--surface)"
    if shape == "square":
        m = f'<rect x="3" y="3" width="10" height="10" fill="{fill}" stroke="{color}" stroke-width="1.6"/>'
    elif shape == "diamond":
        m = f'<path d="M8,1.5 L14.5,8 L8,14.5 L1.5,8 Z" fill="{fill}" stroke="{color}" stroke-width="1.6"/>'
    elif shape == "triangle":
        m = f'<path d="M8,2 L14,13 L2,13 Z" fill="{fill}" stroke="{color}" stroke-width="1.6"/>'
    else:
        m = f'<circle cx="8" cy="8" r="5.5" fill="{fill}" stroke="{color}" stroke-width="1.6"/>'
    return f'<svg width="16" height="16" viewBox="0 0 16 16" aria-hidden="true">{m}</svg>'


# --------------------------------------------------------------------------- data

EXPERIMENT = "qfo_region_midi_plus"
DATA_DIR = Path("data")
TRUTH_FROM_LONG = {"Pfam-A domain instances": "pfam", "Pfam-N domain instances": "pfamn",
                   "Swiss-Prot features": "swissprot"}
HARDWARE = "Stanford Sherlock, SLURM partition hns"
# Search tool -> the tool name its database-build rows carry in compute.parquet, where the
# two differ. Both MMseqs2 searches use the one database built by the `mmseqs2` step.
BUILD_TOOL = {"mmseqs2_seqseq": "mmseqs2", "mmseqs2_iterative": "mmseqs2"}


def load(data: Path):
    """Read the long tables (one row per experiment, tool, setting, target, truth, stratum and
    metric) and shape the leaderboard's experiment into the columns the pages use."""
    long = pl.read_parquet(data / "results.parquet")
    exp = long.filter(pl.col("experiment") == EXPERIMENT)
    res = (
        exp.filter(pl.col("metric").is_in(["fmax", "fmax_precision", "fmax_recall", "auprc", "median_iou_tp"]))
        .with_columns(
            pl.col("truth_source").replace_strict(TRUTH_FROM_LONG).alias("truth_set"),
            (pl.col("dedup_or_raw") == "dedup").alias("dedup"),
            pl.col("dedup_or_raw").alias("calls"),
            pl.when(pl.col("target_db") == "Pfam-A HMM library").then(pl.lit("all"))
            .otherwise(pl.col("target_db").str.replace("^QfO ", "")).alias("species"),
            pl.col("target_mya_to_human").alias("species_mya"),
            pl.col("stratum_bin").alias("stratum"),
            pl.col("n_instances").alias("n_truth_instances"),
        )
        .pivot(on="metric", index=["truth_set", "dedup", "calls", "tool", "setting", "species", "species_mya",
                                   "stratum_axis", "stratum", "n_truth_instances"], values="value")
        .sort("truth_set", "calls", "species", "stratum_axis", "stratum", "tool", "setting")
    )
    # Compute: one row per tool, setting and step, summed over the relatives searched.
    # cpu_h there is reserved core-hours (cores x wall time). The leaderboard shows the search
    # step per relative: cpu_h / n_target_dbs, per 1,000 queries.
    c = pl.read_parquet(data / "compute.parquet").filter(pl.col("experiment") == EXPERIMENT)
    search = c.filter(pl.col("setting") != "database build").filter(~pl.col("step").str.contains("Index|Db|DB"))
    build = c.filter((pl.col("setting") == "database build") | pl.col("step").str.contains("Index|Db|DB"))
    # Every build row must belong to a tool that searched; otherwise its cost would silently
    # drop off every card (the MMseqs2 build is stored as tool "mmseqs2").
    claimed = {BUILD_TOOL.get(t, t) for t in search["tool"].unique().to_list()}
    orphans = sorted(set(build["tool"].unique().to_list()) - claimed)
    if orphans:
        raise SystemExit(f"compute.parquet has database-build rows for {orphans}, which no search tool "
                         "maps to; add them to BUILD_TOOL in build_site.py")
    rows = []
    for r in search.iter_rows(named=True):
        bt = BUILD_TOOL.get(r["tool"], r["tool"])
        b = build.filter((pl.col("tool") == bt) & pl.col("setting").is_in([r["setting"], "database build"]))
        rows.append({
            "tool": r["tool"], "setting": r["setting"], "step": r["step"],
            "n_relatives": r["n_target_dbs"], "n_queries": r["n_queries"], "cores": r["cores"],
            "wall_s": r["wall_s"], "cpu_h_total": r["cpu_h"],
            "cpu_h_per_1000_queries": r["cpu_h"] / r["n_target_dbs"] * 1000 / r["n_queries"],
            "wall_h_per_1000_queries": r["wall_s"] / 3600 / r["n_target_dbs"] * 1000 / r["n_queries"],
            "peak_rss_gb": r["peak_ram_gb"],
            "gpu": "none" if r["gpu_type"] in (None, "none") else r["gpu_type"],
            "gpu_h_per_1000_queries": 0.0 if r["gpu_h"] is None else r["gpu_h"] / r["n_target_dbs"] * 1000 / r["n_queries"],
            "db_cpu_h": b["cpu_h"].sum() if b.height else None,
            "db_n": int(b["n_target_dbs"].sum()) if b.height else None,
            "db_gpu_h": b["gpu_h"].drop_nulls().sum() if b.height and b["gpu_h"].drop_nulls().len() else None,
            "db_gpu_type": next((g for g in b["gpu_type"].to_list() if g and g != "none"), None) if b.height else None,
            "hardware": HARDWARE, "cpu_h_reported": None, "n_completed_runs": 1,
        })
    comp_ts = pl.DataFrame(rows, infer_schema_length=None)
    # One row per relative this tool and setting searched, all carrying the same mean, so the
    # pages can key cost by relative; a relative it did not search gets no row ("not run").
    # Sorted so that two builds from the same data/ write byte-identical files: unique() and
    # join() do not keep row order.
    searched = res.select("tool", "setting", "species").unique().sort("tool", "setting", "species")
    comp = (searched.join(comp_ts, on=["tool", "setting"], how="inner")
            .sort("tool", "setting", "species", "step"))
    src = json.loads((data / "SOURCE.json").read_text())
    prov = json.loads((data / "provenance.json").read_text())
    prov["n_queries"] = int(comp_ts["n_queries"].max())
    cov = (data / "coverage.md").read_text()
    return res, comp, src, prov, cov


def row_order(res: pl.DataFrame) -> list[dict]:
    kinds = list(KINDS)
    pairs = res.select("tool", "setting").unique().to_dicts()
    pairs.sort(key=lambda r: (kinds.index(TOOLS[r["tool"]]["kind"]), list(TOOLS).index(r["tool"]),
                              setting_key(r["tool"], r["setting"])))
    out = []
    for i, p in enumerate(pairs):
        alpha = KMER_SETTING.fullmatch(p["setting"])
        shape = KMERSEEK_SHAPES.get(alpha[1], "circle") if p["tool"] == "kmerseek" and alpha else TOOLS[p["tool"]]["shape"]
        out.append(dict(p, order=i, label=setting_label(p["tool"], p["setting"]),
                        words=setting_words(p["tool"], p["setting"]), shape=shape))
    return out


def results_json(res: pl.DataFrame) -> dict:
    out: dict[str, dict] = {}
    for r in res.iter_rows(named=True):
        key = "|".join([r["truth_set"], r["calls"], r["species"], r["stratum_axis"], r["stratum"]])
        k2 = f'{r["tool"]}|{r["setting"]}'
        if k2 in out.setdefault(key, {}):
            raise SystemExit(f"two rows for {key} {k2}")
        out[key][k2] = [r["fmax"], r["auprc"], r["median_iou_tp"], r["n_truth_instances"]]
    return out


def compute_json(comp: pl.DataFrame) -> dict:
    out = {}
    for r in comp.iter_rows(named=True):
        rep = r["cpu_h_reported"]
        out[f'{r["tool"]}|{r["setting"]}|{r["species"]}'] = {
            "cpu_h_per_1000": r["cpu_h_per_1000_queries"],
            "cpu_h_reported_per_1000": None if rep is None else rep * 1000 / r["n_queries"],
            "peak_rss_gb": r["peak_rss_gb"],
            "gpu": r["gpu"],
            "gpu_h_per_1000": r["gpu_h_per_1000_queries"],
            "hardware": r["hardware"],
            "n_relatives": r["n_relatives"],
            "cores": r["cores"],
            "wall_h_per_1000": r["wall_h_per_1000_queries"],
            # Mean core-hours to build one relative's database, over the builds the run's task
            # records hold (a build reused from an earlier run left no record).
            "db_cpu_h_per_db": None if r["db_cpu_h"] is None else r["db_cpu_h"] / r["db_n"],
            "db_n": r["db_n"],
        }
    return out


# --------------------------------------------------------------------------- pages

def key_line(tools_used: set[str]) -> str:
    """One line above Table 1: the kind colours and the marks inside the table."""
    kinds = "".join(f'<span class="item"><span class="kind k-{k}" aria-hidden="true"></span>{esc(SHORT_KIND[k])}</span>'
                    for k in KINDS if any(TOOLS[t]["kind"] == k for t in tools_used))
    return (f'<div class="legend keyline" aria-label="Key to Table 1">{kinds}'
            '<span class="item"><span class="nr" aria-hidden="true">×</span> not run</span>'
            '<span class="item"><span class="tag warn">graded on its own models</span></span></div>')


def page_leaderboard(meta_json: str, src: dict, n_queries: int, tools_used: set[str], dd: dict) -> str:
    body = f"""
<h1>{SITE}</h1>
<p class="lede">Which tool puts the known domains and features of human proteins in the right
place, and what does it cost?</p>
<p class="verdict" id="verdict" aria-live="polite"></p>
<details class="choices" id="choices">
<summary><span class="choices-label">Change what is ranked</span>
<span class="choices-now" id="filtersummary"></span></summary>
<form class="filters" id="filters-form" onsubmit="return false" aria-label="Filters">
<div id="filters"></div>
<p class="axisnote">The benchmark split the proteins along one of the last four rows at a time,
so picking a value in one of them sets the other three back to "any".</p>
</form>
<div class="prices" aria-label="Prices used for the cost column">
<span class="pricehead">Prices for the estimated cost {term_btn("cost", "prices")}</span>
<label>r_cpu, USD per CPU-hour <input id="r_cpu" type="number" min="0" step="0.001" inputmode="decimal"></label>
<label>r_gpu, USD per GPU-hour <input id="r_gpu" type="number" min="0" step="0.01" inputmode="decimal"></label>
<span class="muted">Default prices are AWS on-demand rates for US East (N. Virginia), Linux, read
2026-10-15 from the <a href="{PRICE_PAGE}">EC2 on-demand price page</a> (<a href="{PRICE_FILE}">its
price file</a>). CPU: a c7i.large costs $0.08925 per hour for 2 vCPU, so $0.044625 per CPU-hour.
GPU: a g5.xlarge with one NVIDIA A10G costs $1.006 per hour.</span>
</div>
</details>
<figure class="lbfig">
<figcaption><b>Table 1.</b> <span id="rankline">Tools ranked by Fmax, highest first, best setting of each tool.</span>
<span id="count"></span>
<span class="tablebtns"><button type="button" class="chip" id="viewbtn">Show every setting</button>
<button type="button" class="chip" id="colbtn" aria-pressed="false">More columns</button></span></figcaption>
{key_line(tools_used)}
<div class="tablebox"><table class="lb" id="lb"><thead><tr></tr></thead><tbody></tbody></table></div>
</figure>
<details class="key"><summary>The marks and colours in Table 1, in full</summary>
{legend_kinds(tools_used)}
<div class="legend">
<span class="item"><span class="nr" aria-hidden="true">×</span> not run: no row in the results table for these filters</span>
<span class="item"><i>not measured</i>: the search ran, but no task record holds its cost</span>
<span class="item"><i>no correct calls</i>: no call was correct, so there is no median IoU</span>
<span class="item"><i>correct calls only on single residues</i>: Fmax is above 0, but every correct
call is on a one-residue annotation, where an IoU says nothing about where the ends are</span>
<span class="item"><span class="tag">needs 3D structure prediction</span> <span class="tag">needs GPU</span>
what the tool needs besides the sequence</span>
<span class="item"><span class="tag warn">graded on its own models</span> hmmscan under the Pfam truth
set: Pfam's domains are defined by the models hmmscan runs, so read its row as an upper limit</span>
<span class="item">Each row shows a tool's best setting for the column it is sorted by;
"+N more settings" opens the rest. Click a column name to sort; rows with no value sort last.</span>
<span class="item">Cost and core-hours are measured on Stanford's Sherlock cluster, from the run's own task records {term_btn("measured", "key")}</span>
</div>
</details>
<p class="prose oneline">In plain words: <b>Fmax</b> asks whether a tool finds the real domains
without adding wrong ones. <b>AUPRC</b> asks whether its correct calls stay at the top of its list.
<b>Median IoU</b> asks how well the ends of its right calls line up with the real domain. All run
0 to 1, and higher is better. Hover or tap a <b>?</b> for more, or see
<a href="explainers.html">How the scores work</a>, which works each through on a made-up protein.</p>
<p class="prose oneline">Each tool searched {n_queries:,} human proteins, most from chromosome 6, against one
relative's proteome and copied its annotations (<a href="methods.html#data">which proteins</a>).</p>
<p class="prose oneline">Estimated cost = core-hours × r_cpu + GPU-hours × r_gpu, per 1,000 queries,
for the search alone (<a href="methods.html#cost">how cost is counted</a>).</p>
<details class="prose"><summary>How to choose between Fmax and AUPRC, and raw against dedup</summary>
<p>Fmax is the score at the single best cutoff on the tool's own score, protein by protein: use it
when you will keep every call above one cutoff. AUPRC sums precision over every cutoff, counting
each call: use it when you will rank calls and read down the list. The two can rank tools in
opposite orders, because AUPRC pools every call of every protein while Fmax averages protein by
protein. A domain reported many times counts against both. Dedup merges those copies; for phmmer against zebrafish (Pfam) it moves Fmax from
{dd['raw_fmax']:.3f} to {dd['dedup_fmax']:.3f} and AUPRC from {dd['raw_auprc']:.3f} to
{dd['dedup_auprc']:.3f}. Definitions: <a href="methods.html#dedup">Methods</a>.</p></details>
<script type="application/json" id="lb-data">{meta_json}</script>
<script src="assets/app.js"></script>
"""
    return layout(SITE, "leaderboard", body, src,
                  description="Protein annotation tools ranked by Fmax, AUPRC and boundary accuracy, with measured cost.")


def page_cost(meta_json: str, src: dict, n_queries: int, tools_used: set[str]) -> str:
    legend = []
    for kind, k in KINDS.items():
        for t, info in TOOLS.items():
            if info["kind"] != kind or t not in tools_used:
                continue
            if t == "kmerseek":
                # One entry: kmerseek's colour, with the shape of each alphabet, smallest first.
                shapes = "".join(f'<span class="sub">{svg_marker(KMERSEEK_SHAPES[a], f"var(--k-{kind})")}{esc(a)}</span>'
                                  for a in ALPHA_ORDER)
                legend.append(f'<span class="item kwrap"><span class="nowrap">kmerseek, by alphabet:</span>{shapes}</span>')
            else:
                legend.append(f'<span class="item">{svg_marker(info["shape"], f"var(--k-{kind})")}{esc(info["name"])}</span>')
    body = f"""
<h1>Compare two columns</h1>
<p class="lede">Any two columns of the leaderboard, drawn against each other: a score against
what the search cost, or one score against another. Pick the column for each axis below; the ?
beside each gives its definition.</p>
<div class="axispick" id="axispick" aria-label="Columns on the axes"></div>
<details class="choices" id="choices">
<summary><span class="choices-label">Change what is plotted</span>
<span class="choices-now" id="filtersummary"></span></summary>
<form class="filters" id="filters-form" onsubmit="return false" aria-label="Filters">
<div id="filters"></div>
<p class="axisnote">The benchmark split the proteins along one of the last four rows at a time,
so picking a value in one of them sets the other three back to "any".</p>
</form>
<div class="prices" aria-label="Prices used for the estimated cost">
<label>r_cpu, USD per CPU-hour <input id="r_cpu" type="number" min="0" step="0.001" inputmode="decimal"></label>
<label>r_gpu, USD per GPU-hour <input id="r_gpu" type="number" min="0" step="0.01" inputmode="decimal"></label>
<span class="muted">Used only when the horizontal axis is the estimated cost. Default prices are AWS
on-demand rates for US East (N. Virginia), Linux, read 2026-10-15 from the
<a href="{PRICE_PAGE}">EC2 on-demand price page</a>: $0.044625 per CPU-hour (a c7i.large, $0.08925 per hour for 2 vCPU) and
$1.006 per GPU-hour (a g5.xlarge with one NVIDIA A10G).</span>
</div>
</details>
<figure>
<figcaption><b>Figure 3.</b> <span id="plotcaption"></span> One marker per tool and setting, for
the filters above. Hover or focus a marker to see its row; click it to open the tool's card.</figcaption>
<div class="legend" aria-label="Key to Figure 3">
{"".join(legend)}
<span class="item" id="frontier-key"><svg width="28" height="10" aria-hidden="true"><line x1="1" y1="5" x2="27" y2="5" stroke="var(--ink)" stroke-width="1.6"/></svg>
<span id="frontier-words">no other tool is both cheaper and better</span></span>
</div>
<div class="plotbox" id="plot"></div>
</figure>
<p id="notplotted" class="prose muted"></p>
<p class="prose muted">Colour is the kind of method; within a kind, marker shape is the tool, and for
kmerseek the alphabet. Every cost here was measured on Sherlock; none was reported by a submitter. Core-hours, cost and wall-clock time are per 1,000 queries for searching
one relative, the mean over the relatives each tool searched, measured on Stanford's Sherlock
cluster. Peak RAM is the largest memory any one search task used. Wall-clock
time is on the cores each tool was given (between 4 and 16; hover a marker for the number), so it
can rank tools differently from core-hours. The database build is the mean over the builds the run's task records
hold; a database reused from an earlier run left no record. Definitions:
<a href="methods.html#cost">Methods</a>.</p>
<script type="application/json" id="lb-data">{meta_json}</script>
<script src="assets/app.js"></script>
"""
    return layout("Compare two columns", "cost", body, src,
                  description="Any two leaderboard columns against each other: Fmax, AUPRC or median IoU "
                              "against measured core-hours, cost, wall-clock time, peak RAM or database build time.")


def fnum(x, d=3):
    return f"{x:.{d}f}"


def best_by(res: pl.DataFrame, **filt) -> pl.DataFrame:
    q = res
    for k, v in filt.items():
        q = q.filter(pl.col(k) == v)
    return q


def strong_weak(comparisons: list) -> tuple[list, list]:
    """From comparisons sorted best margin first: the two strongest, and the two weakest
    (weakest first), never repeating a stratum already listed as strong."""
    return comparisons[:2], list(reversed(comparisons[max(2, len(comparisons) - 2):]))


def card_numbers(tool: str, res: pl.DataFrame, comp: pl.DataFrame) -> dict:
    """Strong and weak subsets, cost and the three questions, all read from the tables."""
    species = "all" if tool == "hmmscan" else "zebrafish"
    base = dict(truth_set="pfam", calls="raw", species=species)
    mine_all = best_by(res, tool=tool, stratum_axis="all", **base)
    out: dict = {"species": species, "rows": []}
    if mine_all.height == 0:
        return out
    best = mine_all.sort(["fmax", "setting"], descending=[True, False]).row(0, named=True)
    out["best_setting"] = best["setting"]
    out["all_settings"] = res.filter(pl.col("tool") == tool).select("setting").unique().to_dicts()
    out["best_all"] = best
    others = res.filter(pl.col("tool") != tool)

    comparisons = []
    for axis in ("identity", "plddt", "feature_length_bin"):
        mine = best_by(res, tool=tool, setting=best["setting"], stratum_axis=axis, **base)
        for r in mine.filter(pl.col("n_truth_instances") >= MIN_N).iter_rows(named=True):
            o = best_by(others, stratum_axis=axis, stratum=r["stratum"], **base)
            if o.height == 0:
                continue
            ob = o.sort(["fmax", "tool", "setting"], descending=[True, False, False]).row(0, named=True)
            comparisons.append({
                "axis": axis, "stratum": r["stratum"], "fmax": r["fmax"],
                "n": r["n_truth_instances"], "other_tool": ob["tool"],
                "other_setting": ob["setting"], "other_fmax": ob["fmax"],
                "margin": r["fmax"] - ob["fmax"],
            })
    comparisons.sort(key=lambda c: (-c["margin"], c["axis"], c["stratum"]))
    out["strong"], out["weak"] = strong_weak(comparisons)

    # Swiss-Prot feature types, same relative, this tool's best setting there.
    sp = best_by(res, tool=tool, truth_set="swissprot", calls="raw", species=species,
                 stratum_axis="feature_type")
    if sp.height:
        s_all = best_by(res, tool=tool, truth_set="swissprot", calls="raw", species=species,
                        stratum_axis="all").sort(["fmax", "setting"], descending=[True, False]).row(0, named=True)
        sp = sp.filter((pl.col("setting") == s_all["setting"]) & (pl.col("n_truth_instances") >= MIN_N)) \
            .sort(["fmax", "stratum"], descending=[True, False])
    if sp.height:
        out["sp_setting"] = s_all["setting"]
        out["sp_best"] = sp.row(0, named=True)
        out["sp_worst"] = sp.row(sp.height - 1, named=True)

    c = comp.filter((pl.col("tool") == tool) & (pl.col("setting") == best["setting"]))
    here = c.filter(pl.col("species") == species)
    out["cost_here"] = here.row(0, named=True) if here.height else None


    ident = {r["stratum"]: r for r in best_by(res, tool=tool, setting=best["setting"],
                                              stratum_axis="identity", **base).iter_rows(named=True)}
    out["identity"] = ident
    length = {r["stratum"]: r for r in best_by(res, tool=tool, setting=best["setting"],
                                               stratum_axis="feature_length_bin", **base).iter_rows(named=True)}
    out["length"] = length
    return out


def stratum_words(axis, stratum):
    v = STRATUM_WORDS[axis].get(stratum, stratum)
    if axis == "identity":
        return "no same-family domain in the relative" if stratum == "no_homolog" else f"identity {v} to the closest same-family domain"
    if axis == "feature_length_bin":
        return f"regions {v} long"
    if axis == "plddt":
        return f"AlphaFold confidence (pLDDT) {v}"
    return v


def card_markdown(tool: str, nums: dict, src: dict) -> tuple[str, str]:
    """The card as Markdown, in two parts: what a reader choosing a tool needs first (the
    one-line description, "Will it work for you?", best Fmax and cost, strong and weak
    subsets), and the rest. The HTML card puts its figure and the claims between them."""
    t = TOOLS[tool]
    sp_name = SPECIES[nums["species"]]
    where = ("Pfam truth set, raw calls, no relative (hmmscan searches the Pfam models)"
             if nums["species"] == "all" else DEFAULT_WORDS)
    first = t["compares"].split(". ")[0].rstrip(".") + "."
    T = [f"# {t['name']}", "", first, "", "## Will it work for you?", ""]
    T += questions(tool, nums)
    T += ["", "## Best Fmax and cost", ""]
    ch = nums.get("cost_here")
    if "best_setting" not in nums:
        T += ["Not run under the Pfam truth set for this relative.", ""]
    else:
        b = nums["best_all"]
        bs = nums["best_setting"]
        line = (f"Under the {where}, its best setting is {setting_label(tool, bs)} (`{bs}`): Fmax "
                f"{fnum(b['fmax'])} over {b['n_truth_instances']:,} true instances.")
        if ch:
            cost = ch["cpu_h_per_1000_queries"] * R_CPU + ch["gpu_h_per_1000_queries"] * R_GPU
            line += (f" Searching costs about ${cost:.3g} per 1,000 queries "
                     f"({ch['cpu_h_per_1000_queries']:.3g} core-hours), "
                     + ("for scanning Pfam's models." if nums["species"] == "all" else "for each relative searched."))
        else:
            line += " Its search cost was not measured."
        T += [line, ""]
    T += ["## Where it is strong and weak", ""]
    if "best_setting" not in nums:
        T += ["Not run under the Pfam truth set for this relative.", ""]
    else:
        T += ["Each line compares this setting with the best setting of any other tool on the same "
              "subset of the proteins, so the other tool's setting can change from line to line. "
              f"Only subsets with at least {MIN_N} true instances are compared.", ""]

        def lines(cs):
            return [f"- {stratum_words(c['axis'], c['stratum'])}: Fmax {fnum(c['fmax'])} "
                    f"(n = {c['n']:,}); best other: {TOOLS[c['other_tool']]['name']}, "
                    f"{setting_label(c['other_tool'], c['other_setting'])}, {fnum(c['other_fmax'])}."
                    for c in cs]
        if nums["strong"]:
            T += ["Where it is closest to, or furthest ahead of, the best other tool:", ""] + lines(nums["strong"]) + [""]
        if nums["weak"]:
            T += ["Where it is furthest behind the best other tool:", ""] + lines(nums["weak"]) + [""]
        if not nums["strong"] and not nums["weak"]:
            T += ["No other tool was scored on the same subsets, so no comparison: not run.", ""]
        if "sp_best" in nums:
            sb, sw = nums["sp_best"], nums["sp_worst"]
            T += [f"Swiss-Prot features, raw calls, {sp_name}, setting {setting_label(tool, nums['sp_setting'])}: "
                  f"best feature type {FEATURE_WORDS.get(sb['stratum'], sb['stratum'])}, Fmax "
                  f"{fnum(sb['fmax'])} (n = {sb['n_truth_instances']:,}); worst "
                  f"{FEATURE_WORDS.get(sw['stratum'], sw['stratum'])}, Fmax {fnum(sw['fmax'])} "
                  f"(n = {sw['n_truth_instances']:,}).", ""]
        else:
            T += ["Swiss-Prot feature types: not run.", ""]

    L = ["## What it compares", "", t["compares"], "", "## What it needs", ""]
    for k, v in t["needs"].items():
        L.append(f"- {k}: {'yes' if v else 'no'}")
    L += ["", "## Reference searched", "", t["reference"], "",
          "## How to read its score", "", t["score"], ""]
    L += ['## Known failure modes {#failure-modes}', ""]
    L += [f"- {f}" for f in t["failure_modes"]]
    L.append("")

    L += ["## Cost per 1,000 queries", ""]
    if ch:
        gpu = ch["gpu_h_per_1000_queries"]
        cost = ch["cpu_h_per_1000_queries"] * R_CPU + gpu * R_GPU
        target = ("scanning Pfam's models" if nums["species"] == "all"
                  else f"searching one relative (mean over the {ch['n_relatives']} relatives it searched)")
        L += [f"Setting `{nums['best_setting']}`, {target}: {ch['cpu_h_per_1000_queries']:.3g} "
              f"core-hours reserved per 1,000 queries ({ch['cores']} cores for the task's wall time), "
              f"measured on Sherlock. Peak memory {ch['peak_rss_gb']:.1f} GB. That is about ${cost:.3g} "
              f"at ${R_CPU:.6g} per CPU-hour."]
        if ch["db_cpu_h"] is not None:
            gtxt = (f", plus {ch['db_gpu_h']:.3g} GPU-hours ({ch['db_gpu_type']})" if ch["db_gpu_h"] else "")
            L += [f"Building its databases took {ch['db_cpu_h']:.3g} core-hours over {ch['db_n']} "
                  f"database{'s' if ch['db_n'] != 1 else ''}{gtxt}. A database is built once and reused, "
                  "so this is not in the per-query cost."]
            if ch["db_n"] != ch["n_relatives"] and nums["species"] != "all":
                L += [f"This covers the {ch['db_n']} database builds the task records hold; the search "
                      f"used {ch['n_relatives']} relatives."]
        else:
            L += [f"Building its databases for setting `{nums['best_setting']}`: not measured (no task "
                  "record holds a database build for this tool and setting)."]
        L += [""]
    elif "best_setting" in nums:
        L += [f"Setting `{nums['best_setting']}`: not measured (no task record holds its search).", ""]
    else:
        L += ["Not run.", ""]
    L += ["The per-query cost counts the search only. The run recorded cost per tool and setting, "
          "summed over all relatives, so the cost for one relative is the mean.", ""]

    L += ["## Licence", "",
          f"- Licence: {t['licence']} ([source]({t['licence_url']}))",
          f"- Commercial use allowed: {t['commercial']}", ""]

    settings = sorted({r["setting"] for r in nums.get("all_settings", [])}, key=lambda st: setting_key(tool, st))
    if settings:
        L += ["## Settings in the table", ""]
        L += [f"- {setting_label(tool, st)} (`{st}`): {setting_words(tool, st)}" for st in settings]
        L += [""]
    L += ["## Version, container, command", "",
          f"- Version: {t['version']}", f"- Container: `{t['container']}`", "",
          "```bash", t["command"], "```", "", t["command_note"], "",
          f"Tool id `{tool}`. Numbers from {src['repo']} commit `{src['commit'][:8]}`.", ""]
    return "\n".join(T) + "\n", "\n".join(L) + "\n"


def card_figure(tool: str, nums: dict, res: pl.DataFrame) -> str:
    """Figure 2 of a card: Fmax by identity bin and by region length bin, this tool's best
    setting against the best other tool in each bin, on one shared 0-to-max axis."""
    if "best_setting" not in nums:
        return ""
    base = dict(truth_set="pfam", calls="raw", species=nums["species"])
    others = res.filter(pl.col("tool") != tool)
    groups = [("identity", "% identity to the closest same-family domain", [v for v, _ in IDENTITY]),
              ("feature_length_bin", "length of the annotated region",
               [v for v, _ in LENGTH_SHORT + LENGTH_DOMAIN] + ["whole_protein"])]
    rows, few = [], []
    for axis, title, order in groups:
        mine = {r["stratum"]: r for r in best_by(res, tool=tool, setting=nums["best_setting"],
                                                stratum_axis=axis, **base).iter_rows(named=True)}
        g = []
        for st in order:
            r = mine.get(st)
            if r is None:
                continue
            if r["n_truth_instances"] < MIN_N:
                few.append(STRATUM_WORDS[axis].get(st, st))
                continue
            o = best_by(others, stratum_axis=axis, stratum=st, **base)
            ob = o.sort(["fmax", "tool", "setting"], descending=[True, False, False]).row(0, named=True) if o.height else None
            g.append((st, r, ob))
        if g:
            rows.append((title, axis, g))
    if not rows:
        return ""
    vals = [r["fmax"] for _, _, g in rows for _, r, _ in g] + [ob["fmax"] for _, _, g in rows for _, _, ob in g if ob]
    xmax = max(0.1, math.ceil(max(vals) * 10) / 10)
    # Drawn at about phone width, so its 11-12 px text stays that size on a phone.
    W, L, R, top = 360, 104, 122, 8
    px = lambda v: L + v / xmax * (W - L - R)  # noqa: E731
    color = f"var(--k-{TOOLS[tool]['kind']})"
    parts, y = [], top
    for title, axis, g in rows:
        y += 14
        parts.append(f'<text x="0" y="{y}" class="cf-h">{esc(title)}</text>')
        y += 2
        for st, r, ob in g:
            y += 20
            lab = "none in relative" if st == "no_homolog" else STRATUM_WORDS[axis].get(st, st)
            parts.append(f'<line x1="{L}" x2="{W - R}" y1="{y}" y2="{y}" class="cf-row"/>')
            parts.append(f'<text x="{L - 6}" y="{y + 4}" text-anchor="end" class="cf-l">{esc(lab)}</text>')
            if ob:
                oname = TOOLS[ob["tool"]]["name"]
                parts.append(f'<path d="M{px(ob["fmax"]):.1f},{y - 5} l5,5 l-5,5 l-5,-5 z" class="cf-o">'
                             f'<title>best other tool: {esc(oname)}, {esc(setting_label(ob["tool"], ob["setting"]))}, Fmax {ob["fmax"]:.3f}</title></path>')
                parts.append(f'<text x="{W - R + 6}" y="{y + 4}" class="cf-n">{esc(oname)}</text>')
            parts.append(f'<circle cx="{px(r["fmax"]):.1f}" cy="{y}" r="5" fill="{color}" stroke="var(--surface)" stroke-width="1">'
                         f'<title>{esc(TOOLS[tool]["name"])}: Fmax {r["fmax"]:.3f}, n = {r["n_truth_instances"]:,}</title></circle>')
        y += 14
    axis_y = y
    step = 0.1 if xmax <= 0.6 else 0.2
    ticks = [round(k * step, 1) for k in range(int(round(xmax / step)) + 1)]
    for tk in ticks:
        parts.append(f'<line x1="{px(tk):.1f}" x2="{px(tk):.1f}" y1="{top + 20}" y2="{axis_y}" class="cf-grid"/>')
        parts.append(f'<text x="{px(tk):.1f}" y="{axis_y + 14}" text-anchor="middle" class="cf-t">{tk:g}</text>')
    parts.append(f'<text x="{(L + W - R) / 2}" y="{axis_y + 30}" text-anchor="middle" class="cf-l">Fmax (higher is better)</text>')
    H = axis_y + 38
    other_mark = '<svg width="14" height="14" aria-hidden="true"><path d="M7,2 l5,5 l-5,5 l-5,-5 z" class="cf-o"/></svg>'
    mine_mark = f'<svg width="14" height="14" aria-hidden="true"><circle cx="7" cy="7" r="5" fill="{color}"/></svg>'
    note = (f' Left out, under {MIN_N} true instances: {esc(", ".join(few))}.' if few else "")
    where = ("Pfam truth set, raw calls, no relative" if nums["species"] == "all"
             else "Pfam truth set, raw calls, zebrafish as the relative")
    return f"""<figure class="cardfig">
<figcaption><b>Figure 2.</b> Where {esc(TOOLS[tool]['name'])} is ahead and behind: its Fmax on each subset of the
proteins, against the best other tool on that subset ({where}; its best setting,
{esc(setting_label(tool, nums['best_setting']))}).{note}</figcaption>
<div class="legend"><span class="item">{mine_mark}{esc(TOOLS[tool]['name'])}</span>
<span class="item">{other_mark}best other tool on that subset, named at the right</span></div>
<svg class="cf" viewBox="0 0 {W} {H}" role="img" aria-label="Fmax by subset for {esc(TOOLS[tool]['name'])} and the best other tool">{"".join(parts)}</svg>
</figure>"""


def questions(tool: str, nums: dict) -> list[str]:
    """Three yes/no questions. Each answer is the Fmax range the table shows for it."""
    t = TOOLS[tool]
    q = []

    def answer(rows: dict, yes: list[str], no: list[str], words: dict) -> str | None:
        ok = {k: r for k, r in rows.items() if r["n_truth_instances"] >= MIN_N}
        y = [k for k in yes if k in ok]
        n = [k for k in no if k in ok]
        if not y or not n:
            return None
        show = lambda ks: ", ".join(  # noqa: E731
            f"{words[k]} {fnum(ok[k]['fmax'])} (n = {ok[k]['n_truth_instances']:,})" for k in ks)
        return f"If yes, Fmax here is {show(y)}. If no, {show(n)}."

    base_words = ("Pfam truth set, raw calls, no relative searched" if nums["species"] == "all"
                  else DEFAULT_WORDS)
    bs = nums.get("best_setting")
    where = f"{base_words}, setting {setting_label(tool, bs) if bs else 'none'}"
    a1 = answer(nums.get("identity", {}), ["40-60%", "60-100%"], ["0-20%", "20-30%", "30-40%"],
                STRATUM_WORDS["identity"])
    if a1:
        q.append("1. **Is the closest same-family domain in the relative you will search at least "
                 f"40% identical, for most of your proteins?** {a1} ({where}; subsets with at least "
                 f"{MIN_N} true instances.)")
    elif nums.get("identity"):
        q.append("1. **Is the closest same-family domain in the relative you will search at least "
                 f"40% identical?** Too few true instances (under {MIN_N}) in the identity subsets to say.")
    else:
        q.append("1. **Is every family you expect to find already in Pfam?** If no, this tool "
                 "cannot find it: it only reports families that have a Pfam model.")
    if t["requirement_question"]:
        a, b = t["requirement_question"]
        q.append(f"2. **{a}** {b}")
    else:
        a2 = answer(nums.get("length", {}), ["61-120", "121-250", "251+"],
                    ["2-15", "16-30", "31-60"], STRATUM_WORDS["feature_length_bin"])
        q.append("2. **Are the features you care about longer than 60 residues?** "
                 + (f"{a2} ({where}.)" if a2 else "Fmax by region length for this tool: not run."))
    ch = nums.get("cost_here")
    if ch:
        gpu = ch["gpu_h_per_1000_queries"]
        extra = f" and {gpu:.3g} GPU-hours" if gpu else ""
        if nums["species"] == "all":
            q.append(f"3. **Can you spend about {ch['cpu_h_per_1000_queries']:.3g} core-hours per "
                     "1,000 proteins?** That is the measured cost of scanning them against Pfam's "
                     "models.")
        else:
            q.append(f"3. **Can you spend about {ch['cpu_h_per_1000_queries']:.3g} core-hours{extra} "
                     "per 1,000 proteins for each relative you search?** That is the measured mean "
                     "per relative; see \"Cost per 1,000 queries\" below.")
    else:
        q.append("3. **Can you afford the search?** Its cost was not measured in this run.")
    return q


def page_card(tool: str, md: tuple[str, str], src: dict, fig: str, cl: dict | None = None) -> str:
    ext = ["fenced_code", "attr_list", "tables"]
    top, rest = (markdown.markdown(x, extensions=ext) for x in md)
    crumb = f'<nav class="crumb" aria-label="Breadcrumb"><a href="index.html">Tool cards</a> › {esc(TOOLS[tool]["name"])}</nav>'
    body = (f'{crumb}<div class="card">{top}{fig}</div>\n'
            + (claims_card(cl, up="../") if cl else "")
            + f'<div class="card card-rest">{rest}</div>\n<p><a href="index.html">All tool cards</a> · '
            f'<a href="{tool}.md">this card as Markdown</a> · <a href="{tool}.json">as JSON</a></p>')
    return layout(f"{TOOLS[tool]['name']} card", "tools", body, src, depth=1,
                  description=f"What {TOOLS[tool]['name']} compares, what it needs, where it is strong and weak, and what it costs.")


# Licence categories a tool card can be filtered by; each tool's is set in scripts/tools.py.
LICENCE_TAGS = ("permissive licence", "GPL-3.0")
SHORT_KIND = {"letters": "aligns letters", "profiles": "family profiles", "shapes": "3-D shape letters",
              "kmers": "reduced-alphabet k-mers", "models": "Pfam models"}


def tool_filters(t: str) -> dict:
    """The three things a reader filters tool cards by, from the tool facts."""
    info = TOOLS[t]
    if info["licence_tag"] not in LICENCE_TAGS:
        raise SystemExit(f"{t}: licence_tag {info['licence_tag']!r} is not one of {LICENCE_TAGS}")
    return {"structure": "needed" if info["needs"]["3D structure prediction"] else "not needed",
            "gpu": "needed" if info["needs"]["GPU"] else "not needed",
            "licence": info["licence_tag"]}


def best_by_relative(res: pl.DataFrame) -> dict:
    """tool -> {relative: best Fmax over its settings}, Pfam domains, raw calls, all proteins."""
    q = res.filter((pl.col("truth_set") == "pfam") & (pl.col("calls") == "raw") & (pl.col("stratum_axis") == "all"))
    out: dict = {}
    for r in q.group_by("tool", "species").agg(pl.col("fmax").max()).iter_rows(named=True):
        out.setdefault(r["tool"], {})[r["species"]] = r["fmax"]
    return out


def fmax_strip(t: str, best: dict, rel_order: list[str], ymax: float) -> str:
    """One bar per relative, nearest to farthest, height = the tool's best Fmax there on the
    0-to-ymax scale every card shares; a grey x where it was not run."""
    mine = best.get(t, {})
    if set(mine) <= {"all"}:
        return ('<p class="strip-none">Searches Pfam\'s models, not a relative, so it has one number '
                'and no bar per relative.</p>')
    W, H, L, B = 300, 64, 30, 30
    bw = (W - L) / len(rel_order)
    sy = lambda v: 4 + (1 - v / ymax) * (H - B - 4)  # noqa: E731
    parts = [f'<line x1="{L}" x2="{W}" y1="{sy(0):.1f}" y2="{sy(0):.1f}" class="st-base"/>',
             f'<line x1="{L}" x2="{W}" y1="{sy(ymax):.1f}" y2="{sy(ymax):.1f}" class="st-grid"/>',
             f'<text x="{L - 4}" y="{sy(ymax) + 4:.1f}" text-anchor="end" class="st-t">{ymax:.2f}</text>',
             f'<text x="{L - 4}" y="{sy(0) + 4:.1f}" text-anchor="end" class="st-t">0</text>']
    labels = []
    for k, sp in enumerate(rel_order):
        x = L + k * bw
        v = mine.get(sp)
        if v is None:
            cx, cy = x + bw / 2, sy(0) - 6
            parts.append(f'<path d="M{cx - 3.5:.1f},{cy - 3.5:.1f} l7,7 M{cx + 3.5:.1f},{cy - 3.5:.1f} l-7,7" class="st-nr"/>')
            labels.append(f"{SPECIES[sp]} ({rel_age(sp)}): not run")
        else:
            parts.append(f'<rect x="{x + 1.5:.1f}" y="{sy(v):.1f}" width="{bw - 3:.1f}" height="{sy(0) - sy(v):.1f}" '
                         f'fill="var(--k-{TOOLS[t]["kind"]})"><title>{esc(SPECIES[sp])}, {rel_age(sp)}: Fmax {v:.3f}</title></rect>')
            labels.append(f"{SPECIES[sp]} ({rel_age(sp)}): {v:.3f}")
    f, l = rel_order[0], rel_order[-1]
    parts.append(f'<text x="{L}" y="{H - 12}" class="st-l">{esc(SPECIES[f])}, {rel_age(f)}</text>')
    parts.append(f'<text x="{W}" y="{H - 12}" text-anchor="end" class="st-l">{esc(SPECIES[l])}, {rel_age(l)}</text>')
    return (f'<svg class="strip" viewBox="0 0 {W} {H}" role="img" aria-label="best Fmax at each relative: {esc("; ".join(labels))}">'
            f'{"".join(parts)}</svg>')


FILTER_GROUPS = [("structure", "3D structure prediction"), ("gpu", "GPU"), ("licence", "Licence")]


def page_tools_index(src: dict, tools_used: set[str], res: pl.DataFrame, comp: pl.DataFrame,
                     claims_by_tool: dict) -> str:
    rel_order = (res.filter(pl.col("species") != "all").select("species", "species_mya").unique()
                 .sort("species_mya", "species", nulls_last=True)["species"].to_list())
    best = best_by_relative(res)
    ymax = math.ceil(max(v for t in best.values() for sp, v in t.items() if sp != "all") * 20) / 20
    entries = []
    for t in tools_used:
        nums = card_numbers(t, res, comp)
        b = nums.get("best_all")
        ch = nums.get("cost_here")
        cost = (ch["cpu_h_per_1000_queries"] * R_CPU + ch["gpu_h_per_1000_queries"] * R_GPU) if ch else None
        entries.append((t, nums, b["fmax"] if b else None, cost))
    # Best Fmax first; hmmscan, which searches no relative, last.
    entries.sort(key=lambda e: (e[0] == "hmmscan", -(e[2] or 0), e[0]))
    cards = []
    for t, nums, fmax, cost in entries:
        info = TOOLS[t]
        f = tool_filters(t)
        cl = claims_by_tool.get(t)
        tags = [f'<span class="tc-tag">needs 3D structure prediction</span>' if f["structure"] == "needed" else "",
                f'<span class="tc-tag">needs GPU</span>' if f["gpu"] == "needed" else "",
                f'<span class="tc-tag">{esc(f["licence"])}</span>',
                '<span class="tag warn">graded on its own models under Pfam</span>' if t == "hmmscan" else "",
                f'<span class="tc-tag">claims: {cl["counts"]["alive"]} hold{"s" if cl["counts"]["alive"] == 1 else ""}, {cl["counts"]["retired"]} disproved</span>' if cl else ""]
        where = "no relative" if nums["species"] == "all" else "zebrafish"
        nums_html = (f'<div class="tc-nums"><span><b>{fmax:.3f}</b> best Fmax</span>'
                     f'<span><b>{"$" + format(cost, ".3g") if cost is not None else "not measured"}</b> per 1,000 queries</span>'
                     f'<span class="muted">Pfam, raw calls, {where}; {esc(setting_label(t, nums["best_setting"]))}</span></div>'
                     if fmax is not None else '<div class="tc-nums"><span><span class="nr">×</span> not run</span></div>')
        first = info["compares"].split(". ")[0].rstrip(".") + "."
        data = " ".join(f'data-{k}="{esc(v)}"' for k, v in f.items())
        cards.append(f"""<li class="tc" {data}>
<h2 class="tc-name"><span class="kind k-{info['kind']}" aria-hidden="true"></span><a href="{t}.html">{esc(info["name"])}</a></h2>
{nums_html}
<p class="tc-desc">{esc(first)}</p>
{fmax_strip(t, best, rel_order, ymax)}
<div class="tc-tags">{"".join(x for x in tags if x)}</div></li>""")
    rows = []
    for key, label in FILTER_GROUPS:
        vals = sorted({tool_filters(t)[key] for t in tools_used})
        chips = '<button type="button" class="chip" data-f="" aria-pressed="true">any</button>' + "".join(
            f'<button type="button" class="chip" data-f="{esc(v)}" aria-pressed="false">{esc(v)}</button>' for v in vals)
        rows.append(f'<div class="tc-frow" role="group" aria-label="{esc(label)}" data-key="{key}">'
                    f'<span class="tc-flab">{esc(label)}</span>{chips}</div>')
    x_mark = '<svg width="12" height="12" aria-hidden="true"><path d="M2.5,2.5 l7,7 M9.5,2.5 l-7,7" class="st-nr"/></svg>'
    body = f"""<h1>Tool cards</h1>
<p class="lede">One card per tool, best Fmax first: what it compares, whether it will work for your
proteins, where it is strong and weak, what it costs, and the exact command that was run.</p>
<details class="choices tc-filters"><summary><span class="choices-label">Filter the cards</span>
<span class="choices-now" id="tc-count"></span></summary>
<div class="tc-fbox">{"".join(rows)}</div></details>
<figure class="tcfig">
<figcaption><b>Figure 1.</b> Each tool's best Fmax at each of the {len(rel_order)} relatives, one bar per
relative {esc(rel_order_words(rel_order))}: Pfam domains, raw calls, all proteins, best
over the tool's settings. Every card uses the same scale, 0 to {ymax:.2f}. Hover a bar for its relative and number.</figcaption>
{legend_kinds(tools_used)}
<div class="legend"><span class="item">{x_mark} not run at that relative</span></div>
<ul class="tc-grid" id="tc-grid">{"".join(cards)}</ul>
</figure>
<script>
(function () {{
  const cards = [...document.querySelectorAll(".tc")], count = document.getElementById("tc-count");
  const want = {{}};
  function apply() {{
    let n = 0;
    for (const c of cards) {{
      const ok = Object.entries(want).every(([k, v]) => !v || c.dataset[k] === v);
      c.hidden = !ok; n += ok;
    }}
    count.textContent = "showing " + n + " of " + cards.length + " tools";
  }}
  document.querySelectorAll(".tc-frow").forEach((row) => row.querySelectorAll(".chip").forEach((b) => b.addEventListener("click", () => {{
    want[row.dataset.key] = b.dataset.f;
    row.querySelectorAll(".chip").forEach((x) => x.setAttribute("aria-pressed", x === b));
    apply();
  }})));
  apply();
}})();
</script>"""
    return layout("Tool cards", "tools", body, src, depth=1,
                  description="One card per protein annotation tool, best Fmax first, with its cost and a bar per relative searched.")


def parse_tables(md: str) -> list[tuple[list[str], list[list[str]]]]:
    """Every markdown table in coverage.md as (header, rows)."""
    out, cur = [], None
    for ln in md.splitlines():
        if ln.startswith("|"):
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            if cur is None:
                cur = (cells, [])
            elif not set("".join(cells)) <= set("-: "):
                cur[1].append(cells)
        elif cur is not None:
            out.append(cur)
            cur = None
    if cur is not None:
        out.append(cur)
    return out


def parse_coverage(md: str):
    """First table: tool by stratum axis. Returns (axes, rows) where rows are [tool, *counts]."""
    head, rows = parse_tables(md)[0]
    return head[1:], rows


AXIS_TITLES = {"target species": "whole set, per relative", "identity": "% identity to nearest relative",
               "plddt": "AlphaFold pLDDT", "feature_type": "Swiss-Prot feature type",
               "feature_length_bin": "region length"}
# One-line titles for the experiment ids in coverage.md, from data/inventory.md (the id is the
# hover text and stays in the page's JSON).
EXPERIMENT_TITLES = {
    "alphabet_ranking_three_cases": "Alphabets ranked on three test proteins (notebook 241)",
    "dark_set_accuracy": "Invertebrate dark set: accuracy of the best kmerseek setting",
    "dark_set_gain": "Invertebrate dark set: what kmerseek reaches",
    "dark_set_pfam_fmax": "Invertebrate dark set: Fmax against Pfam (notebook 291)",
    "deep_twilight": "Deep-twilight controls, 16 proteins (notebook 252)",
    "elm_cover": "Short linear motifs (ELM) covered (notebook 253)",
    "qfo_region_midi": "Domain finding, 9 relatives",
    "qfo_region_midi_plus": "Domain finding, 15 relatives (this leaderboard)",
    "qfo_region_mini": "Domain finding, test set of 2 relatives",
    "scope40_v040": "SCOPe40 structural classes, kmerseek 0.4 (notebook 277)",
    "scope40_v040_below_kmin": "SCOPe40, kmerseek 0.4, k below the usual minimum",
    "scope40_v040_shuffled_k2_seed1": "SCOPe40 control: alphabet letters shuffled",
}
# Rows of coverage.md with no tool card: a plain name instead of the id.
NO_CARD_NAMES = {"mmseqs2": "MMseqs2, search mode not recorded", "folddisco": "Folddisco", "reseek": "Reseek"}
NULL_CONTROL = "scope40_v040_shuffled_k2_seed1"


def does_not_apply(tool: str, col: str) -> bool:
    """A square no run could fill: hmmscan searches no relative, so it has no identity to a
    relative; the shuffled-alphabet SCOPe run is a control for kmerseek only."""
    return (tool == "hmmscan" and col == "identity") or (col == NULL_CONTROL and tool != "kmerseek")


def page_coverage(cov_md: str, src: dict, species_mya: dict) -> str:
    tables = parse_tables(cov_md)
    kinds = list(KINDS)

    def name(t):
        return TOOLS[t]["name"] if t in TOOLS else NO_CARD_NAMES.get(t, t)

    def tool_cell(t):
        if t in TOOLS:
            return (f'<span class="kind k-{TOOLS[t]["kind"]}" aria-hidden="true"></span>'
                    f'<a href="tools/{t}.html">{esc(TOOLS[t]["name"])}</a>')
        return f'<span title="{esc(t)}">{esc(name(t))}</span> <span class="muted">(no card)</span>'

    def order(t):
        return (kinds.index(TOOLS[t]["kind"]), list(TOOLS).index(t)) if t in TOOLS else (len(kinds), name(t))

    figs, todo, n_on, n_off, n_na = [], {}, 0, 0, 0
    for k, (head, rows) in enumerate(tables[:2]):
        cols = head[1:]
        labels = [(AXIS_TITLES.get(c, c) if k == 0 else EXPERIMENT_TITLES.get(c, c), c) for c in cols]
        th = "".join(f'<th scope="col" title="{esc(c)}">{esc(x)}</th>' for x, c in labels)
        body = []
        for r in sorted(rows, key=lambda r: order(r[0])):
            tds = []
            for c, v in zip(cols, r[1:]):
                on = v != ""
                na = not on and does_not_apply(r[0], c)
                n_on += on
                n_na += na
                n_off += not on and not na
                if not on and not na and k == 1:
                    todo.setdefault(c, []).append(r[0])
                if on:
                    mark = f'<span class="sq on" role="img" aria-label="measured"></span><br><span class="mono">{esc(v)}</span>'
                elif na:
                    mark = '<span class="na" role="img" aria-label="does not apply">–</span>'
                else:
                    mark = '<span class="nr" role="img" aria-label="not run">×</span>'
                tds.append(f'<td data-measured="{"1" if on else "na" if na else "0"}">{mark}</td>')
            body.append(f'<tr><td class="txt">{tool_cell(r[0])}</td>{"".join(tds)}</tr>')
        title = ("Table 2a. Tool by subset of the proteins" if k == 0 else "Table 2b. Tool by experiment")
        summed = ("experiments, query sets, relatives and truth sets" if k == 0
                  else "the query sets, relatives, truth sets and subsets of that experiment")
        figs.append(f"""<figure>
<figcaption><b>{title}.</b> The number under each filled square is the count of true instances scored,
summed over {summed}. Hover a column name for its id in the tables.</figcaption>
<div class="tablebox"><table class="cov{" cov-exp" if k == 1 else ""}"><thead><tr><th class="txt" scope="col">Tool</th>{th}</tr></thead>
<tbody>{"".join(body)}</tbody></table></div>
</figure>""")
    todo_html = "".join(
        f'<li><b title="{esc(c)}">{esc(EXPERIMENT_TITLES.get(c, c))}</b>: {esc(", ".join(name(t) for t in sorted(ts, key=order)))}</li>'
        for c, ts in todo.items())
    body = f"""<h1>Coverage</h1>
<p class="lede">What has been measured, across every experiment in the tables, and what has not.
Every × is a run still to do.</p>
<div class="legend" aria-label="Key to Tables 2a and 2b">
<span class="item"><span class="sq on"></span> measured, with the number of true instances scored</span>
<span class="item"><span class="nr">×</span> not run</span>
<span class="item"><span class="na">–</span> does not apply: hmmscan searches no relative, and the shuffled-alphabet run is a control for kmerseek only</span>
</div>
{legend_kinds(set(TOOLS))}
{"".join(figs)}
<p class="prose muted">{n_on} squares are filled, {n_off} are not run and {n_na} do not apply. Where each
experiment's output came from, and which outputs were left out and why:
<a href="data/inventory.md">inventory.md</a> (copied with the tables).</p>
<h2>What to run next</h2>
<p class="prose">One line per experiment: the tools it has not run yet.</p>
<ul class="prose todo">{todo_html}</ul>"""
    return layout("Coverage", "coverage", body, src,
                  description="Which tool, subset and experiment combinations have been measured and which have not been run.")


def page_methods(res: pl.DataFrame, comp: pl.DataFrame, src: dict, prov: dict) -> str:
    ex = best_by(res, tool="hmmer3_phmmer", setting="default", truth_set="pfam", calls="raw",
                 species="zebrafish", stratum_axis="all", stratum="all").row(0, named=True)
    P, R = ex["fmax_precision"], ex["fmax_recall"]
    f_check = 2 * P * R / (P + R)
    assert abs(f_check - ex["fmax"]) < 1e-9, (f_check, ex["fmax"])
    cx = comp.filter((pl.col("tool") == "hmmer3_phmmer") & (pl.col("species") == "zebrafish")).row(0, named=True)
    cpu_h = cx["cores"] * cx["wall_s"] / 3600
    assert abs(cpu_h - cx["cpu_h_total"]) < 1e-9, (cpu_h, cx["cpu_h_total"])
    per1k = cpu_h / cx["n_relatives"] * 1000 / cx["n_queries"]
    cost = per1k * R_CPU
    ded = best_by(res, tool="hmmer3_phmmer", setting="default", truth_set="pfam", calls="dedup",
                  species="zebrafish", stratum_axis="all", stratum="all").row(0, named=True)
    n_rows = res.height
    pr5 = comp.filter(pl.col("tool") == "prostt5").row(0, named=True)
    long = pl.read_parquet(Path(DATA_DIR) / "results.parquet")
    n_long = long.height
    n_exp_long = long.filter(pl.col("experiment") == EXPERIMENT).height
    exps = sorted(long["experiment"].unique().to_list())
    pfn_all = best_by(res, tool="hmmer3_phmmer", truth_set="pfamn", calls="raw", species="zebrafish",
                      stratum_axis="all").row(0, named=True)["n_truth_instances"]
    pfn_nh = best_by(res, tool="hmmer3_phmmer", truth_set="pfamn", calls="raw", species="zebrafish",
                     stratum_axis="identity", stratum="no_homolog").row(0, named=True)["n_truth_instances"]
    pfn = {"all": pfn_all, "nohom": pfn_nh}
    n_rel = res.filter(pl.col("species") != "all")["species"].n_unique()

    body = rf"""<h1>Methods</h1>
<p class="lede">How each score and cost on this site is computed, which proteins were searched, and
how the truth sets were built.</p>
<nav class="toc" aria-label="Contents"><b>Contents</b>
<ol><li><a href="#data">The benchmark run</a></li><li><a href="#no-winner">No tool wins everywhere</a></li>
<li><a href="#calls">From a search hit to a call</a></li>
<li><a href="#scores">The scores</a>: <a href="#fmax">Fmax</a>, <a href="#auprc">AUPRC</a>, <a href="#iou">median IoU</a>,
<a href="#n">true instances</a>, <a href="#cpu">core-hours</a>, <a href="#cost">estimated cost</a></li>
<li><a href="#truth">Truth sets</a></li><li><a href="#dedup">Dedup against raw</a></li>
<li><a href="#subsets">Subsets of the proteins</a></li><li><a href="#checks">Checks on this site</a></li>
<li><a href="#experiments">Experiments in the tables</a></li></ol></nav>
<div class="prose">
<h2 id="data">The benchmark run</h2>
<p>Everything on this site except the coverage page comes from one experiment,
<span class="mono">{EXPERIMENT}</span>, a single benchmark run. It searched {prov['n_queries']:,} human
proteins against the proteomes of {n_rel} relatives from the Quest for Orthologs reference proteomes,
release 2020_04. Most of the human proteins are encoded on chromosome 6; the rest are immune
proteins (KIR and LILR receptors, CD1 and MR1, and β2-microglobulin). The copied tables do not
record how many are in each group. The run finished on 2026-09-13 on Stanford's Sherlock
cluster. It has {n_exp_long:,} of the {n_long:,} rows of <span class="mono">results.parquet</span>,
one per score. The leaderboard groups them into {n_rows:,} combinations of tool, setting, relative,
truth set, raw or dedup, and subset. Its cost is in <span class="mono">compute.parquet</span>. Both
tables were built by
<span class="mono">make leaderboard-tables</span> in {src['repo']} at commit
<a class="mono" href="https://github.com/{src['repo']}/commit/{src['commit']}">{src['commit'][:8]}</a>.
These {prov['n_queries']:,} proteins are not a random sample of the human proteome.</p>

<h2 id="no-winner">No tool wins everywhere</h2>
<p>No tool wins everywhere. So this board ranks tools within a subset of the proteins (by
relative, identity, region length, disorder or feature type), and gives no overall score. Two
earlier community benchmarks found the same. The Quest for Orthologs benchmark ran 15 orthology
methods on 20 benchmarks. Among the balanced methods, "no method consistently obtained a better
balance than the other methods across all data sets". Its authors chose not to compute a combined
score, because the right trade-off depends on the application (Altenhoff et al. 2016). The Massive Text Embedding Benchmark (MTEB), for text embedding models,
reported that "no particular text embedding method dominates across all tasks" (Muennighoff et al.
2023).</p>
<p class="refs">Altenhoff, Adrian M., Brigitte Boeckmann, Salvador Capella-Gutierrez, et al.
"Standardized Benchmarking in the Quest for Orthologs." <i>Nature Methods</i> 13, no. 5 (2016):
425–30. <a href="https://doi.org/10.1038/nmeth.3830">https://doi.org/10.1038/nmeth.3830</a>.<br>
Muennighoff, Niklas, Nouamane Tazi, Loic Magne, and Nils Reimers. "MTEB: Massive Text Embedding
Benchmark." In <i>Proceedings of the 17th Conference of the European Chapter of the Association for
Computational Linguistics</i>, 2014–37. Dubrovnik, Croatia, 2023.
<a href="https://doi.org/10.18653/v1/2023.eacl-main.148">https://doi.org/10.18653/v1/2023.eacl-main.148</a>.</p>

<h2 id="calls">From a search hit to a call</h2>
<p>Each tool reports aligned regions: a stretch of the human protein matched to a stretch of a
protein in the relative. A region becomes a <b>call</b> when it covers at least half of an annotated
domain or feature on the relative's protein. The call names that annotation and sits on the human
protein where the region is. A call is <b>correct</b> when the human protein carries that
annotation and the call overlaps it with an IoU of at least 0.5. <b>IoU</b> (intersection over
union) is the number of residues the call and the true region share, divided by the number of
residues covered by either:</p>
<p>\[\text{{IoU}} = \frac{{|\text{{call}} \cap \text{{true region}}|}}{{|\text{{call}} \cup \text{{true region}}|}}\]</p>
<p>For example, a call on residues 30–129 and a true domain on residues 10–109 share 80 residues
(30–109). Together they cover 120 (10–129). So IoU = 80 / 120 = 0.67, and the call counts as
correct. A call of the right family that misses the region counts as wrong. The IoU is drawn as
bars on a protein line, with ends you can move, on <a href="explainers.html#iou">How the scores work</a>.</p>
<p>A true annotation is <b>reachable</b> when its family occurs somewhere in the relative's
proteome. A family the relative lacks cannot be copied from it by any tool.</p>

<h2 id="scores">The scores</h2>
<h3 id="fmax">Fmax</h3>
<p><b>Fmax</b> is the best balance of precision and recall a tool reaches at any score cutoff,
counted protein by protein as in the CAFA function-prediction challenge. At a cutoff \(t\), the
precision of protein \(i\) is its correct calls divided by all its calls scoring at least \(t\). Its
recall is the true annotations found divided by its true annotations. Precision is averaged over
the \(m(t)\) proteins with at least one call at or above \(t\). Recall is averaged over all \(n\)
proteins with an annotation, and counts every true annotation, whether or not the relative could
supply it:</p>
<p>\[\text{{Pr}}(t) = \frac{{1}}{{m(t)}}\sum_{{i}} \frac{{\text{{correct}}_i(t)}}{{\text{{calls}}_i(t)}}\]</p>
<p>\[\text{{Rc}}(t) = \frac{{1}}{{n}}\sum_{{i}} \frac{{\text{{found}}_i(t)}}{{\text{{true}}_i}}\]</p>
<p>\[F_{{\max}} = \max_t \frac{{2\,\text{{Pr}}(t)\,\text{{Rc}}(t)}}{{\text{{Pr}}(t)+\text{{Rc}}(t)}}\]</p>
<p>The cutoffs tried are 101 evenly spaced quantiles of the tool's call scores. Worked example:
phmmer searching zebrafish, Pfam truth set, raw calls, all proteins. At the best cutoff (bit score
on the tool's bit score) precision is {P:.4f} and recall is {R:.4f}. So
F = 2 × {P:.4f} × {R:.4f} / ({P:.4f} + {R:.4f}) = {f_check:.4f}, the Fmax in Table 1.</p>

<h3 id="auprc">AUPRC</h3>
<p><b>AUPRC</b>, the area under the precision-recall curve, sums precision over every cutoff,
weighted by the recall gained at that cutoff:</p>
<p>\[\text{{AUPRC}} = \sum_k (R_k - R_{{k-1}})\,P_k\]</p>
<p>Unlike Fmax, AUPRC pools all proteins. \(P_k\) is correct calls ÷ all judged calls at cutoff
\(k\). \(R_k\) is reachable annotations found ÷ all reachable annotations. A tool that makes many calls per protein pays for each wrong one. A call in a
part of a protein that no annotation covers is left out of \(P_k\), since nobody has looked there;
Fmax counts it. Worked
example, same row: AUPRC = {ex['auprc']:.4f}. The curve itself is in the source run
(<span class="mono">all_domain_curves.parquet</span>), not in these tables, so the sum cannot be
recomputed here.</p>

<h3 id="iou">Median IoU</h3>
<p><b>Median IoU</b> is the middle IoU of all correct calls, with single-residue annotations left
out (an IoU against one residue says nothing about where the ends are). Worked example, same row:
{ex['median_iou_tp']:.3f}, so half of phmmer's correct calls on zebrafish overlap their true domain
by at least that much. A correct call on a region longer than one residue has an IoU of at least
0.5. So the run stores 0 only when there is no such call. Either there are no correct calls at
all, or the correct calls are on single residues only (the 1 aa length bin and the ACT_SITE,
BINDING and SITE features). The pages show the first case as <i>no correct calls</i> and the
second as <i>correct calls only on single residues</i>, never as 0. They tell the two apart by
Fmax: an Fmax above 0 means at least one call was correct.</p>

<h3 id="n">True instances</h3>
<p><b>True instances</b> (n) is the number of annotations on the human proteins in the subset,
reachable or not. Same row: {ex['n_truth_instances']:,}.</p>

<h3 id="cpu">Core-hours per 1,000 queries</h3>
<p><b>Core-hours</b> are taken from the run's Nextflow task records: the cores a task reserved
multiplied by its wall time, so a task that sat partly idle still counts all its cores. The tables
record them per tool and setting, summed over every relative searched, with building the databases
as separate rows. The leaderboard divides the search's core-hours by the number of relatives
searched and by the {cx['n_queries']} queries, times 1,000: the mean cost of searching one relative.
Worked example: phmmer reserved {cx['cores']} cores for {cx['wall_s']:,.0f} s over its
{cx['n_relatives']} searches. That is {cx['cores']} × {cx['wall_s']:,.0f} / 3,600 = {cpu_h:.4f}
core-hours, and {cpu_h:.4f} / {cx['n_relatives']} × 1,000 / {cx['n_queries']} = {per1k:.4f}
core-hours per 1,000 queries per relative. Peak RAM is the largest resident memory of any of
those tasks. The tool cards give the database builds separately.</p>
<p>Numbers marked <span class="prov measured">measured</span> come from these records, on Sherlock
(partition hns). No search used a GPU. ProstT5 has to turn sequences into 3Di shape letters
before it can search; that step ran on one GPU as part of building its databases, queries
included, and is on the ProstT5 card ({pr5['db_gpu_h']:.3g} GPU-hours over {pr5['db_n']}
databases), not in the per-query cost. Numbers marked
<span class="prov reported">reported</span> would come from a tool's submitter, on hardware they
name; none have been submitted yet.</p>

<h3 id="cost">Estimated cost</h3>
<p>\[\text{{cost\_usd\_est}} = \text{{CPU-hours}} \times r_\text{{cpu}} + \text{{GPU-hours}} \times r_\text{{gpu}}\]</p>
<p>per 1,000 queries. The default \(r_\text{{cpu}}\) = ${R_CPU:.6g} per CPU-hour is half of AWS's
on-demand price for a c7i.large (2 vCPU, $0.08925 per hour). The default \(r_\text{{gpu}}\) =
${R_GPU} per GPU-hour is a g5.xlarge (one NVIDIA A10G). Both are for US East (N. Virginia), Linux,
from the <a href="{PRICE_PAGE}">EC2 on-demand price page</a> (<a href="{PRICE_FILE}">price
file</a>, published 2026-10-12, read 2026-10-15). A core-hour is priced as one CPU-hour. Memory is not
priced. Worked example: {per1k:.4f} × {R_CPU:.6g} + 0 ×
{R_GPU} = ${cost:.5f} for phmmer against zebrafish. Change either price on the leaderboard page.</p>

<h2 id="truth">Truth sets</h2>
<h3>Pfam</h3>
<p>Every Pfam-A domain with coordinates on the human proteins, from the Pfam regions file for
UniProt reference proteomes. The same file supplies the domains on each relative's proteins, which
are what a call copies. Pfam families are defined by profile models, so tools that build profiles
(jackhmmer, HH-suite) are graded against the kind of model they build, and hmmscan against its own
models. Release of the regions file: not recorded in the run; the hmmscan models are Pfam-A 38.2.
hmmscan was scored on raw calls only, so it shows "not run" under dedup.</p>
<h3>Pfam-N</h3>
<p>Pfam-N is an extension of Pfam made by a deep-learning model: it labels regions the Pfam
profile models missed. It is only published for Pfam 35.0. This is the one truth set the profile
models did not produce. Its regions were filtered to the human proteins queried here. Few of them
can be copied from any relative: against zebrafish, {pfn['nohom']:,} of the {pfn['all']:,} Pfam-N
regions have no domain of the same family in the zebrafish proteome. So every tool scores low on
Pfam-N, and the set says little about which tool is better.</p>
<h3>Swiss-Prot features</h3>
<p>Curated features from Swiss-Prot entries, of 12 types (the feature keys on the leaderboard,
such as TRANSMEM for a membrane-spanning segment and ACT_SITE for an active site residue). They
come from the literature, not from any sequence or structure model. Single-residue features
(ACT_SITE, BINDING, SITE) are scored by whether the call covers the residue, since an IoU against
one residue cannot reach 0.5. A Swiss-Prot feature counts as reachable when the relative has a
protein that carries the same feature type and shares a Pfam family with the human protein.</p>

<h2 id="dedup">Dedup against raw</h2>
<p><b>Raw</b> keeps every call a tool makes. One human region that matches 300 mouse kinases
becomes 300 kinase calls on nearly the same residues, and all but one count as wrong.
<b>Dedup</b> merges calls of the same family on the same human protein that overlap each other with
an IoU of 0.5 or more, keeping the highest-scoring one. Two copies of a domain side by side
barely overlap, so both are kept. Comparing the two separates finding a domain from reporting it many times.
For phmmer on zebrafish (Pfam, all proteins): Fmax {ex['fmax']:.3f} raw and {ded['fmax']:.3f} dedup.
The same calls drawn raw and after dedup, on a made-up protein: <a href="explainers.html#dedup">How the scores work</a>.</p>

<h2 id="subsets">Subsets of the proteins</h2>
<p>The tables call each subset a stratum. There are four ways to cut the proteins:</p>
<ul>
<li><b>Nearest reviewed relative, % identity</b>: for each human domain, the percent identity to the
closest domain of the same family in the relative. Bins: 0–20, 20–30, 30–40, 40–60 and 60–100%. A
domain whose family is absent from the relative is "no same-family domain".</li>
<li><b>Nearest reviewed relative, Mya</b>: the relative searched, labelled with how many million
years ago its line split from ours: the pairwise divergence time from human in
<a href="https://timetree.org/">TimeTree 5</a>, rounded to the nearest million years. TimeTree
gives no time for human and a bacterium, so E. coli has none.</li>
<li><b>Region type and length</b>: the length of the annotated region, in residues. Bins up to 60
aa are short features; 61 aa and up is the size of a folded domain. Whole-protein annotations
were not scored in this run.</li>
<li><b>Disorder</b>: the mean AlphaFold confidence (pLDDT, 0–100) of the human protein. Low pLDDT
marks regions without a fixed shape.</li>
<li><b>Feature type</b>: the Swiss-Prot feature key, for the Swiss-Prot truth set only.</li>
</ul>
<p>The run computed each subset on its own, so two subsets cannot be combined; a combination with
no row shows × and "not run".</p>

<h2 id="checks">Checks on this site</h2>
<p>A headless browser renders every filter combination and confirms that each × on the
leaderboard is a missing row in results.parquet, and that every number matches the tables. The
source table checksums are in <span class="mono">data/SOURCE.json</span> of the site repository.</p>
<h2 id="experiments">Experiments in the tables</h2>
<p>The tables hold {n_long:,} rows from {len(exps)} experiments; the coverage page counts all of
them. Their ids: {", ".join(f'<span class="mono">{esc(e)}</span>' for e in exps)}.</p>
</div>
<link rel="stylesheet" href="assets/katex/katex.min.css">
<script defer src="assets/katex/katex.min.js"></script>
<script defer src="assets/katex/contrib/auto-render.min.js"
 onload="renderMathInElement(document.body,{{delimiters:[{{left:'\\[',right:'\\]',display:true}},{{left:'\\(',right:'\\)',display:false}}]}})"></script>
"""
    return layout("Methods", "methods", body, src,
                  description="How each score is defined and computed, how the truth sets were built, and what dedup means.")


def page_pileup(src: dict, res: pl.DataFrame, prov: dict, data: dict, check: dict, ro: dict) -> str:
    regions = data["regions"]
    import statistics
    med_len = statistics.median(r["len"] for r in regions)
    n_all = len(regions)
    n_past = sum(r["past_end"] for r in regions)
    n_in = n_all - sum(1 for r in regions if r["s"] > r["plen"])
    n_fam = len({r["fam"] for r in regions})
    n_prot = len({r["acc"] for r in regions})
    past = "".join(f'<li>{r["acc"]} {r["fam"]}: {r["s"]}–{r["e"]} on a {r["plen"]} aa sequence</li>'
                   for r in regions if r["past_end"])
    km = best_by(res, tool="kmerseek", truth_set="pfam", calls="raw", species="zebrafish",
                 stratum_axis="all").sort(["fmax", "setting"], descending=[True, False]).row(0, named=True)
    # One panel per tool by default: its only setting, or kmerseek's setting with the highest
    # Fmax against zebrafish. Figure 10 draws the same panels.
    default = [(q["tool"], q["setting"]) for q in data["panels"]
               if q["tool"] != "kmerseek" or q["setting"] == km["setting"]]
    hist = [{"tool": t, "setting": st} for t, st in default]
    early = {}
    for h in hist:
        pnl = next(q for q in data["panels"] if q["tool"] == h["tool"] and q["setting"] == h["setting"])
        ds = [row["d"] for row in pnl["rows"] if row["d"] is not None]
        h["n"] = len(ds)
        early[h["tool"]] = sum(d < -100 for d in ds) / len(ds)
    aligners = [t for t, _ in default if TOOLS[t]["kind"] in ("letters", "profiles", "shapes")]
    lo_t = min(aligners, key=lambda t: early[t]); hi_t = max(aligners, key=lambda t: early[t])
    # The page's lead sentence, from the same readouts the page computes (scripts/pileup.py).
    R = {k: v["all_domains"] for k, v in ro["panels"].items()}
    al = [R[f"{t}|{st}"] for t, st in default if t in aligners]
    kr = R[f"kmerseek|{km['setting']}"]
    hs = R.get("hmmscan|direct")
    lead = (f"The {len(al)} aligners cover {100 * min(a['mean_coverage_inside'] for a in al):.0f}% to "
            f"{100 * max(a['mean_coverage_inside'] for a in al):.0f}% of each domain on average, and also "
            f"{100 * min(a['mean_coverage_flanks'] for a in al):.0f}% to {100 * max(a['mean_coverage_flanks'] for a in al):.0f}% "
            "of the 100 residues either side: their calls run well past the domain's ends. "
            f"kmerseek ({esc(setting_label('kmerseek', km['setting']))}, its best setting by Fmax) misses "
            f"{100 * kr['missed_fraction']:.0f}% of the domains.")
    if hs:
        lead += (f" hmmscan, graded on its own models, covers {100 * hs['mean_coverage_inside']:.0f}% inside "
                 f"and {100 * hs['mean_coverage_flanks']:.0f}% of the flanks.")
    scores = {t: ("region_enrichment (shared k-mers over the number expected by chance)" if t == "kmerseek" else "bit score")
              for t in TOOLS}
    words = {f'{q["tool"]}|{q["setting"]}': setting_words(q["tool"], q["setting"]) for q in data["panels"]}
    labels = {f'{q["tool"]}|{q["setting"]}': setting_label(q["tool"], q["setting"]) for q in data["panels"]}
    meta = {
        "kinds": {k: v["color"] for k, v in KINDS.items()},
        "names": {t: TOOLS[t]["name"] for t in TOOLS}, "words": words, "labels": labels, "scores": scores,
        "median_len": med_len, "focus_default": "hmmer3_phmmer", "hist": hist,
        "default_panels": [f"{t}|{st}" for t, st in default],
    }
    curve_key = ('<svg width="28" height="10" aria-hidden="true"><line x1="1" y1="5" x2="27" y2="5" stroke="var(--k-letters)" stroke-width="2"/></svg>'
                 ' share of all domains covered, in the tool\'s colour')
    cap_key = ('<svg width="28" height="10" aria-hidden="true"><line x1="1" y1="5" x2="27" y2="5" stroke="var(--ink)" stroke-width="1.2"/></svg>'
               ' thin line: the same, one domain per Pfam family')
    body = f"""
<h1>Region pile-up</h1>
<p class="lede">Where each tool's calls fall around {n_all:,} known Pfam domains on {n_prot:,} human
proteins, searching the zebrafish proteome ({mya_label(SPECIES_MYA["zebrafish"]):,} million years since the split).</p>
<p class="verdict" id="pu-lead">{lead}</p>
<p id="pu-loading" class="muted">Loading {len(data["panels"])} panels…</p>
<figure>
<figcaption><b>Figure 9.</b> Calls piled up around each domain, one panel per tool. Top of each panel:
the share of domains covered at each position, from 100 residues before the domain's start to 100
after its end. Below: every domain as a row; stronger colour means more of the rows there are covered.
Click a panel to draw it full size below.
<button type="button" class="chip" id="pu-all" aria-pressed="false">Show all {len(data["panels"])} settings</button></figcaption>
<div class="legend" aria-label="Key to Figure 9">
<span class="item">{curve_key}</span><span class="item">{cap_key}</span>
<span class="item"><svg width="14" height="12" aria-hidden="true"><line x1="7" y1="0" x2="7" y2="12" stroke="var(--muted)" stroke-dasharray="3 2"/></svg> the domain's start and end</span>
<span class="item"><span class="sw" style="background:var(--surface);border:1px solid var(--hair)"></span> not covered</span>
<span class="item"><span class="sw" style="background:var(--null)"></span> past the protein's end</span>
</div>
<div id="pu-grid" class="pu-grid"></div>
</figure>
<details class="prose pu-howto"><summary>How to read a panel: covered, missed, inside, flanks, edge sharpness</summary>
<p>A <b>call</b> is a stretch of a human protein that a tool matched to a zebrafish protein carrying a
Pfam domain. It names that domain and spans the whole matched stretch, which can be much wider than
the domain. A residue is <b>covered</b> when a call naming that row's own Pfam family spans it, at any
score the tool reported (E-value cutoff 10). Covered does not mean correctly placed: a call that runs
far past the domain still covers it.</p>
<p>Each row shows 100 residues before the domain, then the domain stretched or squeezed into 100
equal steps, then 100 residues after it, as deepTools does for genes. The median domain is
{med_len:.0f} aa, so one step is about {med_len / 100:.1f} aa.</p>
<p>The numbers under each panel: <b>missed</b> is the share of domains with no covered residue inside
(lower is better). <b>Inside</b> is the mean coverage inside the domain (higher is better).
<b>Flanks</b> is the mean coverage in the 100 residues either side; for placing domain ends lower is
better, though a second copy of the same domain nearby counts as covered too. <b>Edge sharpness</b>
is how many residues the curve takes to fall from 0.8 to 0.2 at each end of the domain: fewer
residues means a sharper edge. "None" means the curve never reaches 0.8 inside the domain, or never
falls to 0.2 within 100 residues outside it; the chosen panel below says which.</p>
<p>Rows from one protein or one family are not independent, so the thin line repeats the curve
with one domain per Pfam family ({n_fam:,} families).</p>
</details>
<figure>
<figcaption><b>Figure 9, chosen panel:</b> <span id="pu-ftitle"></span>. <span id="pu-caption-numbers"></span>
Hover a row to see its protein, domain, coordinates, calls and residues.</figcaption>
<div class="metricpick" role="group" aria-label="Row order">
<span>Order rows (in every panel) by:</span>
<button type="button" class="chip" data-sort="coverage">coverage in this panel</button>
<button type="button" class="chip" data-sort="identity">% identity to zebrafish</button>
<button type="button" class="chip" data-sort="length">domain length</button>
<button type="button" class="chip" data-sort="plddt">mean pLDDT</button>
</div>
<div class="legend" aria-label="Key to the chosen panel">
<span class="item"><span class="sw" style="background:var(--strip-a)"></span><span class="sw" style="background:var(--strip-b)"></span>
left strip: rows from one protein share a shade</span>
<span class="item">Colour by: <button type="button" class="chip" data-color="coverage">covered or not</button>
<button type="button" class="chip" data-color="score">call score</button></span>
<span class="item" id="pu-scorebar" hidden><span class="name"></span><span class="lo mono"></span><span class="bar"></span><span class="hi mono"></span>
highest score of the calls over the residue, 5th to 95th percentile of this tool's scores</span>
</div>
<div class="plotbox pu-focus">
<svg id="pu-fcurve" class="pu-fcurve" role="img" aria-label="share of domains covered at each position"></svg>
<div class="pu-scroll"><canvas id="pu-fheat" class="pu-fheat" aria-label="one row per domain"></canvas>
<div id="pu-tip" class="tip pu-tip" hidden></div></div>
</div>
<p id="pu-fro" class="prose muted"></p>
<p class="axisnote">Identity is the best pairwise identity between the human domain and any zebrafish
domain of its family. Domains with no such domain, or no matching AlphaFold model, sort last.</p>
</figure>
<h2>Where calls start</h2>
<figure>
<figcaption><b>Figure 10.</b> Among the aligners, {100 * early[lo_t]:.0f}% ({esc(TOOLS[lo_t]["name"])}) to
{100 * early[hi_t]:.0f}% ({esc(TOOLS[hi_t]["name"])}) of calls start more than 100 residues before the
domain, because a call spans the whole matched stretch. hmmscan's calls start on the domain, whose
ends its own models define. Each panel is one tool. For each domain with a call, it takes the call
that overlaps the domain most (ties to the higher score, then the earlier start). It plots that
call's first residue minus the domain's first residue, as a share of the tool's domains with a call
(n on the panel), in 5-residue steps. The end bins hold everything beyond 100 residues. kmerseek is drawn at
{esc(setting_label("kmerseek", km["setting"]))} (<span class="mono">{esc(km["setting"])}</span>). All panels share one
vertical scale except hmmscan's, whose peak near 1 would flatten the rest. Where calls end is not drawn.</figcaption>
<div class="legend" aria-label="Key to Figure 10"><span class="item"><svg width="28" height="10" aria-hidden="true"><line x1="1" y1="5" x2="27" y2="5" stroke="var(--k-letters)" stroke-width="2"/></svg> this panel's tool, in its colour</span>
<span class="item"><svg width="28" height="10" aria-hidden="true"><line x1="1" y1="5" x2="27" y2="5" stroke="var(--hair)" stroke-width="1.5"/></svg> the other tools, for comparison</span>
<span class="item"><svg width="14" height="12" aria-hidden="true"><line x1="7" y1="0" x2="7" y2="12" stroke="var(--muted)" stroke-dasharray="3 2"/></svg> 0: the call starts on the domain's first residue</span></div>
<div id="pu-hist" class="pu-hist" role="group" aria-label="call start minus domain start, one panel per tool"></div>
</figure>
<h2>Counts and sources</h2>
<div class="prose muted">
<p>{n_in:,} of the {n_all:,} domains have at least one residue inside their protein's sequence; the
readouts use those. With one domain per Pfam family ({n_fam:,} families), the domain kept is the
first by protein accession, then start. {n_past} Pfam domains run partly or wholly past the end of
the protein sequence this run searched, so their Pfam coordinates belong to a different version of
the sequence.</p>
<details><summary>The {n_past} domains that run past their protein's end</summary><ul class="mono">{past}</ul></details>
<p>The pile-up's own tables (domains and calls) come from {src['repo']} at commit
<a class="mono" href="https://github.com/{src['repo']}/commit/{src['pileup_commit']}">{src['pileup_commit'][:8]}</a>
(<span class="mono">results/leaderboard/pileup</span>).</p>
</div>
<script type="application/json" id="pu-meta">{json.dumps(meta)}</script>
<script src="assets/pileup.js"></script>
"""
    return layout("Region pile-up", "pileup", body, src,
                  description="Where each tool's calls fall around known Pfam domains, residue by residue.")
from html.parser import HTMLParser


class _Text(HTMLParser):
    """Headings, paragraphs, captions, list items and table cells of <main>, in order."""
    KEEP = {"h1", "h2", "h3", "p", "figcaption", "li", "dt", "dd", "summary", "pre", "td", "th"}

    def __init__(self):
        super().__init__()
        self.out, self.stack, self.in_main, self.skip = [], [], False, 0

    def handle_starttag(self, tag, attrs):
        if tag == "main":
            self.in_main = True
        if tag in ("script", "style", "svg"):
            self.skip += 1
        if self.in_main and tag in self.KEEP:
            self.stack.append([tag, ""])

    def handle_endtag(self, tag):
        if tag in ("script", "style", "svg"):
            self.skip -= 1
        if tag == "main":
            self.in_main = False
        if self.stack and self.stack[-1][0] == tag:
            t, text = self.stack.pop()
            text = " ".join(text.split())
            if text:
                self.out.append({"tag": t, "text": text})

    def handle_data(self, data):
        if self.stack and not self.skip:
            self.stack[-1][1] += data


def page_json(out: Path, rel: str, src: dict, data) -> None:
    """Write <page>.json beside <page>: the page's text in reading order plus its data, so an
    agent gets what a person gets."""
    html_text = (out / rel).read_text()
    parser = _Text()
    parser.feed(html_text)
    title = re.search(r"<title>(.*?)</title>", html_text).group(1)
    doc = {"page": rel, "title": html.unescape(title), "source": src, "text": parser.out, "data": data}
    (out / (rel + ".json")).write_text(json.dumps(doc, indent=1, ensure_ascii=False, allow_nan=False) + "\n")
    name = rel.split("/")[-1]
    link = f'<p>This page as JSON, for scripts and agents: <a class="mono" href="{name}.json">{name}.json</a></p>\n'
    (out / rel).write_text(html_text.replace("</div></footer>", link + "</div></footer>", 1))



def fmt_date(d):
    if not d:
        return ""
    y, m, dd = d[:10].split("-")
    return f"{int(dd)} {['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'][int(m) - 1]} {y}"


# Claim status marks. The grey x is kept for "not run" across the site, so a disproved
# claim gets its own mark.
MARK = {"alive": "●", "untested": "○", "retired": "⊘"}
STATUS_WORDS = {"alive": "Holds", "untested": "Not yet tested", "retired": "Disproved"}


def plain_refs(text: str) -> str:
    """"nb 239" -> "notebook 239" in text shown to readers."""
    return re.sub(r"\bnb (\d{3})", r"notebook \1", text)


def claims_card(cl: dict, up: str = "") -> str:
    """One tool's claims, as resolved by scripts/claims.py, grouped by where they stand: the
    claim, the result that would disprove it, the verdict with its number, and the evidence."""
    c = cl["counts"]
    pills = (f'<span class="cl-pill alive">{MARK["alive"]} {c["alive"]} hold{"s" if c["alive"] == 1 else ""}</span>'
             f'<span class="cl-pill untested">{MARK["untested"]} {c["untested"]} not yet tested</span>'
             f'<span class="cl-pill retired">{MARK["retired"]} {c["retired"]} disproved</span>'
             + (f'<span class="cl-pill plain">last disproved {fmt_date(cl["last_kill"])}</span>' if cl["last_kill"] else ""))
    groups = []
    for status in ("alive", "untested", "retired"):
        rows = []
        for x in [x for x in cl["claims"] if x["status"] == status]:
            ev = x.get("evidence") or {}
            href = ev.get("href", "")
            if href and not href.startswith("http"):
                href = up + href
            when = f' {fmt_date(x["date"])}' if status == "retired" and x.get("date") else ""
            verdict = f'<span class="cl-st {status}">{STATUS_WORDS[status].lower()}{when}</span> · {esc(plain_refs(x["verdict"]))}'
            if href:
                verdict += f' · <a href="{esc(href)}">{esc(plain_refs(ev["label"]))}</a>'
            ctx = "".join(f'<li>{esc(plain_refs(k["label"]))}: <span class="mono">{esc(k["text"])}</span></li>' for k in x["context"])
            more = (f'<details class="cl-more"><summary>Scope and related numbers</summary>'
                    f'<p>{esc(plain_refs(x.get("scope") or ""))}</p><ul>{ctx}</ul>'
                    + (f'<p class="muted">Source: <span class="mono">{esc(x["source"])}</span></p>' if x.get("source") else "")
                    + '</details>')
            rows.append(f'<li class="cl-row" data-status="{status}"><span class="cl-mark" aria-hidden="true">{MARK[status]}</span>'
                        f'<div><div class="cl-claim">{esc(plain_refs(x["claim"]))}</div>'
                        f'<div class="cl-kill">Disproved if: {esc(plain_refs(x["kill"]))}</div>'
                        f'<div class="cl-verdict">{verdict}</div>{more}</div></li>')
        if rows:
            groups.append(f'<h3 class="cl-group">{MARK[status]} {STATUS_WORDS[status]} ({len(rows)})</h3>'
                          f'<ul class="cl-rows">{"".join(rows)}</ul>')
    return f"""<section class="cl-card" aria-label="Claims about {esc(cl['title'])}">
<div class="cl-head"><h2 class="cl-title">{esc(cl['title'])}</h2>
<span class="cl-file"><a class="mono" href="{up}claims/{cl['tool']}.yaml">claims/{cl['tool']}.yaml</a>, updated {fmt_date(cl['updated'])}</span></div>
<div class="cl-pills">{pills}</div>
<p class="cl-sub">each row: the claim, the result that would disprove it, the verdict with its number, the evidence</p>
{"".join(groups)}
<div class="cl-buttons"><a class="cl-btn" href="{up}claims-schema.html">claims.yaml schema ↗</a>
<a class="cl-btn" href="{up}claims-schema.html#badge">badge and CI check ↗</a></div>
</section>"""


def page_claims(src: dict, cards: list[dict]) -> str:
    who = ("Only kmerseek has filed claims so far." if [c["tool"] for c in cards] == ["kmerseek"]
           else f"{len(cards)} tools have filed claims so far.")
    body = f"""<h1>Claims</h1>
<p class="lede">What each tool's authors claim, the result that would disprove each claim, and where
it stands. {who}</p>
<div class="legend" aria-label="Key to the claims"><span class="item">{MARK["alive"]} holds: its test ran and did not disprove it</span>
<span class="item">{MARK["untested"]} not yet tested: the test has not run, or the data cannot decide it</span>
<span class="item">{MARK["retired"]} disproved: the test disproved it</span></div>
<p class="prose muted">The build recomputes every verdict from the tables or checks it against the
recorded notebook output, so nothing on a card is typed in by hand.</p>
{"".join(claims_card(c) for c in cards)}"""
    return layout("Claims", "claims", body, src, description="Each tool's claims, the result that would disprove each one, and its verdict.")


# The made-up example on explainers.html. The page's widgets compute every number from it, and
# tests/check_explainers.py recomputes them in Python. Coordinates are residues, inclusive;
# scores are the tool's score, higher meaning more confident. Family names are only labels.
TOY = {
    "note": "A made-up example for explaining the scores. Not benchmark data.",
    "proteins": {
        "A": {"name": "protein A", "short": "A", "length": 400,
              "true": [{"family": "PH", "start": 20, "end": 120},
                       {"family": "SH2", "start": 150, "end": 230},
                       {"family": "SH3", "start": 280, "end": 340}],
              "calls": [{"id": "a1", "family": "PH", "start": 25, "end": 118, "score": 95},
                        {"id": "a2", "family": "SH2", "start": 155, "end": 228, "score": 80},
                        {"id": "a3", "family": "PH", "start": 30, "end": 115, "score": 70},
                        {"id": "a4", "family": "PDZ", "start": 160, "end": 225, "score": 55},
                        {"id": "a5", "family": "SH3", "start": 285, "end": 338, "score": 40},
                        {"id": "a6", "family": "SH3", "start": 200, "end": 260, "score": 25},
                        {"id": "a7", "family": "PH", "start": 22, "end": 120, "score": 15}]},
        "B": {"name": "protein B", "short": "B", "length": 600,
              "true": [{"family": "Ig", "start": s, "end": s + 90} for s in (20, 130, 240, 350, 460)],
              "calls": [{"id": f"b{i + 1}", "family": "Ig", "start": s + 2, "end": s + 88, "score": sc}
                        for i, (s, sc) in enumerate(zip((20, 130, 240, 350, 460), (60, 58, 52, 50, 45)))]},
    },
    "iou": {"length": 200, "true": {"start": 10, "end": 109}, "call": {"start": 30, "end": 129}},
}


def legend_svg(kind: str) -> str:
    """Small key marks for the explainer figures, drawn with the same classes as the figures."""
    inner = {
        "true": '<line x1="0" y1="8" x2="30" y2="8" class="ex-backbone"/><rect x="7" y="2" width="16" height="12" rx="2" class="ex-true"/>',
        "kept": '<rect x="2" y="3" width="26" height="10" rx="2" class="ex-kept"/>',
        "drop": '<rect x="2" y="3" width="26" height="10" rx="2" class="ex-drop"/>',
        "cut": '<line x1="15" y1="0" x2="15" y2="16" class="ex-cut"/>',
        "ring": '<circle cx="15" cy="8" r="6" class="ex-ring"/>',
        "dot": '<circle cx="15" cy="8" r="4" class="ex-dot"/>',
        "dot-off": '<circle cx="15" cy="8" r="4" class="ex-dot-off"/>',
        "area": '<rect x="2" y="2" width="26" height="12" class="ex-area"/>',
        "area-off": '<rect x="2" y="2" width="26" height="12" class="ex-area-off"/>',
        "line": '<path d="M1,12 H10 V6 H20 V3 H29" class="ex-line"/>',
        "bracket": '<path d="M2,4 V11 H28 V4" class="ex-bracket"/>',
        "your": '<circle cx="15" cy="8" r="4.5" class="ex-ring-small"/>',
        "median": '<line x1="15" y1="0" x2="15" y2="16" class="ex-median"/>',
    }[kind]
    return f'<svg width="30" height="16" viewBox="0 0 30 16" aria-hidden="true">{inner}</svg>'


def page_explainers(src: dict, comp: pl.DataFrame) -> str:
    cx = comp.filter((pl.col("tool") == "hmmer3_phmmer") & (pl.col("species") == "zebrafish")).row(0, named=True)
    per1k = cx["cores"] * cx["wall_s"] / 3600 / cx["n_relatives"] * 1000 / cx["n_queries"]
    li = lambda kind, words: f'<span class="item">{legend_svg(kind)} {words}</span>'  # noqa: E731
    glyphs = '<span class="item"><b>✓</b> correct</span><span class="item"><b>✗</b> wrong</span>'
    toy_json = json.dumps(TOY, separators=(",", ":")).replace("</", "<\\/")
    body = r"""<h1>How the scores work</h1>
<p class="lede">The leaderboard ranks tools by three scores and one cost. Each section below
works one of them through on a made-up protein that you can change. Every number in Figures 4 to
7 is computed in your browser from that made-up example, not from the benchmark.</p>
<nav class="toc" aria-label="Contents"><b>Contents</b>
<ol><li><a href="#example">The made-up example</a></li>
<li><a href="#fmax">Fmax: the best single score cutoff</a> (Figure 4)</li>
<li><a href="#auprc">AUPRC: every cutoff at once</a> (Figure 5)</li>
<li><a href="#iou">Median IoU: are the ends in the right place?</a> (Figure 6)</li>
<li><a href="#dedup">Raw or dedup: one domain, many calls</a> (Figure 7)</li>
<li><a href="#cost">Cost: core-hours, GPU-hours and the price of an hour</a> (Figure 8)</li></ol></nav>

<h2 id="example">The made-up example</h2>
<div class="prose">
<p>Protein A is 400 residues (aa) long and has three true domains: the known annotations a tool
should find. A tool searched it and made seven <b>calls</b>. A call is a region the tool says
carries a domain of a given family, copied from a protein of the relative it searched. Each call
has the tool's <b>score</b>: higher means the tool is more sure. The family names are only labels.</p>
<p>A call is <b>correct</b> when it names the family of a true domain and overlaps that domain
with an IoU of at least 0.5 (Figure 6 shows IoU). Each true domain can be found once. If two calls
find the same domain, the higher-scoring one is correct and the other counts as wrong, as on the
leaderboard.</p>
<p>In this example every domain can be found in the relative, and every call lies on an annotated
part of the protein. So the small differences in what Fmax and AUPRC count (see
<a href="methods.html">Methods</a>) do not arise here.</p>
</div>

<h2 id="fmax">Fmax: the best single score cutoff</h2>
<div class="prose">
<p>To use a tool's calls you pick a <b>score cutoff</b>: you keep every call scoring at least that
much and drop the rest. At a cutoff, <b>precision</b> is the share of kept calls that are correct.
<b>Recall</b> is the share of true domains that a kept, correct call found. A high cutoff keeps
few calls, mostly correct; a low cutoff finds more domains but keeps more wrong calls. <b>F</b>
combines the two into one number, high only when both are high:</p>
<p>\[F = \frac{2 \times \text{precision} \times \text{recall}}{\text{precision} + \text{recall}}\]</p>
<p><b>Fmax</b> is the highest F over all cutoffs: the score of the tool at its best single cutoff.
Move the slider in Figure 4 to try each cutoff.</p>
</div>
<figure class="ex-fig" id="fig1">
<figcaption><span class="ex-made-up">made-up example</span><b>Figure 4.</b> Fmax on protein A. <details class="capmore"><summary>How to read Figure 4</summary>Top: the
protein as a line with its three true domains as boxes. Each call is a bar below it, at the
residues it covers, labelled with its family and score. Bottom: F at each cutoff. Move the cutoff
and read off which calls are kept and what F they give; the ring marks the highest F, which is the
Fmax.</details></figcaption>
<div class="ex-legend" aria-label="Key to Figure 4">
@@LEG1@@
</div>
<div class="ex-controls">
<label for="f1-cut"><span>Score cutoff: keep calls scoring at least <span class="v mono" id="f1-t"></span></span>
<input type="range" id="f1-cut" min="0" max="100" step="1" value="50"></label>
<button type="button" class="chip" id="f1-go">Go to the best cutoff</button>
</div>
<div class="ex-tiles" id="f1-read" aria-live="polite">
<div class="ex-tile"><span class="ex-tl">Calls kept</span><span class="ex-tv" id="f1-kept"></span>
<span class="ex-ts"><span class="v" id="f1-correct"></span> of them correct</span></div>
<div class="ex-tile"><span class="ex-tl">Precision</span><span class="ex-tv" id="f1-p"></span>
<span class="ex-ts">correct ÷ kept: how many of the kept calls are right</span></div>
<div class="ex-tile"><span class="ex-tl">Recall</span><span class="ex-tv" id="f1-r"></span>
<span class="ex-ts">domains found ÷ true domains = <span class="v" id="f1-found"></span> ÷ <span class="v" id="f1-ntrue"></span></span></div>
<div class="ex-tile now"><span class="ex-tl">F at your cutoff</span><span class="ex-tv" id="f1-f"></span>
<span class="ex-ts">the balance of precision and recall</span></div>
<div class="ex-tile best"><span class="ex-tl">○ Highest F (Fmax)</span><span class="ex-tv" id="f1-fmax"></span>
<span class="ex-ts">at any cutoff in <span class="v" id="f1-best"></span>; it does not move with the slider</span></div>
</div>
<div class="ex-box" id="f1-track"></div>
<div class="ex-box" id="f1-curve"></div>
</figure>
<details class="prose ex-more"><summary>How the leaderboard does this over many proteins</summary>
<p>The leaderboard has many proteins, not one. At each cutoff it works out precision and recall
for each protein, then averages them over proteins; F comes from the two averages. This is the
protein-by-protein Fmax of the CAFA function-prediction challenge. Precision is averaged over the
proteins with at least one call kept; recall over every protein with a true annotation. The
<b>true instances</b> (n) on the leaderboard are those true annotations, the denominator of
recall. The leaderboard tries 101 cutoffs, evenly spaced through the tool's own scores; here
every whole-number cutoff from 0 to 100 is tried.</p>
</details>

<h2 id="auprc">AUPRC: every cutoff at once</h2>
<div class="prose">
<p>Fmax keeps one cutoff. <b>AUPRC</b>, the area under the precision-recall curve, looks at every
cutoff. Put every call in one list, best score first. Take the calls one at a time; after each,
write down the precision and recall so far. Each call that finds a new domain adds a strip of
area: its width is the recall gained and its height is the precision at that point.</p>
<p>\[\text{AUPRC} = \sum_k (R_k - R_{k-1})\,P_k\]</p>
<p>Here \(P_k\) and \(R_k\) are the precision and recall after the \(k\)-th call of the list.
A wrong call adds no area itself, but it lowers the precision, and so the height of every strip
after it. AUPRC runs from 0 to 1, and higher is better.</p>
</div>
<figure class="ex-fig" id="fig2">
<figcaption><span class="ex-made-up">made-up example</span><b>Figure 5.</b> The precision-recall curve
of protein A's calls, and what happens when protein B is added. <details class="capmore"><summary>How to read Figure 5</summary>Each point is one call, taken in
score order; its position is the precision and recall after taking it. The shaded strips add up to
the AUPRC. Lower the cutoff to add calls one by one, and switch protein B on to see Fmax and
AUPRC move in opposite directions.</details></figcaption>
<div class="ex-legend" aria-label="Key to Figure 5">
@@LEG2@@
</div>
<div class="ex-controls">
<label for="f2-cut"><span>Score cutoff: take calls scoring at least <span class="v mono" id="f2-t"></span></span>
<input type="range" id="f2-cut" min="0" max="100" step="1" value="0"></label>
<span class="chipgroup" id="f2-which" role="group" aria-label="Proteins in the list">
<button type="button" class="chip" data-b="0" aria-pressed="true">protein A only</button>
<button type="button" class="chip" data-b="1" aria-pressed="false">add protein B</button></span>
</div>
<div class="ex-box" id="f2-curve"></div>
<p class="ex-read">Area of the calls taken so far: <span class="v" id="f2-area"></span>. AUPRC, with every
call taken: <span class="v" id="f2-auprc"></span>. Fmax of the same proteins: <span class="v" id="f2-fmax"></span>.</p>
<div class="ex-scroll"><table class="ex-t" id="f2-table"><thead><tr><th>score</th><th data-col="a">protein A only</th><th data-col="ab">A and B</th></tr></thead>
<tbody><tr><td>Fmax</td><td class="mono" data-col="a" id="f2-a-only-fmax"></td><td class="mono" data-col="ab" id="f2-ab-fmax"></td></tr>
<tr><td>AUPRC</td><td class="mono" data-col="a" id="f2-a-only-auprc"></td><td class="mono" data-col="ab" id="f2-ab-auprc"></td></tr></tbody></table></div>
<p class="ex-read">The list behind the curve, best score first. Grey rows score below the cutoff.</p>
<div class="ex-scroll"><table class="ex-t" id="f2-list"><thead><tr><th>rank</th><th>protein</th><th>call</th><th></th>
<th>precision</th><th>recall</th><th>area added</th></tr></thead><tbody></tbody></table></div>
</figure>
<details class="prose ex-more"><summary>Why adding protein B moves Fmax up and AUPRC down</summary>
<p>Protein B has five Ig domains, and its five calls are all correct, but they score 45 to 60.
Fmax averages protein by protein, so protein B counts as one protein of two, whatever its number
of domains. On its own protein B is perfect, so Fmax goes up. AUPRC pools every call of every
protein into one list. Protein B's domains are five of the eight true domains. In the list its
calls come after one or two of protein A's wrong calls (scores 70 and 55), so AUPRC goes down.
Use Fmax when
you will keep every call above one cutoff. Use AUPRC when you will rank calls and read down the
list.</p>
</details>

<h2 id="iou">Median IoU: are the ends in the right place?</h2>
<div class="prose">
<p><b>IoU</b> (intersection over union) is the number of residues the call and the true domain
share, divided by the number of residues covered by either. It is 1 when the call starts and ends
where the domain does, and 0 when they do not overlap. A call needs an IoU of 0.5 or more to be correct.</p>
<p>\[\text{IoU} = \frac{|\text{call} \cap \text{true domain}|}{|\text{call} \cup \text{true domain}|}\]</p>
</div>
<figure class="ex-fig" id="fig3">
<figcaption><span class="ex-made-up">made-up example</span><b>Figure 6.</b> IoU of one call with one true
domain. <details class="capmore"><summary>How to read Figure 6</summary>Top: a 200-residue protein with one true domain, and the call below it; the brackets
underneath show the shared residues and the residues covered by either. It starts at the example on
the Methods page (call 30–129, domain 10–109). Drag the call, or its handles, to see IoU change. Bottom: the IoU
of each correct call of protein A in Figure 4, with your call added when it is correct; the dashed
line is their median.</details></figcaption>
<div class="ex-legend" aria-label="Key to Figure 6">
@@LEG3@@
</div>
<p class="ex-hint">Drag the black call to move it, or drag its handles to move one end.</p>
<div class="ex-box ex-drag" id="f3-track"></div>
<div class="ex-tiles" aria-live="polite">
<div class="ex-tile"><span class="ex-tl">Shared residues</span><span class="ex-tv" id="f3-inter"></span>
<span class="ex-ts">in both the call and the true domain</span></div>
<div class="ex-tile"><span class="ex-tl">Covered by either</span><span class="ex-tv" id="f3-union"></span>
<span class="ex-ts">in the call, the true domain, or both</span></div>
<div class="ex-tile now"><span class="ex-tl">IoU of your call</span><span class="ex-tv" id="f3-iou"></span>
<span class="ex-ts">shared ÷ covered by either: <span id="f3-verdict"></span></span></div>
<div class="ex-tile best"><span class="ex-tl">Median IoU</span><span class="ex-tv" id="f3-median"></span>
<span class="ex-ts">the middle of the <span class="v" id="f3-n"></span> correct calls below</span></div>
</div>
<div class="ex-controls ex-exact">
<label for="f3-start"><span>Call starts at residue <span class="v mono" id="f3-start-v"></span></span>
<input type="range" id="f3-start" min="1" max="200" step="1" value="@@S@@"></label>
<label for="f3-end"><span>Call ends at residue <span class="v mono" id="f3-end-v"></span></span>
<input type="range" id="f3-end" min="1" max="200" step="1" value="@@E@@"></label>
</div>
<div class="ex-box" id="f3-strip"></div>
<p class="ex-read"><span id="f3-mine"></span></p>
</figure>
<details class="prose ex-more"><summary>Median IoU on the leaderboard</summary>
<p><b>Median IoU</b> on the leaderboard is the middle IoU of all of a tool's correct calls. Half of
the correct calls match their domain at least that well. Since a correct call has an IoU of at least
0.5, the median runs from 0.5 to 1. Single-residue annotations, such as an active site, are left
out: an IoU against one residue says nothing about where the ends are.</p>
</details>

<h2 id="dedup">Raw or dedup: one domain, many calls</h2>
<div class="prose">
<p>A tool often finds the same domain through several proteins of the relative, and reports one
call for each. <b>Raw</b> keeps every call: the first copy can be correct, and every other copy
counts as wrong. <b>Dedup</b> merges calls of the same family on the same protein that overlap each
other with an IoU of 0.5 or more, keeping the highest-scoring one. Two copies of a domain side by
side barely overlap, so both are kept.</p>
</div>
<figure class="ex-fig" id="fig4">
<figcaption><span class="ex-made-up">made-up example</span><b>Figure 7.</b> Protein A's seven calls, raw
and after dedup. <details class="capmore"><summary>How to read Figure 7</summary>Switch between the two: dedup merges the two lower-scoring PH calls into the PH
call scoring 95, and the table shows how the scores change.</details></figcaption>
<div class="ex-legend" aria-label="Key to Figure 7">
@@LEG4@@
</div>
<div class="ex-controls">
<span class="chipgroup" id="f4-which" role="group" aria-label="Raw or dedup">
<button type="button" class="chip" data-mode="raw" aria-pressed="true">raw</button>
<button type="button" class="chip" data-mode="dedup" aria-pressed="false">dedup</button></span>
</div>
<div class="ex-box" id="f4-track"></div>
<div class="ex-scroll"><table class="ex-t" id="f4-table"><thead><tr><th>score</th><th data-col="raw">raw</th><th data-col="dedup">dedup</th></tr></thead>
<tbody><tr><td>calls</td><td class="mono" data-col="raw" id="f4-raw-n"></td><td class="mono" data-col="dedup" id="f4-dedup-n"></td></tr>
<tr><td>precision, every call kept</td><td class="mono" data-col="raw" id="f4-raw-p0"></td><td class="mono" data-col="dedup" id="f4-dedup-p0"></td></tr>
<tr><td>Fmax</td><td class="mono" data-col="raw" id="f4-raw-fmax"></td><td class="mono" data-col="dedup" id="f4-dedup-fmax"></td></tr>
<tr><td>AUPRC</td><td class="mono" data-col="raw" id="f4-raw-auprc"></td><td class="mono" data-col="dedup" id="f4-dedup-auprc"></td></tr></tbody></table></div>
</figure>
<details class="prose ex-more"><summary>When dedup raises a score, and when it lowers one</summary>
<p>Here both merged calls were wrong copies, so recall does not change, and precision, Fmax and
AUPRC go up. Dedup can also lower a score: when the higher-scoring call misses the domain's ends
and a lower-scoring copy matches them, dedup keeps the one that misses. Comparing raw with dedup
separates finding a domain from reporting it many times.</p>
</details>

<h2 id="cost">Cost: core-hours, GPU-hours and the price of an hour</h2>
<div class="prose">
<p><b>Core-hours</b> are the cores a search reserved multiplied by its wall time, so a task that sat
partly idle still counts all its cores. They come from the run's own task records on Stanford's
Sherlock cluster, so they are marked <span class="prov measured">measured</span>. Numbers a tool's
submitter gives, on hardware they name, would be marked <span class="prov reported">reported</span>
(CPU-h reported); none have been submitted yet. The leaderboard divides the search's core-hours by
the number of relatives searched and scales them to 1,000 human proteins (queries). <b>GPU-hours</b>
are the same for GPU time; no search here used a GPU. <b>Peak RAM</b> is the largest memory any of
the search's tasks held, in GB; it is not priced.</p>
<p>The estimated cost, cost_usd_est, puts a price on each hour:</p>
<p>\[\begin{aligned}\text{cost\_usd\_est} &= \text{core-hours} \times r_\text{cpu}\\ &+ \text{GPU-hours} \times r_\text{gpu}\end{aligned}\]</p>
<p>per 1,000 queries, where \(r_\text{cpu}\) is the price of one CPU-hour and \(r_\text{gpu}\) of one
GPU-hour, in US dollars. A core-hour is priced as one CPU-hour.</p>
</div>
<figure class="ex-fig" id="fig5">
<figcaption><b>Figure 8.</b> The cost formula with your own numbers. <details class="capmore"><summary>How to read Figure 8</summary>It starts at phmmer's measured
core-hours per 1,000 queries against zebrafish (from Table 1) and the site's default prices. These
are AWS on-demand rates for a c7i.large CPU (half its 2-vCPU hourly price) and a g5.xlarge GPU. Drag any underlined
number to see the cost change.</details></figcaption>
<p class="ex-hint">Drag any <span class="scrub-demo">underlined number</span> left or right to change it
(hold Shift for small steps), or click it and use the arrow keys.</p>
<p class="ex-formula">cost per 1,000 queries =
<span class="scrub" data-for="f5-cpu" data-step="0.005" data-name="core-hours per 1,000 queries"></span> core-hours
× $<span class="scrub" data-for="f5-rcpu" data-step="0.0005" data-name="price per CPU-hour"></span> per CPU-hour
+ <span class="scrub" data-for="f5-gpu" data-step="0.01" data-name="GPU-hours per 1,000 queries"></span> GPU-hours
× $<span class="scrub" data-for="f5-rgpu" data-step="0.01" data-name="price per GPU-hour"></span> per GPU-hour</p>
<div class="ex-tiles" aria-live="polite">
<div class="ex-tile"><span class="ex-tl">CPU part</span><span class="ex-tv" id="f5-a"></span>
<span class="ex-ts">core-hours × price per CPU-hour</span></div>
<div class="ex-tile"><span class="ex-tl">GPU part</span><span class="ex-tv" id="f5-b"></span>
<span class="ex-ts">GPU-hours × price per GPU-hour</span></div>
<div class="ex-tile best"><span class="ex-tl">Estimated cost</span><span class="ex-tv" id="f5-out"></span>
<span class="ex-ts">US dollars per 1,000 queries</span></div>
</div>
<details class="ex-exact"><summary>Type exact numbers</summary>
<div class="ex-controls">
<label for="f5-cpu">core-hours per 1,000 queries <input type="number" id="f5-cpu" min="0" step="any" inputmode="decimal" value="@@CPU@@"></label>
<label for="f5-gpu">GPU-hours per 1,000 queries <input type="number" id="f5-gpu" min="0" step="any" inputmode="decimal" value="0"></label>
<label for="f5-rcpu">r_cpu, USD per CPU-hour <input type="number" id="f5-rcpu" min="0" step="any" inputmode="decimal" value="@@RCPU@@"></label>
<label for="f5-rgpu">r_gpu, USD per GPU-hour <input type="number" id="f5-rgpu" min="0" step="any" inputmode="decimal" value="@@RGPU@@"></label>
</div></details>
</figure>
<p class="prose muted">The made-up example is also in this page's JSON copy. Definitions of every
score as the benchmark computes them, with a worked example from the real data, are on
<a href="methods.html">Methods</a>.</p>
<script type="application/json" id="toy-data">@@TOY@@</script>
<script defer src="assets/explainers.js"></script>
<link rel="stylesheet" href="assets/katex/katex.min.css">
<script defer src="assets/katex/katex.min.js"></script>
<script defer src="assets/katex/contrib/auto-render.min.js"
 onload="renderMathInElement(document.body,{delimiters:[{left:'\\[',right:'\\]',display:true},{left:'\\(',right:'\\)',display:false}]})"></script>
"""
    leg1 = "".join([li("true", "true domain on the protein line"), li("kept", "call kept at this cutoff"),
                    li("drop", "call below the cutoff: not counted"), glyphs, li("line", "F at each cutoff"),
                    li("cut", "your cutoff"), li("ring", "highest F (Fmax)")])
    leg2 = "".join([li("dot", "call taken: precision and recall after it"), li("dot-off", "call below the cutoff"),
                    li("area", "area of the calls taken"), li("area-off", "area the lower-scoring calls would add")])
    leg3 = "".join([li("true", "true domain"), li("kept", "the call you set"),
                    li("bracket", "shared residues, and residues covered by either"),
                    li("dot", "IoU of a correct call of protein A"), li("your", "IoU of your call"),
                    li("median", "median")])
    leg4 = "".join([li("true", "true domain"), li("kept", "call that counts"),
                    li("drop", "call merged into a higher-scoring one: not counted"), glyphs])
    for k, v in {"@@LEG1@@": leg1, "@@LEG2@@": leg2, "@@LEG3@@": leg3, "@@LEG4@@": leg4,
                 "@@S@@": str(TOY["iou"]["call"]["start"]), "@@E@@": str(TOY["iou"]["call"]["end"]),
                 "@@CPU@@": f"{per1k:.4f}", "@@RCPU@@": f"{R_CPU:.6g}", "@@RGPU@@": f"{R_GPU}",
                 "@@TOY@@": toy_json}.items():
        body = body.replace(k, v)
    return layout("How the scores work", "explainers", body, src,
                  description="Fmax, AUPRC, median IoU, raw against dedup and the cost formula, each worked "
                              "through on a made-up protein you can change.")


def page_claims_schema(src: dict, cards: list[dict], inv: dict) -> str:
    example = (ROOT_DIR / "claims" / "kmerseek.yaml").read_text().split("claims:\n", 1)[1].split("\n\n  - id:", 1)[0]
    badges = "".join(f'<p><img src="badges/claims-{c["tool"]}.svg" alt="claims for {esc(c["tool"])}: '
                     f'{c["counts"]["alive"]} alive, {c["counts"]["retired"]} retired, {c["counts"]["untested"]} untested"></p>'
                     f'<pre><code>![claims](https://seanome.github.io/protein-function-annotation-leaderboard/badges/claims-{c["tool"]}.svg)</code></pre>'
                     for c in cards)
    body = f"""<h1>claims.yaml</h1>
<div class="prose">
<p>One file per tool, <span class="mono">claims/&lt;tool&gt;.yaml</span>, lists the tool's claims. Each claim
has the claim in plain words, the result that would retire it (<b>kill</b>), the data it is tested on
(<b>scope</b>), and a <b>test</b>. The file never states a verdict. <span class="mono">scripts/claims.py</span>
runs each test and sets the status:</p>
<ul>
<li><b>table</b>, <b>table_bins</b>, <b>table_max</b>, <b>table_count_passing</b>: recomputed from
<span class="mono">results.parquet</span> on every build.</li>
<li><b>record</b>, <b>record_compare</b>: a number quoted from the experiment inventory
(<span class="mono">data/experiments.json</span>, copied from {esc(inv['repo'])} at commit
<a class="mono" href="https://github.com/{esc(inv['repo'])}/blob/{esc(inv['commit'])}/{esc(inv['path'])}">{esc(inv['commit'][:8])}</a>,
PR {inv['pr']}; recorded in <span class="mono">data/experiments.SOURCE.json</span>). The quote must
match the recorded value character for character, so a notebook that changes a number fails the check until the claim
is updated.</li>
<li><b>not_run</b>: no test yet; the claim is untested and the reason is shown.</li>
</ul>
<p>A test may carry a <b>verdict</b>: the wording of its verdict line, with the resolved numbers
filled in. Without one, the line is the test's own result.</p>
<p>A claim retired by a record test is dated with the date of that inventory record. A claim retired
by a table test has no date, because the tables are recomputed on every build and the build cannot
tell when the claim was first disproved. The "last disproved" date on a card counts only dated
results.</p>
<p>The full schema: <a class="mono" href="claims/schema.json">claims/schema.json</a>. The first claim of
<a class="mono" href="claims/kmerseek.yaml">claims/kmerseek.yaml</a>:</p>
<pre><code>{esc(example)}</code></pre>
<h2 id="badge">Badge and CI check</h2>
<p>Every push and pull request runs <span class="mono">python scripts/claims.py --check</span> in a GitHub
Action (<a class="mono" href="https://github.com/{REPO}/blob/main/.github/workflows/claims.yml">.github/workflows/claims.yml</a>).
It fails when a file breaks the schema, a quote no longer matches its record, or a table test cannot
find its rows. A second job rebuilds the site from <span class="mono">data/</span> and fails if any
page differs from the committed one. The build also writes a badge per tool from the same verdicts:</p>
{badges}
</div>"""
    return layout("claims.yaml", "claims", body, src, description="The claims.yaml format, the CI check and the badge.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data"))
    ap.add_argument("--out", type=Path, default=Path("docs"))
    a = ap.parse_args()
    global DATA_DIR
    DATA_DIR = a.data
    res, comp, src, prov, cov = load(a.data)

    unknown = set(res["tool"].unique()) - set(TOOLS)
    if unknown:
        raise SystemExit(f"tools with no card facts: {unknown}")
    species_mya = {r["species"]: r["species_mya"] for r in
                   res.select("species", "species_mya").unique().iter_rows(named=True)
                   if r["species_mya"] is not None and not math.isnan(r["species_mya"])}
    global SPECIES_MYA, DEFAULT_WORDS
    SPECIES_MYA = species_mya
    DEFAULT_WORDS = (f"Pfam truth set, raw calls, zebrafish ({mya_label(species_mya['zebrafish'])} Mya) "
                     "as the relative searched")
    tools_used = set(res["tool"].unique())
    rows = row_order(res)
    searched = sorted({f'{r["tool"]}|{r["setting"]}|{r["species"]}'
                       for r in res.select("tool", "setting", "species").unique().iter_rows(named=True)})
    meta = {
        "defaults": DEFAULTS,
        "chiprows": mark_empty(chiprows(species_mya), res),
        "rows": rows,
        "tools": {t: {"name": TOOLS[t]["name"], "kind": TOOLS[t]["kind"], "shape": TOOLS[t]["shape"],
                      "tags": [w for k, w in (("3D structure prediction", "needs 3D structure prediction"), ("GPU", "needs GPU"))
                               if TOOLS[t]["needs"][k]]}
                  for t in sorted(tools_used)},
        "kinds": KINDS,
        "species": {s: {"name": SPECIES[s], "mya": mya_label(species_mya[s]) if s in species_mya else None}
                    for s in SPECIES},
        "truthNames": TRUTH,
        "metrics": METRICS,
        "axes": AXES,
        "searched": searched,
        "source": src,
        # "?" buttons with their definitions, drawn by app.js into Table 1's header and the
        # filter rows (terms.js shows them).
        "terms": {"th": {t: term_btn(t, "th") for t in TERMS}, "short": {t: PLAIN[t][0] for t in PLAIN}, "filters": {"calls": term_btn("calls", "filters")},
                  "axes": {a["id"]: term_btn(a["id"], "ax") for a in AXES if a["id"] in TERMS}},
    }
    payload = json.dumps({"meta": meta, "results": results_json(res), "compute": compute_json(comp)},
                         separators=(",", ":"), allow_nan=False)
    payload = payload.replace("</", "<\\/")
    n_q = prov["n_queries"]

    out = a.out
    (out / "tools").mkdir(parents=True, exist_ok=True)
    (out / "data").mkdir(parents=True, exist_ok=True)
    (out / ".nojekyll").write_text("")
    phm = res.filter((pl.col("tool") == "hmmer3_phmmer") & (pl.col("truth_set") == "pfam")
                     & (pl.col("species") == "zebrafish") & (pl.col("stratum_axis") == "all"))
    dd = {f"{k}_{m}": phm.filter(pl.col("calls") == k)[m][0] for k in ("raw", "dedup") for m in ("fmax", "auprc")}
    (out / "leaderboard.html").write_text(page_leaderboard(payload, src, n_q, tools_used, dd))
    (out / "cost.html").write_text(page_cost(payload, src, n_q, tools_used))
    (out / "coverage.html").write_text(page_coverage(cov, src, species_mya))
    (out / "methods.html").write_text(page_methods(res, comp, src, prov))
    (out / "explainers.html").write_text(page_explainers(src, comp))
    import claims as claims_mod
    schema = json.loads((ROOT_DIR / "claims/schema.json").read_text())
    resolver = claims_mod.Resolver(a.data)
    cards = []
    (out / "claims").mkdir(parents=True, exist_ok=True)
    for f in sorted((ROOT_DIR / "claims").glob("*.yaml")):
        cl = claims_mod.resolve_file(f, resolver, schema)
        cards.append(cl)
        (out / "claims" / f.name).write_text(f.read_text())
        (out / "data/claims").mkdir(parents=True, exist_ok=True)
        (out / f"data/claims/{cl['tool']}.json").write_text(json.dumps(cl, indent=1, ensure_ascii=False) + "\n")
        (out / "badges").mkdir(parents=True, exist_ok=True)
        (out / f"badges/claims-{cl['tool']}.svg").write_text(claims_mod.badge(cl["tool"], cl["counts"]))
    (out / "claims" / "schema.json").write_text((ROOT_DIR / "claims/schema.json").read_text())
    claims_by_tool = {c["tool"]: c for c in cards}
    (out / "claims.html").write_text(page_claims(src, cards))
    (out / "tools" / "index.html").write_text(page_tools_index(src, tools_used, res, comp, claims_by_tool))
    (out / "claims-schema.html").write_text(page_claims_schema(src, cards, resolver.inv_src))
    for t in sorted(tools_used):
        nums = card_numbers(t, res, comp)
        md = card_markdown(t, nums, src)
        (out / "tools" / f"{t}.md").write_text(md[0] + "\n" + md[1])
        card = {
            "tool": t, "name": TOOLS[t]["name"], "kind": KINDS[TOOLS[t]["kind"]]["label"],
            "compares": TOOLS[t]["compares"], "needs": TOOLS[t]["needs"],
            "reference": TOOLS[t]["reference"], "score": TOOLS[t]["score"],
            "failure_modes": TOOLS[t]["failure_modes"], "version": TOOLS[t]["version"],
            "licence": TOOLS[t]["licence"], "commercial_use": TOOLS[t]["commercial"],
            "licence_url": TOOLS[t]["licence_url"],
            "container": TOOLS[t]["container"], "command": TOOLS[t]["command"],
            "numbers_from": {"repo": src["repo"], "commit": src["commit"],
                             "filters": "pfam, raw, " + nums["species"]},
            "best_setting": nums.get("best_setting"),
            "strong": nums.get("strong", []), "weak": nums.get("weak", []),
            "cost_per_1000_queries": None if not nums.get("cost_here") else {
                "core_h_reserved_per_relative": nums["cost_here"]["cpu_h_per_1000_queries"],
                "n_relatives": nums["cost_here"]["n_relatives"],
                "database_build_core_h": nums["cost_here"]["db_cpu_h"],
                "database_build_gpu_h": nums["cost_here"]["db_gpu_h"],
                "peak_rss_gb": nums["cost_here"]["peak_rss_gb"],
                "hardware": nums["cost_here"]["hardware"]},
            "questions": [re.sub(r"\*\*", "", q[3:]) for q in questions(t, nums)],
            "markdown": f"{t}.md",
        }
        (out / "tools" / f"{t}.json").write_text(json.dumps(card, indent=1, allow_nan=False) + "\n")
        (out / "tools" / f"{t}.html").write_text(page_card(t, md, src, card_figure(t, nums, res), claims_by_tool.get(t)))
    (out / "data" / "SOURCE.json").write_text(json.dumps(src, indent=1) + "\n")
    if (a.data / "pileup/regions.parquet").exists():
        # The page credits its tables to this commit; without it, the copy is of unknown origin.
        if not src.get("pileup_commit"):
            raise SystemExit("data/pileup/ exists but data/SOURCE.json has no pileup_commit; "
                             "re-run `make tables` so the pile-up tables are credited to their commit")
        import pileup
        calls = pl.read_parquet(a.data / "pileup/calls.parquet").select("tool", "setting").unique()
        kinds = list(KINDS)
        panels = sorted(({"tool": t, "setting": st, "kind": TOOLS[t]["kind"]} for t, st in calls.iter_rows()),
                        key=lambda q: (kinds.index(q["kind"]), list(TOOLS).index(q["tool"]),
                                       setting_key(q["tool"], q["setting"])))
        data, check = pileup.build(a.data, panels)
        (out / "data" / "pileup.json").write_text(json.dumps(data, separators=(",", ":"), allow_nan=False))
        (out / "data" / "pileup_check.json").write_text(json.dumps(check, indent=1) + "\n")
        prov["pileup"] = json.loads((a.data / "pileup/provenance.json").read_text())
        (out / "pileup.html").write_text(page_pileup(src, res, prov, data, check, pileup.readouts(data)))
    else:
        # No pile-up tables in data/: remove the pages an earlier build wrote from older ones.
        for stale in ("pileup.html", "pileup.html.json", "data/pileup.json", "data/pileup_check.json"):
            (out / stale).unlink(missing_ok=True)

    # Every page also as JSON at its path plus .json.
    full = json.loads(payload.replace("<\\/", "</"))
    page_json(out, "leaderboard.html", src, {"note": "The table is drawn from these rows: results keyed "
              "truth|calls|species|stratum_axis|stratum, then tool|setting -> [fmax, auprc, median_iou_tp, "
              "n_truth_instances]; compute keyed tool|setting|species.", **full,
              "prices": {"r_cpu": R_CPU, "r_gpu": R_GPU, "source": PRICE_PAGE, "price_file": PRICE_FILE}})
    page_json(out, "cost.html", src, {"note": "Same rows as leaderboard.html.json. Figure 1 plots any two of "
              "the columns in 'axes' for every tool and setting that has both: scores from results "
              "[fmax, auprc, median_iou_tp], cpu = cpu_h_per_1000, cost = cpu_h_per_1000 * r_cpu + "
              "gpu_h_per_1000 * r_gpu, wall = wall_h_per_1000, ram = peak_rss_gb, db = db_cpu_h_per_db.",
              "axes": AXES, "prices": {"r_cpu": R_CPU, "r_gpu": R_GPU}, "results": full["results"],
              "compute": full["compute"], "rows": full["meta"]["rows"]})
    page_json(out, "coverage.html", src, {"tables": [
        {"columns": h[1:], "rows": [{"tool": r[0], "instances": {c: (int(v.replace(",", "")) if v else None)
                                                             for c, v in zip(h[1:], r[1:])}} for r in rows]}
        for h, rows in parse_tables(cov)[:2]]})
    if (a.data / "inventory.md").exists():
        (out / "data" / "inventory.md").write_text((a.data / "inventory.md").read_text())
    else:
        (out / "data" / "inventory.md").unlink(missing_ok=True)
    page_json(out, "methods.html", src, None)
    page_json(out, "explainers.html", src, {"made_up_example": TOY, "terms": {
        k: {"name": n, "in_plain_words": PLAIN[k][1], "short": PLAIN[k][0], "explainer": h, "definition": w}
        for k, (n, h, w) in TERMS.items()}})
    page_json(out, "tools/index.html", src, {"cards": [f"{t}.json" for t in sorted(tools_used)]})
    page_json(out, "claims.html", src, {"claims": cards})
    page_json(out, "claims-schema.html", src, {"schema": json.loads((ROOT_DIR / "claims/schema.json").read_text())})
    if (out / "pileup.html").exists():
        page_json(out, "pileup.html", src, {
            "note": "Readouts computed by scripts/pileup.py readouts(), the same rules as the page's "
                    "readouts(); the matrix itself is data/pileup.json.",
            "matrix": "data/pileup.json", "per_residue_recall_check": check, **pileup.readouts(data)})
    for t in sorted(tools_used):
        page_json(out, f"tools/{t}.html", src, json.loads((out / "tools" / f"{t}.json").read_text()))
    # The organism pages: the home page and one page per pre-built species.
    if (a.data / "organism/lookups.jsonl").exists():
        sys.path.insert(0, str(ROOT_DIR / "scripts/organism"))
        import build as organism_build
        n_org = organism_build.build(a.data, out, sys.modules[__name__])
        print(f"organism pages: {n_org} pre-built species")
    else:
        raise SystemExit("data/organism/lookups.jsonl is missing: run scripts/organism/fetch.py and "
                         "lookups.py first (see the README), or the home page would not be written")
    print(f"built {out}: {len(rows)} tool-setting rows, {len(tools_used)} cards, payload {len(payload) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
