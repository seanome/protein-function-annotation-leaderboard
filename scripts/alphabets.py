"""The alphabet page (alphabets.html): what a reduced amino acid alphabet keeps, measured on
aligned pairs with the alignment given.

Reads data/alphabets/, copied by scripts/copy_tables.py from results/leaderboard/alphabets/
of one commit of seanome/2024-kmerseek-analysis (written there by
`make leaderboard-alphabets`):

* classes.parquet: residue classes of the 19 kmerseek alphabets
* hp_borderline_residues.csv: PR 116's table of the seven H/P alphabets, 140 squares
* evidence.parquet: notebook 291's evidence per aligned position, ratio to the 20 amino acids
* kappa.parquet: notebook 230's agreement, chance and kappa per identity bin
* usable_k.parquet, k_reach.parquet: notebook 230's usable seed lengths k

Nothing here computes a result: every number is a value from those tables, rounded for
display. The page is drawn on the server; a small script only switches the identity bin of
the usable-k bars.
"""

from __future__ import annotations

import csv
import html
from pathlib import Path

import polars as pl

ANALYSIS = "https://github.com/seanome/2024-kmerseek-analysis"
PR = {n: f'<a href="{ANALYSIS}/pull/{n}">PR {n}</a>' for n in (99, 100, 101, 116)}
NB = {
    230: f'<a href="{ANALYSIS}/blob/main/notebooks/230_hp_class_conservation_in_aligned_homologs.ipynb">notebook 230</a>',
    291: f'<a href="{ANALYSIS}/blob/main/notebooks/291_evidence_per_position_19_alphabets_pfam_seed_scope40.ipynb">notebook 291</a>',
    274: "notebook 274",
}
MIN_N = 30  # fewer than this many pairs: shown grey

# Alphabet rows, grouped by size (as in notebook 230), with a gap between groups.
SIZE_GROUPS = [
    ("2 or 3 classes", ["hp_thomas_dill2", "hp_pbotc_1st_ed2", "hp_thomas_dill_no_c2", "hp_lehninger2",
                        "hp_lehninger_c_nonpolar2", "hp_kyte_doolittle2", "hp_lehninger_hpc3"]),
    ("4 to 8 classes", ["gbmr4", "polarity4", "wwmj5", "dayhoff6", "gbmr7", "funcgroups8"]),
    ("12 to 18 classes", ["mmseqs12", "sdm12", "wass14", "hsdm17", "uniprot18"]),
    ("20 classes", ["protein20"]),
]
# Residue order of PR 116's figure: hydrophobic in all seven H/P alphabets, the five they
# differ on, polar in all seven.
RES_ORDER = "AFILMVCGPWYDEHKNQRST"
HP_WORDS = {"h": "hydrophobic", "p": "polar", "c": "cysteine in its own class"}
# Evidence columns of the main table: the fair test, homologs only.
EV_GROUPS = [("Pfam seed", "Pfam seed families"),
             ("SCOPe: same family", "SCOPe, same family"),
             ("SCOPe: same superfamily, different family", "SCOPe, same superfamily, different family")]
FAIR = "low_half"
IDENT = ["<20%", "20-30%", "30-40%", "40-60%", ">=60%"]
IDENT_WORDS = {"<20%": "20% or less", "20-30%": "20–30%", "30-40%": "30–40%", "40-60%": "40–60%",
               ">=60%": "over 60%"}
K_BINS = ["20-30%", "30-40%", "40-60%"]
K_AXIS = (3, 40)


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def load(data: Path) -> dict:
    d = data / "alphabets"
    t = {n: pl.read_parquet(d / f"{n}.parquet") for n in ("classes", "evidence", "kappa", "usable_k", "k_reach")}
    with open(d / "hp_borderline_residues.csv") as f:
        t["hp116"] = list(csv.DictReader(f))
    order = [a for _, g in SIZE_GROUPS for a in g]
    have = set(t["classes"]["alphabet"].unique())
    if sorted(order) != sorted(have):
        raise SystemExit(f"alphabets in data/alphabets differ from the page's list: {sorted(set(order) ^ have)}")
    t["order"] = order
    return t


def f3(x: float) -> str:
    return f"{x:.3f}"


def classes_of(t: dict, a: str) -> list[dict]:
    """The alphabet's classes in reading order: residues by RES_ORDER, classes by first residue."""
    rows = t["classes"].filter(pl.col("alphabet") == a)
    out = {}
    for r in rows.iter_rows(named=True):
        out.setdefault(r["class_members"], {"members": r["class_members"], "hp": r["hp_class"]})
    cls = list(out.values())
    for c in cls:
        c["letters"] = sorted(c["members"], key=RES_ORDER.index)
    return sorted(cls, key=lambda c: RES_ORDER.index(c["letters"][0]))


def class_cell(t: dict, a: str) -> str:
    blocks = []
    for c in classes_of(t, a):
        cls = f' hp-{c["hp"]}' if c["hp"] else ""
        squares = "".join(f'<span class="sq">{r}</span>' for r in c["letters"])
        word = f', {HP_WORDS[c["hp"]]}' if c["hp"] else ""
        blocks.append(f'<span class="cls{cls}" data-members="{"".join(c["letters"])}" '
                      f'title="class {"".join(c["letters"])}{word}">{squares}</span>')
    return "".join(blocks)


def ev_scale(ev: pl.DataFrame) -> tuple[float, float]:
    lo, hi = min(ev["lo"].min(), 1.0), max(ev["hi"].max(), 1.0)
    pad = (hi - lo) * 0.04
    return lo - pad, hi + pad


W_EV, W_K = 128, 200


def ev_x(v: float, dom) -> float:
    return 6 + (W_EV - 12) * (v - dom[0]) / (dom[1] - dom[0])


def ev_axis(dom) -> str:
    # Ticks at 1 (the 20 amino acids) and at the round tenths nearest each end of the scale,
    # dropping an end tick that would sit on top of the tick at 1.
    import math
    span = dom[1] - dom[0]
    ends = [math.ceil(dom[0] * 10) / 10, math.floor(dom[1] * 10) / 10]
    ticks = [1.0] + [x for x in ends if abs(x - 1.0) > 0.18 * span]
    s = "".join(f'<line x1="{ev_x(x, dom):.1f}" x2="{ev_x(x, dom):.1f}" y1="0" y2="5" stroke="var(--muted)"/>'
                f'<text x="{ev_x(x, dom):.1f}" y="16" text-anchor="middle">{x:g}</text>' for x in sorted(set(ticks)))
    return f'<svg class="axis" width="{W_EV}" height="18" aria-hidden="true">{s}</svg>'


def ev_cell(r: dict, dom, n_txt: str) -> str:
    grey = r["n_items"] < MIN_N
    x, xl, xh, x1 = (ev_x(v, dom) for v in (r["ratio_to_20aa"], r["lo"], r["hi"], 1.0))
    txt = f'{f3(r["ratio_to_20aa"])} ({f3(r["lo"])}–{f3(r["hi"])})'
    tip = (f'{r["alphabet"]}: {txt} times the evidence per aligned position of the 20 amino acids, '
           f'95% interval in brackets. Rank {r["rank"]} of 19. {n_txt}')
    svg = (f'<svg width="{W_EV}" height="16" aria-hidden="true">'
           f'<line x1="{x1:.1f}" x2="{x1:.1f}" y1="0" y2="16" stroke="var(--muted)" stroke-dasharray="3 2"/>'
           f'<line x1="{xl:.1f}" x2="{xh:.1f}" y1="8" y2="8" stroke="var(--ink)" stroke-width="1.6"/>'
           f'<circle cx="{x:.1f}" cy="8" r="3.4" fill="var(--ink)"/></svg>')
    return (f'<td class="ev{" grey" if grey else ""}" data-alphabet="{esc(r["alphabet"])}" data-group="{esc(r["group"])}" '
            f'data-ratio="{r["ratio_to_20aa"]!r}" data-lo="{r["lo"]!r}" data-hi="{r["hi"]!r}" '
            f'data-n="{r["n_items"]}" data-clusters="{r["n_clusters"]}" title="{esc(tip)}">{svg}'
            f'<span class="v">{txt}</span></td>')


def kappa_cells(t: dict, a: str, kmin: float, kmax: float) -> str:
    out = []
    k = t["kappa"].filter(pl.col("alphabet") == a)
    for b in IDENT:
        r = k.filter(pl.col("identity_bin") == b).row(0, named=True)
        pct = 12 + 70 * (r["kappa"] - kmin) / (kmax - kmin)
        ink = "var(--bg)" if pct > 50 else "var(--ink)"
        grey = r["n_pairs"] < MIN_N
        tip = (f'{a}, Pfam seed pairs at {IDENT_WORDS[b]} identity: kappa {r["kappa"]:.2f}. '
               f'Same class at {100 * r["agreement"]:.0f}% of aligned positions, {100 * r["chance"]:.0f}% expected by chance. '
               f'n = {r["n_pairs"]:,} pairs.')
        bg = "var(--null)" if grey else f"color-mix(in srgb, var(--accent) {pct:.0f}%, var(--surface))"
        out.append(f'<td class="kp{" grey" if grey else ""}" data-alphabet="{esc(a)}" data-bin="{esc(b)}" '
                   f'data-kappa="{r["kappa"]!r}" data-agreement="{r["agreement"]!r}" data-chance="{r["chance"]!r}" '
                   f'data-n="{r["n_pairs"]}" style="background:{bg};color:{ink}" title="{esc(tip)}">'
                   f'{r["kappa"]:.2f}</td>')
    return "".join(out)


def k_x(k: float) -> float:
    return 4 + (W_K - 8) * (k - K_AXIS[0]) / (K_AXIS[1] - K_AXIS[0])


def k_axis() -> str:
    s = "".join(f'<line x1="{k_x(k - 0.5):.1f}" x2="{k_x(k - 0.5):.1f}" y1="0" y2="5" stroke="var(--muted)"/>'
                f'<text x="{k_x(k - 0.5):.1f}" y="16" text-anchor="middle">{k}</text>' for k in (5, 10, 15, 20, 25, 30, 35, 40))
    return f'<svg class="axis" width="{W_K}" height="18" aria-hidden="true">{s}</svg>'


def k_cell(t: dict, a: str) -> str:
    u = t["usable_k"].filter(pl.col("alphabet") == a)
    reach = t["k_reach"].filter(pl.col("alphabet") == a)
    parts = []
    for b in K_BINS:
        r = u.filter(pl.col("identity_bin") == b).row(0, named=True)
        lo, hi, most = r["k_use_lo"], r["k_use_hi"], r["k_most_hi"]
        rb = reach.filter(pl.col("identity_bin") == b)
        s_min = rb.filter(pl.col("k") == r["k_min"])["share_reached"][0]
        s_hi = rb.filter(pl.col("k") == hi)["share_reached"][0]
        rng = str(lo) if lo == hi else f"{lo}–{hi}"
        tip = (f'{a}, Pfam seed pairs at {b.replace("-", "–")} identity (n = {r["n_pairs"]:,}): usable k {rng}. '
               f'{100 * s_min:.0f}% of pairs share a same-class run of k = {r["k_min"]} or longer, '
               f'{100 * s_hi:.0f}% one of k = {hi}. Chance-match range from notebook 274: k {r["k_min"]}–{r["k_max"]}.')
        y = 9
        g = (f'<line x1="{k_x(r["k_min"] - 0.5):.1f}" x2="{k_x(r["k_max"] + 0.5):.1f}" y1="{y}" y2="{y}" stroke="var(--muted)" stroke-width="1.2"/>'
             f'<line x1="{k_x(r["k_max"] + 0.5):.1f}" x2="{k_x(r["k_max"] + 0.5):.1f}" y1="{y - 4}" y2="{y + 4}" stroke="var(--muted)" stroke-width="1.2"/>'
             f'<rect x="{k_x(lo - 0.5):.1f}" y="{y - 5}" width="{k_x(hi + 0.5) - k_x(lo - 0.5):.1f}" height="10" fill="var(--k-use-light)"/>')
        if most is not None:
            g += f'<rect x="{k_x(lo - 0.5):.1f}" y="{y - 5}" width="{k_x(most + 0.5) - k_x(lo - 0.5):.1f}" height="10" fill="var(--k-use)"/>'
        hidden = "" if b == K_BINS[0] else " hidden"
        parts.append(f'<span class="kb" data-bin="{esc(b)}" data-k-min="{r["k_min"]}" data-k-max="{r["k_max"]}" '
                     f'data-k-lo="{lo}" data-k-hi="{hi}" data-k-most="{"" if most is None else most}" data-n="{r["n_pairs"]}" '
                     f'title="{esc(tip)}"{hidden}><svg width="{W_K}" height="18" aria-hidden="true">{g}</svg>'
                     f'<span class="v">{rng}</span></span>')
    return f'<td class="uk" data-alphabet="{esc(a)}">{"".join(parts)}</td>'


def hp_figure(t: dict) -> str:
    rows = t["hp116"]
    alphs = list(dict.fromkeys(r["alphabet"] for r in sorted(rows, key=lambda r: int(r["row"]))))
    groups = list(dict.fromkeys(r["residue_group"] for r in rows))
    res_of = {g: list(dict.fromkeys(r["residue"] for r in rows if r["residue_group"] == g)) for g in groups}
    head1 = "<th></th>" + "".join(f'<th colspan="{len(res_of[g])}" class="grp">{esc(g)}</th>' for g in groups) + '<th></th>'
    head2 = '<th class="txt">alphabet</th>' + "".join(
        f'<th class="res{" gstart" if i == 0 else ""}">{r}</th>' for g in groups for i, r in enumerate(res_of[g])) + \
        '<th>residues in h<br>(of 20)</th>'
    body = []
    for a in alphs:
        cells = {r["residue"]: r for r in rows if r["alphabet"] == a}
        tds = "".join(
            f'<td class="sq116 hp-{cells[r]["class"]}{" gstart" if i == 0 else ""}" data-alphabet="{esc(a)}" data-residue="{r}" '
            f'data-class="{cells[r]["class"]}" title="{esc(a)}: {r} is {HP_WORDS[cells[r]["class"]]}">{cells[r]["class"]}</td>'
            for g in groups for i, r in enumerate(res_of[g]))
        n_h = next(iter(cells.values()))["n_residues_in_h"]
        body.append(f'<tr><th class="txt mono" scope="row">{esc(a)}</th>{tds}<td class="num" data-n-h="{n_h}">{n_h}</td></tr>')
    return f'<table class="hp116"><thead><tr>{head1}</tr><tr>{head2}</tr></thead><tbody>{"".join(body)}</tbody></table>'


def build(t: dict, src: dict) -> tuple[str, dict]:
    ev = t["evidence"]
    fair = ev.filter(pl.col("mode") == FAIR)
    kmin, kmax = t["kappa"]["kappa"].min(), t["kappa"]["kappa"].max()
    n_kappa = {r["identity_bin"]: r["n_pairs"] for r in t["kappa"].filter(pl.col("alphabet") == "protein20").iter_rows(named=True)}
    ev_n, ev_dom, ev_head = {}, {}, []
    for g, words in EV_GROUPS:
        e = fair.filter(pl.col("group") == g)
        ev_dom[g] = ev_scale(e)
        n = e["n_items"].unique().to_list()
        c = e["n_clusters"].unique().to_list()
        assert len(n) == 1 and len(c) == 1, (g, n, c)
        ev_n[g] = (f'n = {n[0]:,} pair-halves from {c[0]:,} '
                   f'{"Pfam families" if g == "Pfam seed" else "SCOPe folds"}.')
        ev_head.append(f'<th class="evh" data-group="{esc(g)}">{esc(words)}<span class="key">{n[0]:,} pair-halves, '
                       f'{c[0]:,} {"families" if g == "Pfam seed" else "folds"}</span>{ev_axis(ev_dom[g])}</th>')
    k_head = "".join(f'<th class="kph" data-bin="{esc(b)}">{IDENT_WORDS[b]}<span class="key">{n_kappa[b]:,}</span></th>' for b in IDENT)
    n_k = {r["identity_bin"]: r["n_pairs"] for r in t["usable_k"].filter(pl.col("alphabet") == "protein20").iter_rows(named=True)}
    chips = "".join(f'<button type="button" class="chip" data-kbin="{esc(b)}" aria-pressed="{"true" if i == 0 else "false"}">'
                    f'{b.replace("-", "–")} identity</button>' for i, b in enumerate(K_BINS))

    rows = []
    for gname, alphs in SIZE_GROUPS:
        rows.append(f'<tr class="sizegrp"><td colspan="{2 + len(EV_GROUPS) + len(IDENT) + 1}">{gname}</td></tr>')
        for a in alphs:
            n_cls = t["classes"].filter(pl.col("alphabet") == a)["n_classes"][0]
            evc = "".join(ev_cell(fair.filter((pl.col("group") == g) & (pl.col("alphabet") == a)).row(0, named=True),
                                  ev_dom[g], ev_n[g]) for g, _ in EV_GROUPS)
            rows.append(f'<tr data-alphabet="{esc(a)}"><th class="txt mono" scope="row">{esc(a)}'
                        f'<span class="key">{n_cls} classes</span></th>'
                        f'<td class="txt classes">{class_cell(t, a)}</td>{evc}{kappa_cells(t, a, kmin, kmax)}{k_cell(t, a)}</tr>')

    sd = fair.filter((pl.col("group") == "Pfam seed") & (pl.col("alphabet") == "sdm12")).row(0, named=True)
    above = fair.filter((pl.col("group") == "Pfam seed") & (pl.col("lo") > 1) & (pl.col("alphabet") != "protein20"))
    above_names = ", ".join(above.sort("rank")["alphabet"].to_list())
    pf = fair.filter(pl.col("group") == "Pfam seed")
    two = pf.filter(pl.col("classes") == 2).sort("rank")
    two_below = bool((two["hi"] < 1).all())
    best_two = two.row(0, named=True)
    count_words = {0: "no alphabet keeps", 1: "one alphabet keeps", 2: "two alphabets keep", 3: "three alphabets keep",
                   4: "four alphabets keep"}.get(above.height, f"{above.height} alphabets keep")

    # All groups and pair selections of notebook 291, for readers who want the rest.
    mode_words = {"low_half": "fair test: 20% identity or less, grouped on one half, measured on the other",
                  "all": "all identities, whole pairs",
                  "low_same": "20% identity or less, grouped and measured on the same positions (biased)"}
    full_rows = []
    for g in ev["group"].unique(maintain_order=True).to_list():
        for m in ("low_half", "all", "low_same"):
            e = ev.filter((pl.col("group") == g) & (pl.col("mode") == m))
            if not e.height:
                continue
            n, c = e["n_items"][0], e["n_clusters"][0]
            hom = "homologs" if e["homologs_by_truth"][0] else "not homologs by SCOPe"
            vals = {r["alphabet"]: r for r in e.iter_rows(named=True)}
            cells = "".join(
                f'<td class="num full" data-alphabet="{esc(a)}" data-group="{esc(g)}" data-mode="{m}" '
                f'data-ratio="{vals[a]["ratio_to_20aa"]!r}" data-n="{n}">{vals[a]["ratio_to_20aa"]:.2f}</td>'
                for a in t["order"])
            full_rows.append(f'<tr><th class="txt" scope="row">{esc(g)} ({hom})<span class="key">{esc(mode_words[m])}; '
                             f'{n:,} {"pair-halves" if m == "low_half" else "pairs"}, {c:,} clusters</span></th>{cells}</tr>')
    full_head = '<th class="txt">group and pairs used</th>' + "".join(f'<th class="mono">{esc(a)}</th>' for a in t["order"])

    lede = ("For anyone building a prefilter on a reduced amino acid alphabet: three measures for each of "
            "kmerseek's 19 alphabets. How much evidence of relatedness one aligned position keeps. How far above "
            "chance the classes of aligned residues agree. Which seed lengths k can reach homologs.")
    body = f"""
<h1>Reduced alphabets</h1>
<p class="lede">{lede}</p>
<p class="prose" id="what-for">This page measures alphabets with the alignment given. Each pair of
proteins was aligned already, by Pfam's curators or by structure; the page asks how much of that
alignment each alphabet's classes keep. It does not rank search tools, and it does not say which
alphabet finds the most homologs in a search, because a search has to find the alignment first. Use it
to choose which alphabets and seed lengths are worth testing in a search, then test them on the
<a href="index.html">leaderboard</a>.</p>
<p class="verdict" id="al-lead" data-above="{esc(above_names)}" data-best-two="{esc(best_two["alphabet"])}">In the fair
test on Pfam seed pairs at 20% identity or less, {count_words} more evidence per aligned position than the
20 amino acids, with the 95% interval above 1 ({esc(above_names)}). sdm12 keeps
<span class="mono">{f3(sd["ratio_to_20aa"])}</span> times as much (95% interval
<span class="mono">{f3(sd["lo"])}–{f3(sd["hi"])}</span>). {"Every two-class alphabet keeps less, with its interval below 1; the best of them" if two_below else "The best two-class alphabet"},
{esc(best_two["alphabet"])}, keeps <span class="mono">{best_two["ratio_to_20aa"]:.2f}</span> times as much.
n = {sd["n_items"]:,} pair-halves from {sd["n_clusters"]:,} Pfam families.</p>

<h2 id="table">One row per alphabet</h2>
<div class="legend" aria-label="Key to the table">
<span class="item"><span class="cls"><span class="sq">A</span><span class="sq">S</span><span class="sq">T</span></span> one class: residues the alphabet treats as the same letter</span>
<span class="item"><span class="cls hp-h"><span class="sq">h</span></span> hydrophobic</span>
<span class="item"><span class="cls hp-p"><span class="sq">p</span></span> polar</span>
<span class="item"><span class="cls hp-c"><span class="sq">c</span></span> cysteine in its own class (H/P alphabets only)</span>
<span class="item"><svg width="44" height="14" aria-hidden="true"><line x1="22" x2="22" y1="0" y2="14" stroke="var(--muted)" stroke-dasharray="3 2"/><line x1="8" x2="38" y1="7" y2="7" stroke="var(--ink)" stroke-width="1.6"/><circle cx="26" cy="7" r="3.4" fill="var(--ink)"/></svg> evidence ratio with its 95% interval; dashed line at 1, the 20 amino acids</span>
<span class="item"><span class="sw" style="background:color-mix(in srgb, var(--accent) 12%, var(--surface))"></span><span class="sw" style="background:color-mix(in srgb, var(--accent) 82%, var(--surface))"></span> kappa, low to high</span>
<span class="item"><svg width="60" height="14" aria-hidden="true"><line x1="2" x2="58" y1="7" y2="7" stroke="var(--muted)" stroke-width="1.2"/><line x1="58" x2="58" y1="3" y2="11" stroke="var(--muted)" stroke-width="1.2"/></svg> k values where a seed matches few proteins by chance (notebook 274)</span>
<span class="item"><span class="sw" style="background:var(--k-use)"></span> usable k, at least 20% of pairs reachable</span>
<span class="item"><span class="sw" style="background:var(--k-use-light)"></span> usable k, 5–20% of pairs reachable</span>
<span class="item"><span class="sw" style="background:var(--null)"></span> fewer than {MIN_N} pairs</span>
</div>
<details class="prose al-howto"><summary>How to read each column</summary>
<p><b>Classes.</b> The squares in one box are residues the alphabet puts in one class, so a seed
cannot tell them apart. In the seven hydrophobic/polar (H/P) alphabets, a box is amber for hydrophobic,
blue for polar, grey for a cysteine class of its own. The classes are kmerseek's, copied into the analysis
repository; the H/P letters are {PR[116]}'s.</p>
<p><b>Evidence per aligned position</b> ({NB[291]}, {PR[99]}). For one aligned pair, an aligned position
is a column where neither protein has a gap. Each residue is replaced by its class. Agreement <i>a</i> is
the share of positions where the two classes are the same; chance <i>c</i> is the share expected from
the two proteins' class make-up alone. The evidence in bits is
<i>a</i>·log₂(<i>a</i>/<i>c</i>) + (1−<i>a</i>)·log₂((1−<i>a</i>)/(1−<i>c</i>)): 0 when the pair agrees
no more than chance, higher when one aligned position says more about whether the proteins are
related. The table shows each alphabet's mean evidence divided by the 20 amino acids' on the same pairs:
right of the dashed line keeps more. The 95% interval comes from resampling Pfam families, or the query
domain's SCOPe fold, 2,000 times. <b>Fair test:</b> pairs are put in the 20% or less group by the
identity of half of their positions (odd or even) and measured on the other half, then the halves are
swapped. Picking pairs by the same positions that are measured favours every alphabet that merges
letters (the biased view is in the full table below).</p>
<p><b>Kappa</b> ({NB[230]}, {PR[100]}). Agreement corrected for chance: (<i>a</i> − <i>c</i>) / (1 − <i>c</i>),
the mean over Pfam seed pairs in each identity bin. 0 means no better than a shuffled partner, 1 means
every aligned position agrees. Identity is the 20 amino acid identity over aligned positions. A higher
kappa for a smaller alphabet does not mean it finds more homologs: merging letters also counts the
common swaps between homologs (I and V, D and E) as matches.</p>
<p><b>Usable k</b> ({NB[230]}, ranges from {NB[274]}, {PR[101]}). An exact seed of length k can find a
pair only if the pair shares an unbroken run of k aligned positions in the same class. A seed length is
usable when it is inside notebook 274's range: no more than 100 human proteins share a seed by chance,
and up to one chance match in Swiss-Prot. It must also reach pairs: at least 5% of the Pfam seed pairs
in that identity bin share such a run. The grey line is notebook 274's range; the bar is the usable part of it.</p>
</details>
<div class="kbins" role="group" aria-label="Identity bin for the usable k column">Usable k for Pfam seed pairs at: {chips}</div>
<div class="tablebox"><table class="al">
<thead><tr><th class="txt" rowspan="2">alphabet</th><th class="txt" rowspan="2">classes</th>
<th colspan="{len(EV_GROUPS)}">evidence per aligned position, times that of the 20 amino acids<span class="key">fair test, 20% identity or less</span></th>
<th colspan="{len(IDENT)}">kappa by identity, Pfam seed pairs<span class="key">n pairs under each bin</span></th>
<th rowspan="2" class="ukh">usable seed length k<span class="key" id="kbin-n">n = {n_k[K_BINS[0]]:,} pairs</span>{k_axis()}</th></tr>
<tr>{"".join(ev_head)}{k_head}</tr></thead>
<tbody>{"".join(rows)}</tbody></table></div>
<p class="muted small">Hover over a cell for its number and n. Each evidence column has its own scale. On a
narrow screen the table scrolls sideways; the alphabet names stay in view.</p>

<h2 id="hp">The seven H/P alphabets, residue by residue</h2>
<p class="prose">From {PR[116]}: one row per H/P alphabet, one square per amino acid, as kmerseek v0.4.0
classes it. The seven agree on 15 residues and differ only on C, G, P, W and Y.</p>
<div class="legend" aria-label="Key to the H/P table">
<span class="item"><span class="cls hp-h"><span class="sq">h</span></span> hydrophobic</span>
<span class="item"><span class="cls hp-p"><span class="sq">p</span></span> polar</span>
<span class="item"><span class="cls hp-c"><span class="sq">c</span></span> cysteine in its own class</span>
</div>
<div class="tablebox">{hp_figure(t)}</div>

<h2 id="full">Every group and pair selection of notebook 291</h2>
<details><summary>Evidence ratio to the 20 amino acids for all 19 alphabets, five groups of pairs, three ways of picking pairs</summary>
<p class="prose">The two SCOPe groups marked "not homologs by SCOPe" are pairs in different superfamilies.
There, more evidence of relatedness is a flaw, not a gain: two-class alphabets lead in those groups,
which may come from the structural superposition putting buried positions on buried positions.</p>
<div class="tablebox"><table class="lb full"><thead><tr>{full_head}</tr></thead><tbody>{"".join(full_rows)}</tbody></table></div>
</details>

<h2 id="sources">Where the numbers come from</h2>
<p class="prose">All tables are in <span class="mono">results/leaderboard/alphabets/</span> of
<a href="{ANALYSIS}/tree/{src["commit"]}/results/leaderboard/alphabets">seanome/2024-kmerseek-analysis at
commit <span class="mono">{src["commit"][:8]}</span></a>, written by <span class="mono">make leaderboard-alphabets</span>.
Evidence values are {NB[291]}'s stored table; the build adds the pair counts and checks that its own mean
evidence equals the notebook's. Kappa and usable k are computed with {NB[230]}'s code from its per-pair
table, because the notebook prints them but does not store them. Data: Pfam-A 38.2 seed alignments
(pairs with at least 50 aligned positions) and SCOPe 2.08 domains at 40% identity aligned by structure
with USalign (TM-score at least 0.5).</p>
"""
    data = {
        "note": "Values as stored in data/alphabets/*.parquet. evidence: ratio_to_20aa with its 95% interval "
                "(lo, hi), n_items pairs or pair-halves, n_clusters; kappa per identity_bin with n_pairs; "
                "usable_k per identity_bin with n_pairs, k_min and k_max from notebook 274.",
        "classes": {a: ["".join(c["letters"]) for c in classes_of(t, a)] for a in t["order"]},
        "hp_classes_pr116": [{k: r[k] for k in ("alphabet", "residue", "class")} for r in t["hp116"]],
        "evidence": ev.select("alphabet", "classes", "group", "mode", "pairs_used", "homologs_by_truth", "ratio_to_20aa",
                              "lo", "hi", "I_bits", "first_in", "rank", "n_items", "item_kind", "n_clusters",
                              "cluster_kind").to_dicts(),
        "kappa": t["kappa"].to_dicts(),
        "usable_k": t["usable_k"].to_dicts(),
    }
    return body, data
