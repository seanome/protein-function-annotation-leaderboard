"""Methods section for SCOPe40: how a hit is labelled, and why scoring waits.

Reads data/scope40_testbed.json (written by scripts/scope40_testbed.py). Every number in the
section comes from that file or from label_hits() on the made-up query of Figure 18, so the
prose and the figure cannot disagree.
"""

import html
import json
from pathlib import Path

FIG = 19  # Figures 11-17 are the leaderboard and pile-up figures; 18 is Known cases (PR 30)


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def split(sccs: str) -> tuple[str, str, str, str]:
    p = sccs.split(".")
    return p[0], ".".join(p[:2]), ".".join(p[:3]), sccs


def label(query: str, target: str, level: str, rule: str) -> str:
    """'tp', 'fp' or 'skip' for one hit at one level.

    rule 'paper' is Foldseek's: family TP = same family; superfamily TP = same superfamily,
    different family; fold TP = same fold, different superfamily. rule 'nested' counts every
    same-superfamily hit as a superfamily TP. Under both, a different fold is a false positive
    and anything else is not counted.
    """
    q, t = split(query), split(target)
    if q[1] != t[1]:
        return "fp"
    same_fam, same_sf = q[3] == t[3], q[2] == t[2]
    if level == "family":
        return "tp" if same_fam else "skip"
    if level == "superfamily":
        if rule == "nested":
            return "tp" if same_sf else "skip"
        return "tp" if same_sf and not same_fam else "skip"
    if level == "fold":
        return "tp" if not same_sf else "skip"
    raise ValueError(level)


def sensitivity(query: str, hits: list[str], level: str, rule: str, n_possible: int) -> dict:
    """True positives ranked above the first false positive, over the possible true positives."""
    labels = [label(query, h, level, rule) for h in hits]
    cut = labels.index("fp") if "fp" in labels else len(labels)
    above = sum(x == "tp" for x in labels[:cut])
    return {"labels": labels, "cut": cut, "above": above, "n_possible": n_possible,
            "sens": above / n_possible if n_possible else None}


# The made-up query of Figure 18: an immunoglobulin-like domain, b.1.1.1, and its eight hits
# in score order, best first, after its own self-hit was dropped. Its relatives in the
# made-up database: 4 others in its family, 5 in its superfamily outside the family, 3 in its
# fold outside the superfamily.
EX_QUERY = "b.1.1.1"
EX_HITS = ["b.1.1.1", "b.1.1.4", "b.1.2.1", "b.1.1.2", "b.1.1.1", "a.4.5.28", "b.1.1.1", "b.2.2.2"]
EX_N = {"family": 4, "superfamily": 5, "fold": 3}
EX_ROWS = [("Family (paper)", "family", "paper", EX_N["family"]),
           ("Superfamily (paper)", "superfamily", "paper", EX_N["superfamily"]),
           ("Fold (paper)", "fold", "paper", EX_N["fold"]),
           ("Superfamily (nested)", "superfamily", "nested", EX_N["family"] + EX_N["superfamily"])]


def worked_example() -> list[dict]:
    out = []
    for name, level, rule, n in EX_ROWS:
        r = sensitivity(EX_QUERY, EX_HITS, level, rule, n)
        # Foldseek's script divides the superfamily and fold counts by one more than this.
        r["n_script"] = n + 1 if (rule == "paper" and level != "family") else None
        out.append({"name": name, "level": level, "rule": rule, **r})
    return out


def figure_svg(rows: list[dict]) -> str:
    LW, BW, BH, GAP, TOP, RH = 150, 52, 26, 5, 58, 40
    x0 = LW
    xr = x0 + len(EX_HITS) * (BW + GAP) + 14
    W, H = xr + 205, TOP + len(rows) * RH + 6
    s = [f'<svg viewBox="0 0 {W} {H}" width="{W}" role="img" aria-label="Eight made-up hits labelled '
         f'under each rule, with the cut at the first false positive">']
    s.append(f'<text class="ex-lab" x="0" y="18">Hit, best score first</text>')
    s.append(f'<text class="ex-lab" x="0" y="40">SCOPe id of the hit</text>')
    for i, h in enumerate(EX_HITS):
        cx = x0 + i * (BW + GAP) + BW / 2
        s.append(f'<text class="ex-tk" x="{cx}" y="18" text-anchor="middle">{i + 1}</text>')
        s.append(f'<text class="ex-tk" x="{cx}" y="40" text-anchor="middle">{esc(h)}</text>')
    s.append(f'<text class="ex-lab" x="{xr}" y="40">sensitivity</text>')
    cut = rows[0]["cut"]
    for j, r in enumerate(rows):
        y = TOP + j * RH
        s.append(f'<text class="ex-lab" x="0" y="{y + BH / 2 + 5}">{esc(r["name"])}</text>')
        for i, lab in enumerate(r["labels"]):
            bx = x0 + i * (BW + GAP)
            dim = ' opacity="0.4"' if i > r["cut"] else ""
            if lab == "tp":
                s.append(f'<g data-label="tp"{dim}><rect class="ex-kept" x="{bx}" y="{y}" width="{BW}" height="{BH}" rx="3"/>'
                         f'<text class="ex-glyph" x="{bx + BW / 2}" y="{y + 18}" text-anchor="middle" '
                         f'style="fill:var(--bg)">✓</text></g>')
            elif lab == "fp":
                s.append(f'<g data-label="fp"{dim}><rect x="{bx + 1}" y="{y + 1}" width="{BW - 2}" height="{BH - 2}" rx="3" '
                         f'style="fill:var(--surface);stroke:var(--ink);stroke-width:2"/>'
                         f'<text class="ex-glyph" x="{bx + BW / 2}" y="{y + 18}" text-anchor="middle">✗</text></g>')
            else:
                s.append(f'<g data-label="skip"{dim}><rect class="ex-drop" x="{bx + 1}" y="{y + 1}" width="{BW - 2}" '
                         f'height="{BH - 2}" rx="3"/></g>')
        txt = f'{r["above"]} of {r["n_possible"]} = {r["sens"]:.2f}'
        s.append(f'<text class="ex-lab mono" x="{xr}" y="{y + 12}">{txt}</text>')
        if r["n_script"]:
            s.append(f'<text class="ex-note" x="{xr}" y="{y + 28}">Foldseek\'s original script: {r["above"]} of '
                     f'{r["n_script"]} = {r["above"] / r["n_script"]:.2f}</text>')
    cx = x0 + cut * (BW + GAP) - GAP / 2
    s.append(f'<line class="ex-cut" x1="{cx}" x2="{cx}" y1="{TOP - 8}" y2="{H - 2}"/>')
    s.append("</svg>")
    return "".join(s)


def legend_html() -> str:
    sw = lambda inner: f'<svg width="30" height="18" viewBox="0 0 30 18" aria-hidden="true">{inner}</svg>'
    tp = sw('<rect class="ex-kept" x="1" y="1" width="28" height="16" rx="3"/>'
            '<text class="ex-glyph" x="15" y="14" text-anchor="middle" style="fill:var(--bg);font-size:12px">✓</text>')
    fp = sw('<rect x="2" y="2" width="26" height="14" rx="3" style="fill:var(--surface);stroke:var(--ink);stroke-width:2"/>'
            '<text class="ex-glyph" x="15" y="14" text-anchor="middle" style="font-size:12px">✗</text>')
    sk = sw('<rect class="ex-drop" x="2" y="2" width="26" height="14" rx="3"/>')
    cut = sw('<line class="ex-cut" x1="15" x2="15" y1="0" y2="18"/>')
    pale = (sw('<g opacity="0.4"><rect class="ex-kept" x="1" y="1" width="28" height="16" rx="3"/>'
               '<text class="ex-glyph" x="15" y="14" text-anchor="middle" style="fill:var(--bg);font-size:12px">✓</text></g>')
            + sw('<g opacity="0.4"><rect x="2" y="2" width="26" height="14" rx="3" style="fill:var(--surface);stroke:var(--ink);stroke-width:2"/>'
                 '<text class="ex-glyph" x="15" y="14" text-anchor="middle" style="font-size:12px">✗</text></g>')
            + sw('<g opacity="0.4"><rect class="ex-drop" x="2" y="2" width="26" height="14" rx="3"/></g>'))
    return (f'<p class="ex-read" style="display:flex;flex-wrap:wrap;gap:4px 16px;align-items:center">'
            f'<span>{tp} true positive</span><span>{fp} false positive (a different fold)</span>'
            f'<span>{sk} not counted</span><span>{cut} first false positive: the list is cut here</span>'
            f'<span>{pale} any mark drawn pale (TP, FP or not counted) is below the cut, so it adds nothing</span></p>')


def section(data_dir: Path) -> str:
    p = Path(data_dir) / "scope40_testbed.json"
    if not p.exists():
        return ('<h2 id="scope40">SCOPe40: ranking a domain\'s relatives</h2>\n'
                f'<p>Not in this build: <span class="mono">data/{p.name}</span> is missing.</p>')
    tb = json.loads(p.read_text())
    n8, n1 = tb["scope208_astral40"], tb["scope201_foldseek"]
    q8, q1 = n8["n_queries_with_a_true_positive"], n1["n_queries_with_a_true_positive"]
    tq = tb["tea_query_set"]
    fs = next(x for x in tb["rocx"] if x["file"] == "foldseek.rocx")["levels"]
    mm = next(x for x in tb["rocx"] if x["file"] == "mmseqs.rocx")["levels"]
    ex = tb["example"]
    n_fs = next(x for x in tb["rocx"] if x["file"] == "foldseek.rocx")["n_queries"]
    rows = worked_example()
    fam, sf, fold, nest = rows
    sha = tb["inputs"]["astral208"]["sha256"]
    cv = tb["rocx_counts_vs_208"]
    nb76 = next(e for e in json.loads((Path(data_dir) / "experiments.json").read_text()) if e["id"] == "076")
    both = next(r for r in nb76["result"] if r["name"] == "benchmark queries both detect")
    only = next(r for r in nb76["result"] if r["name"] == "benchmark queries only kmerseek detects")
    n_det = int(both["value"]) + int(only["value"])
    n_mean = next(int(r["name"].split(", ")[1].split()[0]) for r in nb76["result"]
                  if r["name"].startswith("kmerseek mean per-query family sensitivity"))
    nb76_url = (f"https://github.com/seanome/2024-kmerseek-analysis/blob/{nb76['source']['commit']}/"
                f"{nb76['source']['path']}")

    def f2(x):
        return f"{x:.2f}"

    def f3(x):
        return f"{x:.3f}"

    return f"""<h2 id="scope40">SCOPe40: ranking a domain's relatives</h2>
<p><b>Scored so far: Foldseek, MMseqs2 and kmerseek</b>, on the
<a href="scope40.html">SCOPe40 page</a>. This section gives the rules. The question: given one known domain, does the tool rank the domain's
known relatives above the first unrelated hit? Here a relative is a domain of the same
family, superfamily or fold, one level at a time; an unrelated hit is a domain of a different
fold. The Foldseek (a structure search tool) and MMseqs2 (a sequence search tool) numbers are read
from per-query files that other authors published; kmerseek (the k-mer search tool whose settings
fill most rows of this site) is scored from its own hit lists.</p>
<p>SCOPe (Structural Classification of Proteins, extended) sorts protein domains into four nested
levels. An id such as <span class="mono">b.1.1.1</span> reads class.fold.superfamily.family:
class b (all beta), fold 1, superfamily 1, family 1. Domains in one superfamily are taken to share
an ancestor. Domains in one fold share a shape, with or without an ancestor. A family numbered 0
holds SCOPe's automated matches: domains placed in the superfamily but in no curated family; two
domains in the same family 0 count as the same family, as in the scripts below. Each
domain has its own id, such as <span class="mono">d6jc0b_</span>: d, the PDB entry
(<span class="mono">6jc0</span>), the chain (<span class="mono">b</span>) and a domain number
(<span class="mono">_</span> for the whole chain). SCOPe40 is the set of SCOPe domains with no two
above 40% sequence identity, as distributed by ASTRAL, SCOPe's companion set of domain sequences
and structures (Chandonia et al. 2022).</p>

<h3 id="scope40-paper">What Foldseek's paper counts</h3>
<p>The paper that set this benchmark, van Kempen et al. (2024), defines the labels in its
Methods (quoted under its CC BY 4.0 licence):</p>
<blockquote style="margin:8px 0 8px 1em;padding-left:10px;border-left:3px solid var(--hair)">"Per domain, we
measured the fraction of detected TPs up to the first FP. For family-level, superfamily-level and
fold-level recognition, TPs were defined as same family, same superfamily and not same family and
same fold and not same superfamily, respectively. Hits from different folds are FPs."</blockquote>
<blockquote style="margin:8px 0 8px 1em;padding-left:10px;border-left:3px solid var(--hair)">"After sorting
the alignment result of each query (...), we calculated the sensitivity as the fraction of TPs in
the sorted list up to the first FP, all excluding self-hits. For comparison, we took the mean
sensitivity over all queries for family-level, superfamily-level and fold-level classifications.
We evaluated only SCOPe members with at least one other family, superfamily and fold
member."</blockquote>
<p>Its main text gives the denominator: "the fraction of TPs out of all possible correct matches
until the first FP, a match to a different fold". Its Figure 2a draws the per-query values as a
cumulative curve and says "Sensitivity is the area under the ROC (AUROC) curve up to the first FP".
The paper means each query's own value: the share of that query's TPs ranked above its first FP,
the sensitivity defined below. The AUC on this site is a second area, under those values sorted
across queries (defined below). The paper's numbers come
from a scoring script in Foldseek's analysis repository, <span class="mono">bench.noselfhit.awk</span>.</p>
<p>LoL-align (Reifenrath et al. 2026), a structure aligner from the same two groups (Johannes
Söding's and Martin Steinegger's), reuses this benchmark. Its paper
defines the score as "the fraction of true positives (TPs) found before the first FP in the list of
ranked matches", with "TPs ... from the same family, superfamily or fold and FPs ... from different
folds", over SCOPe entries "with at least one additional family-, superfamily- or fold-level
member": 3,566 queries on SCOPe 2.01. The sentence says "or", but 3,566 is the count of domains
with members at all three levels (the domain-set table below), which is what the script keeps
(line 24). Its copy of the
script (<span class="mono">soedinglab/LoLalign-analysis</span>,
<span class="mono">scope_benchmark/bench.noselfhit.awk</span>) is Foldseek's line for line, except
for the superfamily and fold denominators, lines 26 and 27.</p>
<p>For one query, a <b>true positive</b> (TP) is a hit at the level being scored. A <b>false
positive</b> (FP) is a hit in a different fold. Any other hit is <b>not counted</b>: neither TP
nor FP. Sort the hits by the tool's own score, best first, after dropping the query's hit to
itself. Walk down the list and stop at the first FP. The <b>sensitivity up to the first false
positive</b> is the number of TPs above that cut divided by the number of possible TPs for that
query in the whole set. The walk uses every hit the tool reports at the settings it was run with.
Every query on this site has at least one possible TP at each level, so no sensitivity is 0 of 0.</p>
<p>This site first planned a <b>nested rule</b>. A hit is a superfamily TP when query and target
share the first three fields of the SCOPe id (class letter, fold, superfamily), family included.
It is an FP when their folds differ, and not counted when they share the fold only. <a href="scope40.html#fig-scope40-labels">Figure {FIG}</a>, on the
SCOPe40 page, labels one made-up hit list both ways.</p>

<p>Over many queries, the <b>area under the curve</b> (AUC) sorts the queries by sensitivity,
highest first, and plots sensitivity against the fraction of queries. Its area is the mean
sensitivity over all queries, so the AUC and the mean are the same number. <b>Coverage</b> is the
share of queries with at least one hit at the settings the tool was run with; it goes next to every
AUC. Coverage says nothing about whether the hits are right: the shuffled-query control has
coverage near 100% too.</p>
<p><b>Every query counts.</b> A query that gets no hits scores 0 and stays in every mean and every
area. An earlier kmerseek claim, that it "outperforms Foldseek on SCOPe", was computed over only
part of the queries. That claim and its number are retired. Its source, <a href="{nb76_url}">notebook 076</a>, ran one hydrophobic-polar alphabet
at k = 26 with a strict cutoff ({esc(nb76["dataset"])}), so kmerseek reported hits for few queries:
it found a family member for {n_det} of the {n_fs:,} queries of Foldseek's file
({100 * n_det / n_fs:.1f}%), and its mean was taken over the {n_mean} queries that both tools
returned any hit for (a hit need not be a family member, so {n_mean} is more than {n_det}). The
kmerseek settings on the SCOPe40 page are different ones and report hits for nearly every query.</p>

<h3 id="scope40-differences">The paper's rule against the planned rule</h3>
<p>The two rules differ in the rows below. On 2026-10-17 Olga Botvinnik, who runs this
benchmark, chose the paper's rule as its Methods and its script state it ("Foldseek paper" and
"Foldseek's original script" columns), with LoL-align's denominators (below). The planned rule is
not used.</p>
<div class="tablebox"><table class="lb"><thead><tr><th class="txt" scope="col">Point</th>
<th class="txt" scope="col">Foldseek paper (Methods)</th><th class="txt" scope="col">Foldseek's original script</th>
<th class="txt" scope="col">This site's planned rule</th></tr></thead><tbody>
<tr><td class="txt">Superfamily TP</td><td class="txt">same superfamily, different family</td>
<td class="txt">same</td><td class="txt">same superfamily, family included</td></tr>
<tr><td class="txt">Same-family hit, superfamily level</td><td class="txt">not counted</td><td class="txt">not counted</td>
<td class="txt">TP</td></tr>
<tr><td class="txt">Fold TP</td><td class="txt">same fold, different superfamily</td><td class="txt">same</td>
<td class="txt">no fold level given</td></tr>
<tr><td class="txt">Superfamily denominator</td><td class="txt">all possible TPs: superfamily members outside the family</td>
<td class="txt">one more than that (LoL-align's copy: not)</td><td class="txt">every other superfamily member</td></tr>
<tr><td class="txt">Fold denominator</td><td class="txt">all possible TPs: fold members outside the superfamily</td>
<td class="txt">one more than that (LoL-align's copy: not)</td><td class="txt">no fold level given</td></tr>
<tr><td class="txt">Queries scored</td><td class="txt">domains with another member at all three levels</td>
<td class="txt">same</td><td class="txt">every query with a possible TP at the level</td></tr>
<tr><td class="txt">Hit to a domain missing from the script's domain list (it happens when the
list used for scoring is shorter than the set searched)</td><td class="txt">not stated</td><td class="txt">FP</td>
<td class="txt">not stated</td></tr>
</tbody></table></div>
<p>The paper's Figure 2a legend reads "TPs are matches within the same superfamily", which is the
nested rule; its Methods, quoted above, and its script leave the family out. A third rule, in the
Methods subsection that fits Foldseek's probability of a true match, ignores "hits to the same
family or to another superfamily". That rule only sets Foldseek's reported probability; it
changes no score here.</p>
<p>Where the extra 1 comes from: Foldseek's original script counts the possible family TPs as the family size minus
1, which removes the query itself. For the superfamily it takes the superfamily size minus (the
family size minus 1), and for the fold the fold size minus (the superfamily size minus 1). Those
two leave the query in, so each is one more than the number of domains that could be found.
In the script: lines 7, 10 and 13 count every domain, the query included, into its family,
superfamily and fold; line 17 drops the self-hit from the hit list; line 25 takes 1 off the family
count, but lines 26 and 27 subtract only the other family members (or the other superfamily
members), so the query stays in. It can never be found at those levels: its hit to itself is
dropped, and would count at the family level anyway. So no query can reach a superfamily or fold
sensitivity of 1 there. LoL-align's copy subtracts the whole family size and the whole superfamily
size instead, so its denominators are the domains that can be found. This site uses LoL-align's
denominators (chosen by Olga Botvinnik, 2026-10-17).</p>
<p>The extra 1 is in the per-query files of Foldseek and MMseqs2 on SCOPe 2.08 that other authors
published, as copied into seanome/2024-kmerseek-analysis
(<span class="mono">data/tea_scope40_rocx_files</span>). Query <span class="mono">{esc(ex['name'])}</span>
(<span class="mono">{esc(ex['scop'])}</span>) has {ex['SFAMCNT']} domains in its superfamily and
{ex['FAMCNT']} in its family, itself included in both, so {ex['sfam_possible']} possible
superfamily TP. Foldseek found it and the file says {ex['SFAM']:g}: {ex['sfam_tp']} of
{ex['sfam_denominator_script']}. Over the {n_fs:,} queries of the Foldseek file, the mean
superfamily sensitivity is {f3(fs['SFAM']['mean_as_published'])} as published and
{f3(fs['SFAM']['mean_over_possible_true_positives'])} over the possible TPs; at the fold level
{f3(fs['FOLD']['mean_as_published'])} and {f3(fs['FOLD']['mean_over_possible_true_positives'])}.
For MMseqs2: superfamily {f3(mm['SFAM']['mean_as_published'])} and
{f3(mm['SFAM']['mean_over_possible_true_positives'])}, fold {f3(mm['FOLD']['mean_as_published'])}
and {f3(mm['FOLD']['mean_over_possible_true_positives'])}. The family level has no extra 1
({f3(fs['FAM']['mean_as_published'])} for Foldseek either way). This site rescales those files to
LoL-align's denominators: each row's value times its original denominator is a whole number of TPs,
which is divided by the domains that row could find. These means are over the file's {n_fs:,} rows;
the SCOPe40 page's AUCs are over all {q8['all_three']:,} queries, with the
{q8['all_three'] - n_fs} the files lack scored 0, so they are a little lower. The Foldseek paper's own numbers are on the
original scale. The kmerseek SCOPe40 scorer of PR 107 already divides by the domains that can be
found, so its numbers that divide by every true pair in the set (column
<span class="mono">auc_ref_db</span>) are on the LoL-align scale too.</p>
<p>At the fold level, a TP is not always a relative: two domains can share a fold with no shared
ancestor (Cheng et al. 2008). So a fold-level TP means "same shape", not "homolog".</p>

<h3 id="scope40-set">The domain set</h3>
<p>Foldseek's set is SCOPe 2.01 at 40% identity, {n1['n_domains']:,} domains (Foldseek's list
<span class="mono">scop_lookup.fix.tsv</span>, as copied into seanome/2024-kmerseek-analysis). The
set every kmerseek SCOPe40 run in seanome/2024-kmerseek-analysis searched is a different release:
ASTRAL SCOPe 2.08 at 40% identity, {n8['n_domains']:,} domains
(<span class="mono">astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa</span>, sha256
<span class="mono">{sha[:12]}…</span>). Numbers from the two are not comparable without saying so.</p>
<div class="tablebox"><table class="lb"><thead><tr><th class="txt" scope="col">Domains with a possible TP</th>
<th scope="col">SCOPe 2.08, ASTRAL 40%</th><th scope="col">SCOPe 2.01, Foldseek</th></tr></thead><tbody>
<tr><td class="txt">All domains</td><td class="num">{n8['n_domains']:,}</td><td class="num">{n1['n_domains']:,}</td></tr>
<tr><td class="txt">Family level</td><td class="num">{q8['family']:,}</td><td class="num">{q1['family']:,}</td></tr>
<tr><td class="txt">Superfamily level, paper</td><td class="num">{q8['superfamily']:,}</td><td class="num">{q1['superfamily']:,}</td></tr>
<tr><td class="txt">Superfamily level, nested</td><td class="num">{q8['superfamily_nested']:,}</td><td class="num">{q1['superfamily_nested']:,}</td></tr>
<tr><td class="txt">Fold level, paper</td><td class="num">{q8['fold']:,}</td><td class="num">{q1['fold']:,}</td></tr>
<tr><td class="txt">All three levels (the paper's queries)</td><td class="num">{q8['all_three']:,}</td><td class="num">{q1['all_three']:,}</td></tr>
</tbody></table></div>
<p>The published Foldseek and MMseqs2 per-query files on SCOPe 2.08 hold {tq['n_foldseek_rocx']:,}
of the {q8['all_three']:,} all-three-level domains; {len(tq['missing_from_foldseek_rocx'])} are
missing. Their family, superfamily and fold sizes do not match SCOPe 2.08 on {cv['any']:,} of
the {cv['n_rows']:,} rows (fold size on {cv['FOLDCNT']:,}, superfamily on {cv['SFAMCNT']:,}, family on
{cv['FAMCNT']:,}){", always smaller" if cv['smaller_only'] else ""}. So the runs behind those files
searched a domain list with some 2.08 domains missing, and a tool scored on the full 2.08 set is not
on the same footing as them.</p>

<h3 id="scope40-choices">Choices made (Olga Botvinnik, 2026-10-17)</h3>
<ol>
<li>The rule: the paper's Methods and script. Its Figure 2a legend reads differently (above).</li>
<li>The denominator: LoL-align's, the domains that can be found at each level. Foldseek's
original, one more at the superfamily and fold levels, is kept beside it in the long table
(<span class="mono">sensitivity_foldseek_script</span>), and the published Foldseek and MMseqs2
files are rescaled from it.</li>
<li>The queries: Foldseek's set, domains with another member at all three levels, computed on
SCOPe 2.08: {q8['all_three']:,}.</li>
<li>The set: SCOPe 2.08 at 40% identity, {n8['n_domains']:,} domains, searched all against all.</li>
<li>Hit lists: kmerseek is scored from its own hit lists (one row per protein pair, kept by the
runs of <a href="https://github.com/seanome/2024-kmerseek-analysis/pull/107">analysis PR 107</a>).
Foldseek and MMseqs2 come from per-query files with no hit lists. The pipeline that keeps hit lists
for every other tool is written, not launched
(<a href="https://github.com/seanome/2024-kmerseek-analysis/pull/166">analysis PR 166</a>).</li>
</ol>
<h3 id="scope40-scores">The score each tool is ranked by</h3>
<p>Each tool's hits are sorted by its own score, best first, ties broken by target id. In the
pipeline: Foldseek, ProstT5 + Foldseek, MMseqs2 (both modes), phmmer, jackhmmer and HH-suite by
bit score, higher first; kmerseek by region E-value, smaller first. The two published files do not
record the order their runs used. kmerseek's region E-value is the E-value of the best-matching
region of the pair: about how many regions that good one search would find by chance.</p>
<p class="refs">van Kempen, M., S. S. Kim, C. Tumescheit, M. Mirdita, J. Lee, C. L. M. Gilchrist,
J. Söding, and M. Steinegger. "Fast and Accurate Protein Structure Search with Foldseek."
<i>Nature Biotechnology</i> 42, no. 2 (2024): 243–46.
<a href="https://doi.org/10.1038/s41587-023-01773-0">https://doi.org/10.1038/s41587-023-01773-0</a>.<br>
Chandonia, J.-M., L. Guan, S. Lin, C. Yu, N. K. Fox, and S. E. Brenner. "SCOPe: Improvements to the
Structural Classification of Proteins – Extended Database to Facilitate Variant Interpretation and
Machine Learning." <i>Nucleic Acids Research</i> 50, no. D1 (2022): D553–59.
<a href="https://doi.org/10.1093/nar/gkab1054">https://doi.org/10.1093/nar/gkab1054</a>.<br>
Reifenrath, L., M. van Kempen, G. Kim, S. H. Kim, M. Radnezhad, C. L. M. Gilchrist, M. Mirdita,
M. Steinegger, and J. Söding. "LoL-align: Sensitive and Fast Probabilistic Protein Structure
Alignment." <i>Nature Methods</i> (2026).
<a href="https://doi.org/10.1038/s41592-026-03249-w">https://doi.org/10.1038/s41592-026-03249-w</a>.<br>
Cheng, H., B.-H. Kim, and N. V. Grishin. "MALISAM: A Database of Structurally Analogous Motifs in
Proteins." <i>Nucleic Acids Research</i> 36, Database issue (2008): D211–17.
<a href="https://doi.org/10.1093/nar/gkm698">https://doi.org/10.1093/nar/gkm698</a>.</p>"""
