"""Draw functions for Figures 15 to 17, read from the region pile-up tables (data/pileup/):
the true Pfam regions of the human proteins, and every call a tool made of the same family
over each region's window, with is_tp ("correct") as the benchmark set it.

As in figures.py, each function takes one table from a data step in build_site.py and
returns the figure's HTML and a sentence built from the same rows. `sets` lists the tools
drawn, in order: dicts with key, name, label, kind.
"""

from __future__ import annotations

import math

import polars as pl

from figures import esc


def set_name(s: dict) -> str:
    return s["name"] + (f" ({s['label']})" if s["key"].startswith("kmerseek") else "")


# --------------------------------------------------------------------------- Figure 15

def draw_upset(t: pl.DataFrame, sets: list[dict], max_cols: int = 24) -> tuple[str, str]:
    """Figure 15: which tools found each true region, as an UpSet plot.

    Table columns: row_id, found_by (list of set keys whose tool made a correct call on the
    region). One column per combination of tools, largest first. Every combination of one
    tool alone is kept even when small, and so is "found by none"; combinations of two or
    more fill the remaining columns up to max_cols. A click on a column lists its regions
    (assets/figures.js reads data/upset.json).
    """
    keys = [s["key"] for s in sets]
    combos = (t.with_columns(pl.col("found_by").list.sort().list.join("+").alias("combo"))
              .group_by("combo").agg(pl.len().alias("n")).sort(["n", "combo"], descending=[True, False]))
    rows = combos.to_dicts()
    keep = [r for r in rows if r["combo"] and "+" not in r["combo"]]
    rest = [r for r in rows if "+" in r["combo"]]
    none = [r for r in rows if not r["combo"]]
    room = max(0, max_cols - len(keep))
    cols = sorted(keep + rest[:room], key=lambda r: (-r["n"], r["combo"])) + none
    left_out = rest[room:]
    size = {k: t.filter(pl.col("found_by").list.contains(k)).height for k in keys}
    by_key = {s["key"]: s for s in sets}
    LW, SW, CW, top, rowh = 230, 70, 24, 150, 22
    W = LW + SW + CW * len(cols) + 10
    H = top + rowh * len(keys) + 14
    nmax = max(r["n"] for r in cols)
    smax = max(size.values())
    by = lambda n: top - 18 - n / nmax * (top - 40)  # noqa: E731
    p = []
    for j in range(0, len(sets), 2):
        p.append(f'<rect x="0" y="{top + rowh * j:.1f}" width="{W}" height="{rowh}" fill="var(--chip)"/>')
    p.append(f'<text x="{LW + SW - 4}" y="{top - 22}" text-anchor="end" class="fg-small">regions in</text>')
    p.append(f'<text x="{LW + SW - 4}" y="{top - 8}" text-anchor="end" class="fg-small">the column</text>')
    p.append(f'<text x="{LW + SW - 4}" y="{top + rowh * len(keys) + 12}" text-anchor="end" class="fg-small">'
             'regions each tool found</text>')
    p.append(f'<line x1="{LW + SW}" x2="{W - 10}" y1="{top - 18}" y2="{top - 18}" class="fg-axis"/>')
    for i, r in enumerate(cols):
        x = LW + SW + CW * i + CW / 2
        members = r["combo"].split("+") if r["combo"] else []
        alone = len(members) == 1
        fill = f"var(--k-{by_key[members[0]]['kind']})" if alone else "var(--muted)"
        if not members:
            words = "found by none of these tools"
        elif alone:
            words = "found only by " + set_name(by_key[members[0]])
        else:
            words = "found by " + ", ".join(set_name(s) for s in sets if s["key"] in members) + " and no other tool"
        tip = f"{r['n']:,} regions {words}. Click to list them."
        style = f'fill="{fill}"' if members else 'fill="none" stroke="var(--muted)" stroke-dasharray="3 2"'
        p.append(f'<g class="fg-col" data-combo="{esc(r["combo"] or "none")}" data-tip="{esc(tip)}" tabindex="0" role="button">'
                 f'<rect x="{x - CW / 2:.1f}" y="20" width="{CW}" height="{H - 34}" fill="transparent"/>'
                 f'<rect x="{x - 8:.1f}" y="{by(r["n"]):.1f}" width="16" height="{top - 18 - by(r["n"]):.1f}" {style}/>'
                 f'<text x="{x:.1f}" y="{by(r["n"]) - 4:.1f}" text-anchor="middle" class="fg-tick">{r["n"]:,}</text>')
        ys = [top + rowh * keys.index(m) + rowh / 2 for m in members]
        if len(ys) > 1:
            p.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{min(ys):.1f}" y2="{max(ys):.1f}" stroke="var(--ink)" stroke-width="2"/>')
        for j, k in enumerate(keys):
            y = top + rowh * j + rowh / 2
            on = k in members
            p.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{5 if on else 3.5}" fill="var(--{"ink" if on else "hair"})"/>')
        p.append('</g>')
    for j, s in enumerate(sets):
        y = top + rowh * j + rowh / 2
        p.append(f'<text x="{LW - 40}" y="{y + 4:.1f}" text-anchor="end" class="fg-small">{esc(set_name(s))}</text>')
        w = size[s["key"]] / smax * (SW - 8)
        p.append(f'<rect x="{LW + SW - 4 - w:.1f}" y="{y - 6:.1f}" width="{w:.1f}" height="12" fill="var(--k-{s["kind"]})" '
                 f'data-tip="{esc(set_name(s))} found {size[s["key"]]:,} of the {t.height:,} regions" tabindex="0"/>')
        p.append(f'<text x="{LW + SW - 8 - w:.1f}" y="{y + 4:.1f}" text-anchor="end" class="fg-tick">{size[s["key"]]:,}</text>')
    svg = (f'<div class="figscroll"><svg class="fg fg-upset" viewBox="0 0 {W} {H}" style="min-width:{min(W, 720)}px;'
           f'max-width:{W * 1.2:.0f}px" role="img" aria-label="UpSet plot: true regions found by each combination of tools">'
           f'{"".join(p)}</svg></div>')
    leg = ('<div class="legend">'
           '<span class="item"><svg width="16" height="14" aria-hidden="true"><circle cx="8" cy="7" r="5" fill="var(--ink)"/></svg>'
           'the tool found the region</span>'
           '<span class="item"><svg width="16" height="14" aria-hidden="true"><circle cx="8" cy="7" r="3.5" fill="var(--hair)"/></svg>'
           'it did not</span>'
           '<span class="item"><svg width="14" height="14" aria-hidden="true"><rect x="2" y="1" width="10" height="12" fill="var(--k-letters)"/></svg>'
           'bar in the colour of a kind of tool: found by one tool alone (its dot says which)</span>'
           '<span class="item"><svg width="14" height="14" aria-hidden="true"><rect x="2" y="1" width="10" height="12" fill="var(--muted)"/></svg>'
           'found by two or more</span>'
           '<span class="item"><svg width="14" height="14" aria-hidden="true"><rect x="2" y="1" width="10" height="12" fill="none" '
           'stroke="var(--muted)" stroke-dasharray="3 2"/></svg>found by none</span></div>')
    found = t.filter(pl.col("found_by").list.len() > 0).height
    alone = sorted(keep, key=lambda r: (-r["n"], r["combo"]))
    kms = [s for s in sets if s["key"].startswith("kmerseek")]
    km_alone = next((r["n"] for r in keep if kms and r["combo"] == kms[0]["key"]), 0)
    s = (f"{found:,} of the {t.height:,} true regions were found by at least one of these tools, and "
         f"{t.height - found:,} by none.")
    if alone:
        s += f" The most found by one tool alone: {set_name(by_key[alone[0]['combo']])}, {alone[0]['n']:,} regions."
    if kms:
        s += f" {set_name(kms[0])} alone found {km_alone:,}."
    if left_out:
        s += (f" {len(left_out)} smaller combinations of two or more tools ({sum(r['n'] for r in left_out):,} regions) "
              "are not drawn.")
    panel = ('<div class="fg-list" id="upset-list" aria-live="polite"><p class="muted">Click or tap a column to list '
             'its regions here, with their residues.</p></div>')
    return leg + svg + panel, s


def upset_lists(t: pl.DataFrame) -> dict:
    """data/upset.json: combination -> its regions [symbol, accession, Pfam, start, end, residues]."""
    out: dict = {}
    for r in t.sort("row_id").iter_rows(named=True):
        key = "+".join(sorted(r["found_by"])) or "none"
        out.setdefault(key, []).append([r["symbol"] or "", r["accession"], r["pfam_id"], r["true_start"],
                                        r["true_end"], r["residues"]])
    return out


# --------------------------------------------------------------------------- Figure 16

def protein_line_svg(r0: dict, calls: list[dict], W: int = 640) -> str:
    """One example: the protein window as a thin line, the true region as a box on it, and
    one row per tool under it with that tool's best call as a bar (filled: correct;
    outline: not correct; a dot and words: no call of this family)."""
    lo = r0["window_offset"]
    hi = lo + len(r0["window_seq"]) - 1
    L, R, rowh = 190, 14, 15
    px = lambda v: L + (v - lo) / max(hi - lo, 1) * (W - L - R)  # noqa: E731
    H = 44 + rowh * len(calls) + 22
    p = [f'<line x1="{px(lo):.1f}" x2="{px(hi):.1f}" y1="20" y2="20" stroke="var(--ink)" stroke-width="1.2"/>',
         f'<rect x="{px(r0["true_start"]):.1f}" y="13" width="{max(px(r0["true_end"]) - px(r0["true_start"]), 1.5):.1f}" '
         f'height="14" fill="var(--surface)" stroke="var(--ink)" stroke-width="1.5" data-tip="true {esc(r0["pfam_id"])} '
         f'region, residues {r0["true_start"]:,}-{r0["true_end"]:,} ({r0["true_end"] - r0["true_start"] + 1:,} aa)" tabindex="0"/>',
         f'<text x="{L - 6}" y="24" text-anchor="end" class="fg-small">{esc(r0["symbol"] or r0["accession"])}, '
         f'{esc(r0["pfam_id"])}</text>']
    for i, c in enumerate(calls):
        y = 42 + i * rowh
        p.append(f'<text x="{L - 6}" y="{y + 4}" text-anchor="end" class="fg-small">{esc(c["set_name"])}</text>')
        col = f"var(--k-{c['kind']})"
        if c["call_start"] is None:
            p.append(f'<text x="{L + 2}" y="{y + 4}" class="fg-small" fill="var(--muted)">· no call of this family</text>')
            continue
        a, b = max(c["call_start"], lo), min(c["call_end"], hi)
        tip = (f"{c['set_name']}: residues {c['call_start']:,}-{c['call_end']:,}, overlap with the true region "
               f"(IoU) {c['iou']:.2f}, score {c['score']:,.1f}, "
               + ("finds this region" if c["is_tp"] else "does not find this region"))
        if c["is_tp"]:
            p.append(f'<rect x="{px(a):.1f}" y="{y - 4}" width="{max(px(b) - px(a), 1.5):.1f}" height="8" fill="{col}" '
                     f'data-tip="{esc(tip)}" tabindex="0"/>')
        else:
            p.append(f'<rect x="{px(a):.1f}" y="{y - 3.5}" width="{max(px(b) - px(a), 1.5):.1f}" height="7" '
                     f'fill="var(--surface)" stroke="{col}" stroke-width="1.5" data-tip="{esc(tip)}" tabindex="0"/>')
    ya = 42 + rowh * len(calls)
    p.append(f'<line x1="{px(lo):.1f}" x2="{px(hi):.1f}" y1="{ya}" y2="{ya}" class="fg-axis"/>')
    p.append(f'<text x="{px(lo):.1f}" y="{ya + 14}" class="fg-tick">{lo:,}</text>')
    p.append(f'<text x="{px(hi):.1f}" y="{ya + 14}" text-anchor="end" class="fg-tick">{hi:,}</text>')
    p.append(f'<text x="{(px(lo) + px(hi)) / 2:.1f}" y="{ya + 14}" text-anchor="middle" class="fg-small">residue</text>')
    return (f'<div class="figscroll"><svg class="fg" viewBox="0 0 {W} {H}" style="min-width:520px;max-width:{W * 1.25:.0f}px" '
            f'role="img" aria-label="calls on {esc(r0["symbol"] or r0["accession"])}">{"".join(p)}</svg></div>')


CASE_WORDS = {"best": "Best case: the highest overlap of a correct call",
              "median": "Middle case: the region in the middle when the found regions are ordered by overlap",
              "missed": "Missed: no tool made a correct call"}


def draw_gallery(t: pl.DataFrame) -> tuple[str, str]:
    """Figure 16: three real regions per identity bin.

    Table columns: stratum, stratum_label, stratum_order, n_regions, n_found, median_iou, case, row_id,
    symbol, accession, pfam_id, true_start, true_end, window_seq, window_offset, best_iou,
    set_name, kind, call_start, call_end, iou, is_tp, score: one row per example and tool,
    with the call columns null where the tool made no call of the family. is_tp here means
    the call finds this region: correct, and an IoU with this region of at least 0.5.
    """
    blocks, lines = [], []
    for (_st,), g in t.sort("stratum_order").group_by(["stratum"], maintain_order=True):
        g0 = g.row(0, named=True)
        cards = []
        for case, words in CASE_WORDS.items():
            e = g.filter(pl.col("case") == case)
            if e.height == 0:
                cards.append(f'<div class="fg-ex"><h4>{esc(words)}</h4><p class="fg-none">None in this bin.</p></div>')
                continue
            r0 = e.row(0, named=True)
            seq = r0["window_seq"][r0["true_start"] - r0["window_offset"]: r0["true_end"] - r0["window_offset"] + 1]
            head = (f'{esc(r0["symbol"] or r0["accession"])} (UniProt {esc(r0["accession"])}), Pfam {esc(r0["pfam_id"])}, '
                    f'residues {r0["true_start"]:,}-{r0["true_end"]:,}'
                    + (f'; best correct overlap (IoU) {r0["best_iou"]:.2f}' if r0["best_iou"] is not None else ""))
            cards.append(f'<div class="fg-ex"><h4>{esc(words)}</h4><p class="fg-exhead">{head}</p>'
                         f'{protein_line_svg(r0, e.to_dicts())}<p class="fg-res"><span class="muted">Residues of the true '
                         f'region:</span> <span class="mono">{esc(seq)}</span></p></div>')
        blocks.append(f'<details class="fg-stratum"{" open" if not blocks else ""}><summary>{esc(g0["stratum_label"])} '
                      f'<span class="muted">({g0["n_regions"]:,} regions, {g0["n_found"]:,} found by a tool)</span></summary>'
                      f'{"".join(cards)}</details>')
        m = g.filter(pl.col("case") == "median")
        if m.height:
            lines.append(f"{g0['stratum_label']} {g0['median_iou']:.2f} ({g0['n_found']:,} found)")
    leg = ('<div class="legend">'
           '<span class="item"><svg width="40" height="16" aria-hidden="true"><line x1="0" x2="40" y1="8" y2="8" '
           'stroke="var(--ink)" stroke-width="1.2"/><rect x="12" y="2" width="16" height="12" fill="var(--surface)" '
           'stroke="var(--ink)" stroke-width="1.5"/></svg>protein (line) with its true Pfam region (box)</span>'
           '<span class="item"><svg width="24" height="12" aria-hidden="true"><rect x="1" y="2" width="22" height="8" '
           'fill="var(--k-profiles)"/></svg>a tool\'s best call, finding the region, in the colour of its kind</span>'
           '<span class="item"><svg width="24" height="12" aria-hidden="true"><rect x="1.5" y="2.5" width="21" height="7" '
           'fill="none" stroke="var(--k-profiles)" stroke-width="1.5"/></svg>best call, not finding it</span>'
           '<span class="item">· no call of this family</span></div>')
    s = "Median overlap (IoU) of the best correct call, by identity bin: " + "; ".join(lines) + "."
    return leg + "".join(blocks), s


# --------------------------------------------------------------------------- Figure 17

QUAD = [("letters", 0, 0), ("profiles", 1, 0), ("shapes", 0, 1), ("kmers", 1, 1)]
QUAD_WORDS = {"letters": "top left", "profiles": "top right", "shapes": "bottom left", "kmers": "bottom right"}


def draw_proteome_map(t: pl.DataFrame, kinds: dict) -> tuple[str, str]:
    """Figure 17: every protein a tile, in four quarters, one per kind of tool: a quarter
    takes its kind's colour when a tool of that kind made a correct call on the protein,
    grey when none did.

    Table columns: accession, symbol, n_regions, letters, profiles, shapes, kmers (booleans)
    and detail (hover text: the protein and how many of its regions each tool found).
    `kinds` maps a kind to its label for the legend. An organism page can call this with
    the same columns for its own proteins.
    """
    per_row, s_, gap = 40, 16, 2
    n = t.height
    W = per_row * (s_ + gap)
    H = math.ceil(n / per_row) * (s_ + gap)
    q = s_ / 2
    t = t.with_columns(pl.sum_horizontal([pl.col(k).cast(pl.Int8) for k, _, _ in QUAD]).alias("n_kinds"))
    t = t.sort(["n_kinds", "kmers", "shapes", "profiles", "letters", "symbol", "accession"],
               descending=[True, True, True, True, True, False, False], nulls_last=True)
    p = []
    for i, r in enumerate(t.iter_rows(named=True)):
        x, y = (i % per_row) * (s_ + gap), (i // per_row) * (s_ + gap)
        cells = "".join(
            f'<rect x="{x + dx * q:.1f}" y="{y + dy * q:.1f}" width="{q:.1f}" height="{q:.1f}" '
            f'fill="{("var(--k-" + k + ")") if r[k] else "var(--null)"}"/>' for k, dx, dy in QUAD)
        p.append(f'<g data-tip="{esc(r["detail"])}" tabindex="0">{cells}</g>')
    svg = (f'<div class="figscroll"><svg class="fg" viewBox="0 0 {W} {H}" style="min-width:{W * 0.75:.0f}px;'
           f'max-width:{W * 1.4:.0f}px" role="img" aria-label="one tile per protein">{"".join(p)}</svg></div>')
    big = "".join(f'<rect x="{dx * 12}" y="{dy * 12}" width="11" height="11" fill="var(--k-{k})"/>' for k, dx, dy in QUAD)
    leg = ('<div class="legend"><span class="item"><svg width="24" height="24" aria-hidden="true">' + big
           + '</svg>one protein, four quarters</span>'
           + "".join(f'<span class="item"><svg width="12" height="12" aria-hidden="true"><rect width="11" height="11" '
                     f'fill="var(--k-{k})"/></svg>{QUAD_WORDS[k]}: {esc(kinds[k])}</span>' for k, _, _ in QUAD)
           + '<span class="item"><svg width="12" height="12" aria-hidden="true"><rect width="11" height="11" '
             'fill="var(--null)"/></svg>grey quarter: no tool of that kind found a region of the protein</span></div>')
    none = t.filter(pl.col("n_kinds") == 0).height
    allk = t.filter(pl.col("n_kinds") == len(QUAD)).height
    only_k = t.filter(pl.col("kmers") & (pl.col("n_kinds") == 1)).height
    s = (f"Of the {n:,} proteins, {allk:,} had regions found, between them, by all four kinds of tool, and "
         f"{none:,} had no region found. For {only_k:,}, the k-mer tool was the only kind that found any of its "
         "regions.")
    return leg + svg, s
