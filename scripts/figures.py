"""Draw functions for the leaderboard's figures 11 to 17.

Each function takes one table (a polars DataFrame built by a data step in build_site.py)
and returns the figure as HTML (a legend, then an inline SVG) and a one-line sentence built
from the same rows it draws. No number on a figure or in its sentence is typed here: every
value comes from the table. Hover text sits in data-tip attributes; assets/figures.js shows
it on hover and on tap.

Colours are the site's kind-of-method colours (--k-letters, --k-profiles, --k-shapes,
--k-kmers, --k-models): a tool keeps its kind's colour in every figure, and tools of one
kind differ by marker shape, as in Table 1 and Figure 3.
"""

from __future__ import annotations

import html
import math

import polars as pl


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def fnum(v: float, d: int = 3) -> str:
    return f"{v:.{d}f}"


def text_w(s: str, px: float = 6.6) -> float:
    """Width of a 12 px Source Sans 3 label, about 6.6 px per character."""
    return len(s) * px


def marker(shape: str, x: float, y: float, color: str, r: float = 4.5, filled: bool = True,
           tip: str = "", cls: str = "") -> str:
    """One point mark. Shapes match svg_marker() in build_site.py."""
    fill = color if filled else "var(--surface)"
    a = f' data-tip="{esc(tip)}" tabindex="0"' if tip else ""
    c = f' class="{cls}"' if cls else ""
    st = f'fill="{fill}" stroke="{color}" stroke-width="1.5"'
    if shape == "square":
        return f'<rect x="{x - r:.1f}" y="{y - r:.1f}" width="{2 * r:.1f}" height="{2 * r:.1f}" {st}{a}{c}/>'
    if shape == "diamond":
        k = r * 1.3
        return f'<path d="M{x:.1f},{y - k:.1f} L{x + k:.1f},{y:.1f} L{x:.1f},{y + k:.1f} L{x - k:.1f},{y:.1f} Z" {st}{a}{c}/>'
    if shape == "triangle":
        k = r * 1.2
        return f'<path d="M{x:.1f},{y - k:.1f} L{x + k:.1f},{y + k * .8:.1f} L{x - k:.1f},{y + k * .8:.1f} Z" {st}{a}{c}/>'
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" {st}{a}{c}/>'


def legend_mark(shape: str, color: str, line: bool = True, pale: bool = False) -> str:
    """A legend swatch: a short line with the marker on it, as the series is drawn."""
    op = ' opacity=".5"' if pale else ""
    ln = (f'<line x1="1" x2="27" y1="8" y2="8" stroke="{color}" stroke-width="{1 if pale else 2}"{op}/>'
          if line else "")
    mk = "" if pale else marker(shape, 14, 8, color, r=4)
    return f'<svg width="28" height="16" viewBox="0 0 28 16" aria-hidden="true">{ln}{mk}</svg>'


def tool_legend(t: pl.DataFrame, pale_note: str | None) -> str:
    """Legend of the full-strength series of a table with tool, name, kind, shape and
    highlight columns, in table order, plus the pale kmerseek lines when there are any."""
    seen, items = set(), []
    for r in t.filter(pl.col("highlight")).iter_rows(named=True):
        key = (r["tool"], r["setting"])
        if key in seen:
            continue
        seen.add(key)
        lab = r["name"] + (f", {r['label']}" if r["tool"] == "kmerseek" else "")
        items.append(f'<span class="item">{legend_mark(r["shape"], "var(--k-" + r["kind"] + ")")}{esc(lab)}</span>')
    if pale_note:
        items.append(f'<span class="item">{legend_mark("circle", "var(--k-kmers)", pale=True)}{esc(pale_note)}</span>')
    return '<div class="legend">' + "".join(items) + "</div>"


def place_labels(ys: list[float], gap: float, lo: float, hi: float) -> list[float]:
    """Move end-of-line labels apart so no two are closer than gap, keeping order."""
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    out = [0.0] * len(ys)
    last = lo - gap
    for i in order:
        y = max(ys[i], last + gap)
        out[i] = y
        last = y
    over = last - hi
    if over > 0:
        for i in order:
            out[i] -= over
    return out


# --------------------------------------------------------------------------- Figure 11

def draw_distance(t: pl.DataFrame) -> tuple[str, str]:
    """Figure 11: Fmax against the divergence time of the relative searched.

    Table columns: tool, setting, name, label, kind, shape, highlight, organism, mya, fmax,
    n. One line per tool and setting, one dot per relative searched.
    """
    W, H, L, R, T, B = 760, 414, 52, 168, 58, 50
    mya = sorted(set(t["mya"].to_list()))
    lo, hi = math.log10(min(mya) * 0.8), math.log10(max(mya) * 1.1)
    px = lambda m: L + (math.log10(m) - lo) / (hi - lo) * (W - L - R)  # noqa: E731
    ymax = math.ceil(t["fmax"].max() * 10) / 10
    py = lambda v: T + (1 - v / ymax) * (H - T - B)  # noqa: E731
    p = []
    # y grid from the ticks
    step = 0.1 if ymax <= 0.6 else 0.2
    for k in range(int(round(ymax / step)) + 1):
        v = k * step
        p.append(f'<line x1="{L}" x2="{W - R}" y1="{py(v):.1f}" y2="{py(v):.1f}" class="fg-grid"/>')
        p.append(f'<text x="{L - 6}" y="{py(v) + 4:.1f}" text-anchor="end" class="fg-tick">{v:.1f}</text>')
    p.append(f'<text x="14" y="{(T + H - B) / 2:.0f}" transform="rotate(-90 14 {(T + H - B) / 2:.0f})" '
             f'text-anchor="middle" class="fg-ax">Fmax (higher is better)</text>')
    # x ticks: each decade and its 3x step inside the range
    ticks = [v for e in range(0, 5) for v in (10 ** e, 3 * 10 ** e) if min(mya) * 0.8 <= v <= max(mya) * 1.1]
    for v in ticks:
        p.append(f'<line x1="{px(v):.1f}" x2="{px(v):.1f}" y1="{H - B}" y2="{H - B + 4}" class="fg-axis"/>')
        p.append(f'<text x="{px(v):.1f}" y="{H - B + 17}" text-anchor="middle" class="fg-tick">{v:,}</text>')
    p.append(f'<line x1="{L}" x2="{W - R}" y1="{H - B}" y2="{H - B}" class="fg-axis"/>')
    p.append(f'<text x="{(L + W - R) / 2:.0f}" y="{H - 10}" text-anchor="middle" class="fg-ax">'
             'million years since human and the relative shared an ancestor (Mya)</text>')
    # relatives named along the top in two rows, where their names fit; every one is in the
    # hover text of its dots
    # shortest names are placed first, so more of them fit
    orgs = (t.select("organism", "mya").unique()
            .sort(pl.col("organism").str.split(" (").list.first().str.len_chars(), "mya", "organism"))
    placed: list[list[tuple[float, float]]] = [[], []]
    for o in orgs.iter_rows(named=True):
        x = px(o["mya"])
        name = o["organism"].split(" (")[0]
        w = text_w(name, 6)
        a_, b_ = x - w / 2 - 4, x + w / 2 + 4
        row = next((i for i in (0, 1) if b_ < W and all(b_ < c or a_ > d for c, d in placed[i])), None)
        p.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{T - 6 - (14 if row == 1 else 0)}" y2="{H - B}" class="fg-rel"/>')
        if row is not None:
            p.append(f'<text x="{x:.1f}" y="{T - 12 - 14 * row}" text-anchor="middle" class="fg-small">{esc(name)}</text>')
            placed[row].append((a_, b_))
    # pale kmerseek settings first, full-strength lines on top
    ends = []
    for hl in (False, True):
        for (tool, setting), g in t.filter(pl.col("highlight") == hl).group_by(["tool", "setting"], maintain_order=True):
            # relatives at the same Mya (cow and dog, 89) are joined lowest Fmax first
            g = g.sort("mya", "fmax")
            r0 = g.row(0, named=True)
            col = f"var(--k-{r0['kind']})"
            pts = " ".join(f"{px(m):.1f},{py(v):.1f}" for m, v in zip(g["mya"], g["fmax"]))
            if not hl:
                p.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1" opacity=".5"/>')
                continue
            p.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="2"/>')
            for r in g.iter_rows(named=True):
                tip = (f"{r['name']}, {r['label']}: {r['organism']} ({r['mya']:,.0f} Mya), "
                       f"Fmax {r['fmax']:.3f}, n = {r['n']:,} true domains")
                p.append(marker(r["shape"], px(r["mya"]), py(r["fmax"]), col, tip=tip))
            last = g.row(g.height - 1, named=True)
            ends.append((py(last["fmax"]), px(last["mya"]), r0["name"], col))
    ys = place_labels([e[0] for e in ends], 13, T, H - B)
    for (y0, x0, name, col), y in zip(ends, ys):
        p.append(f'<line x1="{x0 + 7:.1f}" x2="{W - R + 14:.1f}" y1="{y0:.1f}" y2="{y:.1f}" stroke="{col}" stroke-width=".8"/>')
        p.append(f'<text x="{W - R + 18:.1f}" y="{y + 4:.1f}" class="fg-lab">{esc(name)}</text>')
    svg = (f'<div class="figscroll"><svg class="fg" viewBox="0 0 {W} {H}" style="min-width:600px;max-width:950px" role="img" '
           f'aria-label="Fmax against divergence time, one line per tool and setting">{"".join(p)}</svg></div>')
    kmer_pale = t.filter(~pl.col("highlight"))
    note = (f"kmerseek, its other {kmer_pale.select('setting').n_unique()} settings" if kmer_pale.height else None)
    # sentence: best tool at the nearest and the farthest relative, and kmerseek's drawn setting
    near, far = min(mya), max(mya)
    bn = t.filter(pl.col("highlight") & (pl.col("mya") == near)).sort(["fmax", "tool"], descending=[True, False]).row(0, named=True)
    bf = t.filter(pl.col("highlight") & (pl.col("mya") == far)).sort(["fmax", "tool"], descending=[True, False]).row(0, named=True)
    ks = t.filter((pl.col("tool") == "kmerseek") & pl.col("highlight")).sort("mya")
    hl = t.filter(pl.col("highlight"))
    a_near = hl.filter(pl.col("mya") == near).select("tool", "setting", pl.col("fmax").alias("f_near"))
    a_far = hl.filter(pl.col("mya") == far).select("tool", "setting", pl.col("fmax").alias("f_far"))
    both = a_near.join(a_far, on=["tool", "setting"], how="inner")
    lower = both.filter(pl.col("f_far") < pl.col("f_near")).height
    missing = []
    for r in hl.select("tool", "setting", "name").unique(maintain_order=True).iter_rows(named=True):
        have = hl.filter((pl.col("tool") == r["tool"]) & (pl.col("setting") == r["setting"]))["mya"].to_list()
        gone = [bn["organism"] if near not in have else None, bf["organism"] if far not in have else None]
        gone = [g_ for g_ in gone if g_]
        if gone:
            missing.append(f"{r['name']} has no {' or '.join(gone)} search")
    nm = lambda r: r["name"] + (f" ({r['label']})" if r["tool"] == "kmerseek" else "")  # noqa: E731
    who = "Every tool drawn in full" if lower == both.height else f"{lower} of the {both.height} tools drawn in full"
    s = (f"{who} scores lower at {bf['organism']} ({far:,.0f} Mya) than at {bn['organism']} ({near:,.0f} Mya)"
         + (f" ({'; '.join(missing)})" if missing else "") + ". "
         f"At {bn['organism']} the best is {nm(bn)} at Fmax {bn['fmax']:.3f}; at {bf['organism']} it is {nm(bf)} at "
         f"{bf['fmax']:.3f}.")
    if ks.height:
        k0, k1 = ks.row(0, named=True), ks.row(ks.height - 1, named=True)
        s += (f" kmerseek at the setting on its tool card ({k0['label']}) goes from {k0['fmax']:.3f} to "
              f"{k1['fmax']:.3f}.")
    return tool_legend(t, note) + svg, s


# --------------------------------------------------------------------------- Figure 12

def draw_identity(t: pl.DataFrame) -> tuple[str, str]:
    """Figure 12: Fmax in each identity bin, one line per tool and setting.

    Table columns: tool, setting, name, label, kind, shape, highlight, bin, bin_label,
    bin_order, joined (False for the bin outside the identity scale), fmax, n, enough (n large
    enough for the sentence to name the bin).
    """
    bins = t.select("bin", "bin_label", "bin_order", "joined", "n", "enough").unique().sort("bin_order")
    W, H, L, R, T, B = 700, 340, 52, 168, 16, 62
    nb = bins.height
    slot = (W - L - R) / (nb + 0.6)
    xs = {}
    for i, b in enumerate(bins.iter_rows(named=True)):
        xs[b["bin"]] = L + slot * (i + 0.5 + (0.6 if not b["joined"] else 0))
    ymax = math.ceil(t["fmax"].max() * 10) / 10
    py = lambda v: T + (1 - v / ymax) * (H - T - B)  # noqa: E731
    p = []
    step = 0.1 if ymax <= 0.6 else 0.2
    for k in range(int(round(ymax / step)) + 1):
        v = k * step
        p.append(f'<line x1="{L}" x2="{W - R}" y1="{py(v):.1f}" y2="{py(v):.1f}" class="fg-grid"/>')
        p.append(f'<text x="{L - 6}" y="{py(v) + 4:.1f}" text-anchor="end" class="fg-tick">{v:.1f}</text>')
    p.append(f'<text x="14" y="{(T + H - B) / 2:.0f}" transform="rotate(-90 14 {(T + H - B) / 2:.0f})" '
             f'text-anchor="middle" class="fg-ax">Fmax (higher is better)</text>')
    p.append(f'<line x1="{L}" x2="{W - R}" y1="{H - B}" y2="{H - B}" class="fg-axis"/>')
    for b in bins.iter_rows(named=True):
        x = xs[b["bin"]]
        p.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{T}" y2="{H - B}" class="fg-rel"/>')
        words = b["bin_label"].split(" ")
        lines = [b["bin_label"]] if text_w(b["bin_label"]) < slot else [" ".join(words[:len(words) // 2]), " ".join(words[len(words) // 2:])]
        for j, ln in enumerate(lines):
            p.append(f'<text x="{x:.1f}" y="{H - B + 16 + 13 * j}" text-anchor="middle" class="fg-tick">{esc(ln)}</text>')
        p.append(f'<text x="{x:.1f}" y="{H - B + 18 + 13 * len(lines)}" text-anchor="middle" class="fg-small">n = {b["n"]:,}</text>')
    p.append(f'<text x="{(L + W - R) / 2:.0f}" y="{H - 2}" text-anchor="middle" class="fg-ax">'
             '% identity of the closest same-family domain in the relative (n = true domains)</text>')
    ends = []
    for hl in (False, True):
        for _, g in t.filter(pl.col("highlight") == hl).group_by(["tool", "setting"], maintain_order=True):
            g = g.sort("bin_order")
            r0 = g.row(0, named=True)
            col = f"var(--k-{r0['kind']})"
            j = g.filter(pl.col("joined"))
            pts = " ".join(f"{xs[b]:.1f},{py(v):.1f}" for b, v in zip(j["bin"], j["fmax"]))
            if not hl:
                p.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1" opacity=".5"/>')
                continue
            p.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="2"/>')
            for r in g.iter_rows(named=True):
                tip = f"{r['name']}, {r['label']}: {r['bin_label']}, Fmax {r['fmax']:.3f}, n = {r['n']:,} true domains"
                p.append(marker(r["shape"], xs[r["bin"]], py(r["fmax"]), col, filled=r["enough"], tip=tip))
            lastj = j.row(j.height - 1, named=True)
            ends.append((py(lastj["fmax"]), xs[lastj["bin"]], r0["name"], col))
    ys = place_labels([e[0] for e in ends], 13, T, H - B)
    for (y0, x0, name, col), y in zip(ends, ys):
        x1 = W - R + 14
        p.append(f'<line x1="{x0 + 7:.1f}" x2="{x1:.1f}" y1="{y0:.1f}" y2="{y:.1f}" stroke="{col}" stroke-width=".8" opacity=".6"/>')
        p.append(f'<text x="{x1 + 4:.1f}" y="{y + 4:.1f}" class="fg-lab">{esc(name)}</text>')
    svg = (f'<div class="figscroll"><svg class="fg" viewBox="0 0 {W} {H}" style="min-width:560px;max-width:860px" role="img" '
           f'aria-label="Fmax by identity bin, one line per tool">{"".join(p)}</svg></div>')
    pale = t.filter(~pl.col("highlight"))
    note = f"kmerseek, its other {pale.select('setting').n_unique()} settings" if pale.height else None
    # the sentence reads the lowest and highest identity bins with enough true domains
    jb = bins.filter(pl.col("joined") & pl.col("enough"))
    first, last = jb.row(0, named=True), jb.row(jb.height - 1, named=True)
    nm = lambda r: r["name"] + (f" ({r['label']})" if r["tool"] == "kmerseek" else "")  # noqa: E731
    top = lambda b_: t.filter(pl.col("highlight") & (pl.col("bin") == b_)).sort(["fmax", "tool", "setting"], descending=[True, False, False]).row(0, named=True)  # noqa: E731
    a, b = top(first["bin"]), top(last["bin"])
    s = (f"At {first['bin_label']} identity the best is {nm(a)} at Fmax {a['fmax']:.3f} "
         f"(n = {first['n']:,}); at {last['bin_label']} it is {nm(b)} at {b['fmax']:.3f} (n = {last['n']:,}).")
    nh = bins.filter(~pl.col("joined"))
    if nh.height:
        c = top(nh["bin"][0])
        s += (f" Where the relative has no domain of the same family (n = {nh['n'][0]:,}), the best is "
              f"{nm(c)} at {c['fmax']:.3f}.")
    leg = tool_legend(t, note)
    if bins.filter(~pl.col("enough")).height:
        open_mark = marker("circle", 14, 8, "var(--muted)", r=4, filled=False)
        leg = leg[:-len("</div>")] + (
            f'<span class="item"><svg width="28" height="16" aria-hidden="true">{open_mark}</svg>open mark: under '
            f'{t["min_n"][0]} true domains in the bin, too few to rank tools</span></div>')
    return leg + svg, s


# --------------------------------------------------------------------------- Figure 14

def draw_feature_types(t: pl.DataFrame) -> tuple[str, str]:
    """Figure 14: one bar chart per Swiss-Prot feature type, tools as bars.

    Table columns: panel, panel_label, panel_order, measured (False: the tables have no such
    bin), tool, setting, name, label, kind, fmax, n, n_settings (settings the bar is the best of). One row per panel and tool; a panel
    that is not measured has one row with tool null.
    """
    ms = t.filter(pl.col("measured"))
    xmax = math.ceil(ms["fmax"].max() * 10) / 10
    W, L, BW = 390, 164, 112
    px = lambda v: L + v / xmax * BW  # noqa: E731
    panels = []
    for (panel,), g in t.sort("panel_order").group_by(["panel"], maintain_order=True):
        r0 = g.row(0, named=True)
        if not r0["measured"]:
            panels.append(f'<div class="fg-panel"><h4>{esc(r0["panel_label"])}</h4>'
                          f'<p class="fg-none">Not drawn: the results table has no {esc(r0["panel_label"])} subset.</p></div>')
            continue
        n = r0["n"]
        rows = g.height
        H = 26 + rows * 22 + 26
        p = []
        for k in range(int(round(xmax / 0.2)) + 1):
            v = k * 0.2
            p.append(f'<line x1="{px(v):.1f}" x2="{px(v):.1f}" y1="8" y2="{H - 26}" class="fg-grid"/>')
            p.append(f'<text x="{px(v):.1f}" y="{H - 12}" text-anchor="middle" class="fg-tick">{v:.1f}</text>')
        for i, r in enumerate(g.iter_rows(named=True)):
            y = 14 + i * 22
            lab = r["name"] if r["tool"] != "kmerseek" else f"kmerseek, {r['label']}"
            short = lab if text_w(lab, 6) < L - 8 else lab[: int((L - 14) / 6)] + "…"
            tip = f"{lab}: Fmax {r['fmax']:.3f}, n = {r['n']:,} true {r0['panel_label']} features"
            p.append(f'<text x="{L - 6}" y="{y + 11}" text-anchor="end" class="fg-small" data-tip="{esc(tip)}">{esc(short)}</text>')
            p.append(f'<rect x="{L}" y="{y + 2}" width="{max(px(r["fmax"]) - L, 0.8):.1f}" height="14" '
                     f'fill="var(--k-{r["kind"]})" data-tip="{esc(tip)}" tabindex="0"/>')
            p.append(f'<text x="{px(r["fmax"]) + 4:.1f}" y="{y + 13}" class="fg-small">{r["fmax"]:.2f} (n = {r["n"]:,})</text>')
        p.append(f'<text x="{L + BW / 2:.0f}" y="{H}" text-anchor="middle" class="fg-small">Fmax (higher is better)</text>')
        panels.append(f'<div class="fg-panel"><h4>{esc(r0["panel_label"])} <span class="muted">n = {n:,}</span></h4>'
                      f'<svg class="fg" viewBox="0 0 {W} {H}" role="img" aria-label="Fmax per tool, {esc(r0["panel_label"])}">'
                      f'{"".join(p)}</svg></div>')
    best = (ms.sort(["fmax", "tool"], descending=[True, False]).group_by("panel", maintain_order=True).first()
            .sort("panel_order"))
    nm = lambda r: r["name"] + (f" ({r['label']})" if r["tool"] == "kmerseek" else "")  # noqa: E731
    parts = [f"{r['panel_label']}: {nm(r)} {r['fmax']:.2f} (n = {r['n']:,})" for r in best.iter_rows(named=True)]
    missing = t.filter(~pl.col("measured"))["panel_label"].to_list()
    s = "Best tool per feature type: " + "; ".join(parts) + "."
    ns = ms.filter(pl.col("tool") == "kmerseek")["n_settings"]
    if ns.len() and ns.max() > 1:
        s += (f" kmerseek's bar is the best of its {ns.max()} settings in each panel, which favours it; every other "
              "tool ran one setting.")
    if missing:
        s += " Not in the tables: " + ", ".join(missing) + "."
    return '<div class="fg-grid-panels">' + "".join(panels) + "</div>", s


# --------------------------------------------------------------------------- Figure 15

def draw_pr_points(t: pl.DataFrame) -> tuple[str, str]:
    """Figure 15 (tool card): precision against recall at the Fmax cutoff, one point per
    subset of the proteins.

    Table columns: name, label, kind, group ("all", "identity", "length"), stratum_label,
    precision, recall, fmax, n.
    """
    W, H, L, B, T, R = 340, 330, 46, 44, 12, 12
    px = lambda v: L + v * (W - L - R)  # noqa: E731
    py = lambda v: T + (1 - v) * (H - T - B)  # noqa: E731
    col = f"var(--k-{t['kind'][0]})"
    p = []
    for k in range(6):
        v = k / 5
        p.append(f'<line x1="{px(v):.1f}" x2="{px(v):.1f}" y1="{T}" y2="{H - B}" class="fg-grid"/>')
        p.append(f'<line x1="{L}" x2="{W - R}" y1="{py(v):.1f}" y2="{py(v):.1f}" class="fg-grid"/>')
        p.append(f'<text x="{px(v):.1f}" y="{H - B + 15}" text-anchor="middle" class="fg-tick">{v:.1f}</text>')
        if k:  # 0.0 is written once, under the x axis
            p.append(f'<text x="{L - 5}" y="{py(v) + 4:.1f}" text-anchor="end" class="fg-tick">{v:.1f}</text>')
    p.append(f'<text x="{(L + W - R) / 2:.0f}" y="{H - 8}" text-anchor="middle" class="fg-ax">recall at the Fmax cutoff</text>')
    p.append(f'<text x="12" y="{(T + H - B) / 2:.0f}" transform="rotate(-90 12 {(T + H - B) / 2:.0f})" '
             f'text-anchor="middle" class="fg-ax">precision at the Fmax cutoff</text>')
    a = t.filter(pl.col("group") == "all").row(0, named=True)
    F = a["fmax"]
    # every precision-recall pair with the same F as the whole set
    pts = []
    for i in range(200):
        rc = F / 2 + (1 - F / 2) * (i + 1) / 200
        pr = F * rc / (2 * rc - F)
        if pr <= 1:
            pts.append(f"{px(rc):.1f},{py(pr):.1f}")
    p.append(f'<polyline points="{" ".join(pts)}" fill="none" class="fg-iso"/>')
    shapes = {"identity": "circle", "length": "square", "all": "diamond"}
    for r in t.sort("group").iter_rows(named=True):
        tip = (f"{r['stratum_label']}: precision {r['precision']:.3f}, recall {r['recall']:.3f}, "
               f"Fmax {r['fmax']:.3f}, n = {r['n']:,} true domains")
        big = r["group"] == "all"
        p.append(marker(shapes[r["group"]], px(r["recall"]), py(r["precision"]), col, r=6 if big else 4,
                        filled=big or r["group"] == "identity", tip=tip))
    svg = (f'<svg class="fg fg-sq" viewBox="0 0 {W} {H}" role="img" aria-label="precision and recall at the Fmax '
           f'cutoff for each subset">{"".join(p)}</svg>')
    leg = ('<div class="legend">'
           f'<span class="item">{legend_mark("diamond", col, line=False)}all proteins</span>'
           f'<span class="item">{legend_mark("circle", col, line=False)}one identity bin</span>'
           f'<span class="item"><svg width="28" height="16" aria-hidden="true">{marker("square", 14, 8, col, r=4, filled=False)}</svg>one region-length bin</span>'
           '<span class="item"><svg width="28" height="16" aria-hidden="true"><line x1="1" x2="27" y1="8" y2="8" class="fg-iso"/></svg>'
           'every precision and recall with the same Fmax as all proteins</span></div>')
    sub = t.filter(pl.col("group") != "all")
    hi = sub.sort(["fmax", "stratum_label"], descending=[True, False]).row(0, named=True) if sub.height else None
    lo = sub.sort(["fmax", "stratum_label"]).row(0, named=True) if sub.height else None
    s = (f"On all proteins {a['name']} reaches precision {a['precision']:.3f} and recall {a['recall']:.3f} at its "
         f"Fmax cutoff (Fmax {F:.3f}, n = {a['n']:,}).")
    if hi and lo:
        s += f" Its best subset is {hi['stratum_label']} (Fmax {hi['fmax']:.3f}, n = {hi['n']:,})"
        tied = sub.filter(pl.col("fmax") == lo["fmax"])
        if tied.height > 1:
            s += f"; its lowest Fmax, {lo['fmax']:.3f}, is shared by {tied.height} subsets."
        else:
            s += f", its worst {lo['stratum_label']} ({lo['fmax']:.3f}, n = {lo['n']:,})."
    return leg + svg, s
