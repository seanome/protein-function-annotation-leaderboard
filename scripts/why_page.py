"""The "Why this project" page (why.html), for a reader with 30 seconds and no biology.

Every number on it comes from one of these, and is null ("—, not yet measured") when the
file or the row is missing:
  results/why/numbers.parquet          the measured counts and shares (Prompt E)
  results/why/tables/*.parquet         the per-species Swiss-Prot counts and WoRMS habitat flags
  results/why/proteome_labels.parquet  one row per human and oyster protein, with its labels
  results/why/four_proteins_calls.parquet  InterPro matches and Swiss-Prot hits on the four proteins
  data/why/four_sequences.json         their residues (scripts/why/four_sequences.py)
  data/results.parquet                 the leaderboard's own rows (the strata in Figure 4)
Nothing here is estimated.
"""

import html
import json
import math
from pathlib import Path

import polars as pl

DASH = "—"
NYM = "not yet measured"
# The proteome maps: tiles per row, the same for both species, so a tile is one protein
# at the same size in both maps and the map's height is proportional to the protein count.
COLS = 480
ROOT = Path(__file__).resolve().parents[1]

FOUR = {
    "Q25460": {
        "short": "Mussel foot protein",
        "gave": "Blue mussels use it to glue their threads to wet rock. Its repeated 10-residue "
                "stretches carry DOPA, a modified tyrosine, and chemists copy that chemistry to make "
                "glues that set under water.",
        "stratum": "repeats",
    },
    "P24856": {
        "short": "Fish antifreeze glycoprotein",
        "gave": "Antarctic rockcod make it to keep their blood from freezing in sea water below 0 °C. "
                "The chain is cut into 46 short glycopeptides, nearly all repeats of Ala-Ala-Thr, that "
                "stop ice crystals from growing.",
        "stratum": "repeats and short peptides",
    },
    "P05484": {
        "short": "Cone snail toxin (omega-conotoxin MVIIA)",
        "gave": "The 25-residue toxin of a cone snail is the painkiller ziconotide (sold as Prialt), "
                "given to patients whom opiates no longer help.",
        "stratum": "a short peptide",
    },
    "P42212": {
        "short": "Green fluorescent protein (GFP)",
        "gave": "A jellyfish protein that glows green. Labs worldwide fuse it to other proteins to "
                "watch them inside living cells; its discovery and development won the 2008 Nobel "
                "Prize in Chemistry.",
        "stratum": "a folded domain",
    },
}

CAFA_URL = "https://biofunctionprediction.org/cafa/"
OPEN_PROBLEMS_URL = "https://doi.org/10.1038/s41587-025-02694-w"
PROTEINANNOTATOR_URL = "https://nf-co.re/proteinannotator"
LAFA_URL = "https://functionbench.net/"
LAFA_REPO = "https://github.com/anphan0828/CAFA_forever"
SEANOME_URL = "https://seanome.org"
ISSUES_URL = "https://github.com/seanome/protein-function-annotation-leaderboard/issues"
CONTACT_EMAIL = "info@seanome.org"  # listed under "Contact" on seanome.org, read 2026-10-16


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def fint(x) -> str:
    return f"{int(round(x)):,}"


def nym(extra: str = "") -> str:
    return f'<span class="nym">{DASH} <small>{NYM}{extra}</small></span>'


class Why:
    """Loads everything the page reads; a missing file leaves its part as None."""

    def __init__(self, why: Path, data: Path):
        self.why = why
        self.n = {}
        if (why / "numbers.parquet").exists():
            for r in pl.read_parquet(why / "numbers.parquet").iter_rows(named=True):
                self.n[r["id"]] = r
        self.sources = json.loads((why / "sources.json").read_text()) if (why / "sources.json").exists() else {}
        self.labels = pl.read_parquet(why / "proteome_labels.parquet") if (why / "proteome_labels.parquet").exists() else None
        self.calls = (pl.read_parquet(why / "four_proteins_calls.parquet")
                      if (why / "four_proteins_calls.parquet").exists() else None)
        f = data / "why/four_sequences.json"
        self.seqs = json.loads(f.read_text()) if f.exists() else {}
        self.missing = []

    def v(self, key):
        r = self.n.get(key)
        if r is None or r["value"] is None or (isinstance(r["value"], float) and math.isnan(r["value"])):
            self.missing.append(key)
            return None
        return r["value"]

    def release(self, key):
        r = self.n.get(key)
        return r["release"] if r else None


# --------------------------------------------------------------------------- numbers

def ocean_counts(w: Why) -> dict | None:
    """Swiss-Prot entries from ocean animals, two ways (counts.py):
    strict: animal species WoRMS flags as marine, with both its freshwater and land flags
      set to 0 (sp_marine_only_animal_share). A species whose flag is empty is left out.
    broad: every animal species WoRMS flags as marine (sp_marine_animal_share), some of which
      also live in fresh water or on land.
    Recomputed here from the per-species table and checked against the stored values."""
    t = w.why / "tables"
    if not (t / "swissprot_entries_per_species_taxid.parquet").exists():
        return None
    sp = pl.read_parquet(t / "swissprot_entries_per_species_taxid.parquet")
    worms = pl.read_parquet(t / "metazoan_species_worms.parquet")
    broad = worms.filter(pl.col("isMarine") == 1)["species_name"].to_list()
    ocean = (worms.filter((pl.col("isMarine") == 1) & (pl.col("isFreshwater") == 0)
                          & (pl.col("isTerrestrial") == 0))["species_name"].to_list())
    o = sp.filter(pl.col("name").is_in(ocean))
    ob = sp.filter(pl.col("name").is_in(broad))
    total = w.v("sp_entries")
    share = w.v("sp_marine_only_animal_share")
    share_b = w.v("sp_marine_animal_share")
    n, nb = int(o["entries"].sum()), int(ob["entries"].sum())
    if total is not None and share is not None and abs(round(100 * n / total, 2) - share) > 1e-9:
        raise SystemExit(f"ocean animal entries {n} give {100 * n / total:.2f}%, numbers.parquet says {share}%")
    if w.v("sp_marine_animal_entries") is not None and nb != w.v("sp_marine_animal_entries"):
        raise SystemExit(f"broad ocean animal entries {nb}, numbers.parquet says {w.v('sp_marine_animal_entries')}")
    top = sp.sort("entries", descending=True).head(5)
    top5 = int(top["entries"].sum())
    top5_share = w.v("sp_top5_share_by_species")
    if total is not None and top5_share is not None and abs(round(100 * top5 / total, 2) - top5_share) > 1e-9:
        raise SystemExit(f"top five species hold {top5}, numbers.parquet says {top5_share}%")
    in_species = int(sp["entries"].sum())
    rest = sp.filter(~pl.col("name").is_in(top["name"].to_list() + ocean))
    return {
        "entries": n, "species": o.height, "share_pct": share, "total": total,
        "broad_entries": nb, "broad_species": ob.height, "broad_share_pct": share_b,
        "top": [(r["name"], int(r["entries"])) for r in top.iter_rows(named=True)],
        "top5": top5, "top5_share": top5_share,
        "rest_entries": int(rest["entries"].sum()), "rest_species": rest.height,
        "in_species": in_species, "no_species": int(total - in_species) if total is not None else None,
        "n_species": sp.height,
    }


def map_data(w: Why) -> dict | None:
    """One tile per protein. Lit: an InterPro match or a Swiss-Prot hit above 30% identity
    (scripts/why/maps.py, proteome_map_function_labels). Lit tiles first, then dark."""
    if w.labels is None:
        return None
    out = {"cols": COLS, "species": {}}
    for sp in ("human", "oyster"):
        # A protein's hit to its own Swiss-Prot entry does not count: most human proteins are
        # Swiss-Prot entries and would light themselves (numbers.parquet's
        # *_no_interpro_no_sprot_pct counts self-hits; it is checked below, then not shown).
        t = (w.labels.filter(pl.col("species") == sp)
             .with_columns(lit=pl.col("interpro_match") | pl.col("sprot_hit_nonself"),
                           lit_self=pl.col("interpro_match") | pl.col("sprot_hit"))
             .sort(["lit", "accession"], descending=[True, False]))
        n, k = t.height, int(t["lit"].sum())
        stored_n = w.v(f"{sp}_proteins")
        stored_dark_self = w.v(f"{sp}_no_interpro_no_sprot_pct")
        k_self = int(t["lit_self"].sum())
        if stored_n is not None and int(stored_n) != n:
            raise SystemExit(f"{sp}: {n} proteins in proteome_labels, numbers.parquet says {stored_n}")
        if stored_dark_self is not None and abs(round(100 * (n - k_self) / n, 2) - stored_dark_self) > 1e-9:
            raise SystemExit(f"{sp}: {100 * (n - k_self) / n:.2f}% dark with self-hits, numbers.parquet says {stored_dark_self}%")
        prot = [[r["accession"], r["gene"] or "", (r["protein_name"] or "")[:70],
                 int(bool(r["interpro_match"])), int(bool(r["sprot_hit_nonself"]))]
                for r in t.select("accession", "gene", "protein_name", "interpro_match", "sprot_hit_nonself").iter_rows(named=True)]
        dark_pct = round(100 * (n - k) / n, 2)
        stored_dark = w.v(f"{sp}_no_interpro_no_sprot_nonself_pct")
        if stored_dark is not None and abs(dark_pct - stored_dark) > 1e-9:
            raise SystemExit(f"{sp}: {dark_pct}% dark without self-hits, numbers.parquet says {stored_dark}%")
        out["species"][sp] = {"n": n, "lit": k, "dark": n - k, "dark_pct": dark_pct, "lit_pct": round(100 - dark_pct, 2),
                              "dark_counting_self_hits": n - k_self, "dark_pct_counting_self_hits": stored_dark_self,
                              "proteins": prot}
    return out


def strata(res: pl.DataFrame, tools: dict) -> list[dict]:
    """The leaderboard's rows behind Figure 4: Swiss-Prot features, raw calls, zebrafish as the
    relative searched; for each stratum, the tool and setting that finds the largest share of
    the known features at its own best score cut-off (fmax_recall in the tables)."""
    base = res.filter((pl.col("truth_set") == "swissprot") & (pl.col("calls") == "raw")
                      & (pl.col("species") == "zebrafish") & pl.col("fmax_recall").is_not_null())
    want = [("all", "all", "All known features"),
            ("feature_length_bin", "2-15", "Features 2 to 15 residues long"),
            ("feature_type", "REPEAT", "Repeats"),
            ("plddt", "50-70", "Features where AlphaFold's confidence is 50 to 70")]
    rows = []
    for ax, b, label in want:
        t = base.filter((pl.col("stratum_axis") == ax) & (pl.col("stratum") == b)).sort(
            ["fmax_recall", "tool", "setting"], descending=[True, False, False])
        if t.height == 0:
            rows.append({"label": label, "value": None})
            continue
        r = t.row(0, named=True)
        rows.append({"label": label, "value": r["fmax_recall"], "tool": tools[r["tool"]]["name"],
                     "tool_key": r["tool"], "setting": r["setting"], "n": int(r["n_truth_instances"]),
                     "axis": ax, "stratum": b})
    rows.append({"label": "Features where AlphaFold's confidence is below 50 (often disordered)", "value": None,
                 "why": "the leaderboard does not yet group features with confidence below 50"})
    return rows


# --------------------------------------------------------------------------- figures

def fig_concentration(oc: dict) -> str:
    """Figure 2: two panels on the same rows, reviewed proteins and reviewed proteins per species."""
    common = {"Escherichia coli": "E. coli (bacterium)", "Homo sapiens": "Human",
              "Mus musculus": "Mouse", "Arabidopsis thaliana": "Thale cress (plant)",
              "Salmonella enterica": "Salmonella (bacterium)"}
    rows = [(common.get(n, n), e, 1) for n, e in oc["top"]]
    rows.append((f"Ocean animals ({oc['species']:,} species)", oc["entries"], oc["species"]))
    # "All other species" below holds the rest, including the broad-rule ocean animals.
    rows.append((f"All other species ({oc['rest_species']:,})", oc["rest_entries"], oc["rest_species"]))

    def panel(title, unit, vals, fmt):
        W, LBL, BAR, H0, RH = 560, 200, 270, 30, 26
        vmax = max(vals)
        y = H0
        parts = [f'<text class="ptitle" x="0" y="16">{esc(title)}</text>']
        for i, ((label, _, _), v) in enumerate(zip(rows, vals)):
            if i == 5:
                y += 10  # gap between the five species and the two groups
            w_ = BAR * v / vmax
            ocean = i == 5
            parts.append(f'<text class="rlab{" strong" if ocean else ""}" x="{LBL - 8}" y="{y + 13}" text-anchor="end">{esc(label)}</text>')
            parts.append(f'<rect class="bar" x="{LBL}" y="{y + 3}" width="{max(w_, 1):.1f}" height="{RH - 8}"/>')
            parts.append(f'<text class="val" x="{LBL + max(w_, 1) + 5:.1f}" y="{y + 13}">{fmt(v)}</text>')
            y += RH
        parts.append(f'<line class="axis" x1="{LBL}" x2="{LBL}" y1="{H0}" y2="{y}"/>')
        parts.append(f'<text class="unit" x="{LBL}" y="{y + 16}">{esc(unit)}, bar length to scale from 0</text>')
        return (f'<svg viewBox="0 0 {W} {y + 24}" role="img" aria-label="{esc(title)}">' + "".join(parts) + "</svg>")

    a = panel("a. Reviewed proteins", "reviewed proteins", [r[1] for r in rows], lambda v: f"{v:,}")
    b = panel("b. Reviewed proteins per species", "reviewed proteins per species",
              [r[1] / r[2] for r in rows], lambda v: f"{v:,.0f}" if v >= 10 else f"{v:.1f}")
    return f'<div class="figpair">{a}{b}</div>'


def merge_spans(spans):
    out = []
    for s, e in sorted(spans):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def four_calls(w: Why, acc: str) -> dict | None:
    if w.calls is None:
        return None
    c = w.calls.filter(pl.col("accession") == acc)
    ipr = c.filter(pl.col("tool") == "InterPro")
    names = sorted({lab.split(" (")[0].split(" ", 1)[1] for lab in ipr["label"].to_list()
                    if lab.startswith("IPR")})
    ipr_spans = merge_spans([(r["start"], r["end"]) for r in ipr.iter_rows(named=True)])
    hits = c.filter(pl.col("tool").str.starts_with("MMseqs2") & ~pl.col("self_hit").fill_null(False)
                    & pl.col("above_30pct").fill_null(False)).sort("score", descending=True)
    hit_spans = merge_spans([(r["start"], r["end"]) for r in hits.iter_rows(named=True)])
    best = hits.row(0, named=True) if hits.height else None
    af = c.filter(pl.col("tool") == "AlphaFold DB")
    return {
        "ipr_names": names, "ipr_spans": ipr_spans, "ipr_residues": sum(e - s + 1 for s, e in ipr_spans),
        "n_hits": hits.height, "hit_spans": hit_spans, "hit_residues": sum(e - s + 1 for s, e in hit_spans),
        "best": None if best is None else {"target": best["label"].split(" ")[0], "identity": best["score"],
                                           "start": best["start"], "end": best["end"]},
        "af_mean_plddt": af["score"][0] if af.height else None,
    }


def fig_four_track(length: int, lmax: int, calls: dict | None) -> str:
    """The protein as a line at a shared residue scale, with boxes under it for where InterPro
    matches and where other Swiss-Prot proteins align (the InterPro/Pfam convention)."""
    W, X0, PW = 560, 112, 360
    sx = PW / lmax
    end = X0 + length * sx
    p = [f'<text class="tlab" x="{X0 - 8}" y="16" text-anchor="end">Protein</text>',
         f'<line class="backbone" x1="{X0}" x2="{end:.1f}" y1="12" y2="12"/>',
         f'<text class="val" x="{end + 5:.1f}" y="16">{length:,} residues</text>']
    rows = [("InterPro match", "ipr", calls["ipr_spans"] if calls else None),
            ("Swiss-Prot hit", "hit", calls["hit_spans"] if calls else None)]
    y = 30
    for lab, cls, spans in rows:
        p.append(f'<text class="tlab" x="{X0 - 8}" y="{y + 10}" text-anchor="end">{lab}</text>')
        if spans is None:
            p.append(f'<text class="none" x="{X0}" y="{y + 10}">{DASH} {NYM}</text>')
        elif not spans:
            p.append(f'<circle class="nomark" cx="{X0 + 3}" cy="{y + 6}" r="2"/>'
                     f'<text class="none" x="{X0 + 10}" y="{y + 10}">none</text>')
        for s, e in spans or []:
            p.append(f'<rect class="{cls}" x="{X0 + (s - 1) * sx:.1f}" y="{y}" width="{max((e - s + 1) * sx, 1.5):.1f}" height="12"/>')
        y += 20
    return f'<svg class="track" viewBox="0 0 {W} {y + 2}" role="img" aria-label="Where the labels fall on the protein">{"".join(p)}</svg>'


def residues_block(seq: str) -> str:
    lines = []
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60]
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        lines.append(f'<span class="pos">{i + 1:>4}</span> {groups}')
    return '<pre class="residues">' + "\n".join(lines) + "</pre>"


# --------------------------------------------------------------------------- page

def page(w: Why, res: pl.DataFrame, tools: dict, n_tools: int, up: str = "") -> tuple[str, dict]:
    oc = ocean_counts(w)
    mp = map_data(w)
    st = strata(res, tools)
    total = w.v("sp_entries")
    # The oyster and human shares come from the per-protein table with self-hits left out
    # (map_data), so the headline, the card and Figure 1 are one number.
    oy_dark = mp["species"]["oyster"]["dark_pct"] if mp else None
    oy_n = mp["species"]["oyster"]["n"] if mp else w.v("oyster_proteins")
    hu_dark = mp["species"]["human"]["dark_pct"] if mp else None
    animal = w.v("sp_metazoa_entries")

    def src_line(key, extra=""):
        r = w.n.get(key)
        if not r:
            return f'<p class="src">Source: {NYM}</p>'
        s = w.sources.get(r["source"], {})
        url = s.get("url")
        name = {"swissprot_dat": "Swiss-Prot", "worms_matches": "WoRMS",
                "sherlock_mmseqs": "MMseqs2 against Swiss-Prot"}.get(r["source"], r["source"])
        link = f'<a href="{esc(url)}">{esc(name)}</a>' if url and url.startswith("http") else esc(name)
        return f'<p class="src">{link}, {esc(r["release"])}{extra}</p>'

    # Screen 1 numbers, written once and used in the headline and the cards.
    n1 = fint(total) if total is not None else DASH
    n2 = (f"{oc['share_pct']:.2f}% to {oc['broad_share_pct']:.2f}%"
          if oc and oc["share_pct"] is not None and oc["broad_share_pct"] is not None else DASH)
    n3 = f"{oy_dark:.1f}%" if oy_dark is not None else DASH
    head3 = (f"<b>{n3}</b> of a Pacific oyster's proteins match no known family or described protein"
             if oy_dark is not None else
             "the share of a Pacific oyster's proteins that match nothing known is not yet measured")
    headline = (f'Curators have reviewed <b>{n1}</b> proteins by hand, only <b>{n2}</b> '
                f'from ocean animals, and {head3}.')

    oy_count = None
    if mp:
        oy_count = mp["species"]["oyster"]["dark"]
    card3_sub = (f"{oy_count:,} of {fint(oy_n)} have no family in InterPro (the catalogue of known protein "
                 f"families) and no other Swiss-Prot protein at least 30% identical"
                 if oy_count is not None and oy_n is not None else
                 "no InterPro family match and no described protein at least 30% identical")
    cards = f"""<div class="cards">
<div class="card"><div class="big">{n1 if total is not None else nym()}</div>
<p>proteins whose entry a curator has reviewed by hand (Swiss-Prot, the reviewed part of UniProt)</p>
{src_line("sp_entries")}</div>
<div class="card"><div class="big">{n2 if oc else nym()}</div>
<p>come from ocean animals: {oc['entries']:,} proteins from species the World Register of Marine Species (WoRMS) lists as marine only, up to {oc['broad_entries']:,} counting species that also live in fresh water or on land</p>
{src_line("sp_marine_only_animal_share")}</div>
<div class="card"><div class="big">{n3 if oy_dark is not None else nym()}</div>
<p>of a Pacific oyster's proteins get no label: {card3_sub}</p>
<p class="src">InterPro as UniProt 2026_03 lists it; MMseqs2, E-value 0.001 or less; proteome UP000005408</p></div>
</div>"""

    if mp:
        hu, oy = mp["species"]["human"], mp["species"]["oyster"]
        maps = f"""<figure class="maps" id="figure-1">
<div class="legend" aria-label="Key to Figure 1">
<span class="item"><span class="tile lit" aria-hidden="true"></span>lit: the protein has an InterPro family, or another Swiss-Prot protein at least 30% identical</span>
<span class="item"><span class="tile unlit" aria-hidden="true"></span>dark: neither</span>
<span class="item muted">one tile = one protein, {COLS} tiles per row in both maps</span></div>
<div class="mapgrid">
<div class="map"><h3>Human <small>{hu['n']:,} proteins</small></h3>
<p class="frac"><b>{hu['lit_pct']:.1f}%</b> lit, {hu['dark']:,} dark</p>
<canvas style="aspect-ratio:{COLS}/{-(-hu['n'] // COLS)}" data-species="human" aria-label="Human proteome map: {hu['lit']:,} lit tiles, {hu['dark']:,} dark"></canvas></div>
<div class="map"><h3>Pacific oyster <small>{oy['n']:,} proteins</small></h3>
<p class="frac"><b>{oy['lit_pct']:.1f}%</b> lit, {oy['dark']:,} dark</p>
<canvas style="aspect-ratio:{COLS}/{-(-oy['n'] // COLS)}" data-species="oyster" aria-label="Pacific oyster proteome map: {oy['lit']:,} lit tiles, {oy['dark']:,} dark"></canvas></div>
</div>
<div class="tip" id="maptip" hidden></div>
<figcaption><b>Figure 1.</b> Every protein of two animals, one tile each. Like Earth seen from space at night,
lit tiles show where knowledge reaches, and ocean animals are the dark continents. The oyster has {oy['dark'] / hu['dark']:.0f} times
as many dark proteins as human ({oy['dark']:,} against {hu['dark']:,}). A match to a protein's own Swiss-Prot entry does not count.
Hover or tap a tile to see the protein. Proteomes: UP000005640 and UP000005408, UniProt 2026_03, one protein per gene.</figcaption>
</figure>"""
    else:
        maps = f"""<figure class="maps" id="figure-1"><div class="placeholder">{DASH} The two proteome maps are {NYM}:
results/why/proteome_labels.parquet (Prompt E, step 4) is not in this build.</div>
<figcaption><b>Figure 1.</b> Every protein of human and Pacific oyster, one tile each, lit if it has a label. {NYM.capitalize()}.</figcaption></figure>"""

    screen1 = f"""<section class="screen s1" id="argument">
<h1 class="headline">{headline}</h1>
{cards}
{maps}
<p class="does">This site measures, for every software tool and every kind of organism, which proteins get a label
and what the run costs, and tells each lab which tool to run. Today it scores {n_tools} tools on human proteins.
The first organism page will be the Pacific oyster's; it is not built yet, so there is no link to it.</p>
</section>"""

    # Screen 2
    if oc:
        conc = f"""<figure id="figure-2">
<div class="legend" aria-label="Key to Figure 2"><span class="item"><span class="swatch bar" aria-hidden="true"></span>reviewed proteins (Swiss-Prot), counted by species</span></div>
{fig_concentration(oc)}
<figcaption><b>Figure 2.</b> The five best-studied species hold {oc['top5']:,} reviewed proteins
({oc['top5_share']:.1f}% of all), about {oc['top5'] / 5:,.0f} each. The {oc['species']:,} ocean animal species
(the low end of Card 2) hold {oc['entries']:,} together, {oc['entries'] / oc['species']:.1f} each on average.
The other {oc['rest_species']:,} species hold {oc['rest_entries']:,}. {oc['no_species']:,} entries have no species
and are left out. Human has {dict(oc['top']).get('Homo sapiens', 0):,} reviewed entries here; Figure 1 counts the {fint(w.v('human_proteins')) if w.v('human_proteins') else DASH} proteins of the
human reference proteome, one per gene, reviewed or not. Source: Swiss-Prot, {esc(w.release('sp_entries'))}; World Register of Marine Species habitat flags.</figcaption>
</figure>"""
    else:
        conc = f'<figure id="figure-2"><div class="placeholder">{DASH} {NYM}</div><figcaption><b>Figure 2.</b> Reviewed proteins by species.</figcaption></figure>'

    lmax = max((s["length"] for s in w.seqs.values()), default=1)
    blocks = []
    for acc, info in FOUR.items():
        s = w.seqs.get(acc)
        c = four_calls(w, acc)
        if s is None:
            blocks.append(f'<article class="prot"><h3>{esc(info["short"])}</h3><p>{esc(info["gave"])}</p>{nym()}</article>')
            continue
        if c is None:
            says = f"<li>{NYM}: InterPro and Swiss-Prot hits on this protein (results/why/four_proteins_calls.parquet is not in this build)</li>"
        else:
            ipr = (f"matches <i>{esc('; '.join(c['ipr_names']))}</i> on {c['ipr_residues']:,} of {s['length']:,} residues"
                   if c["ipr_spans"] else "no match")
            hit = (f"{c['n_hits']} other Swiss-Prot protein{'s' if c['n_hits'] != 1 else ''} at least 30% identical; "
                   f"the closest, {esc(c['best']['target'])}, is {int(100 * c['best']['identity'] + 0.5)}% identical over residues "
                   f"{c['best']['start']}-{c['best']['end']}" if c["best"] else "no other Swiss-Prot protein at least 30% identical")
            labelled = bool(c["ipr_spans"] or c["best"])
            verdict = ("A tool labels it today." if labelled else
                       "No tool here labels it today.")
            af = (f"<li>AlphaFold DB 3D model: mean confidence {c['af_mean_plddt']:.0f} of 100 "
                  f"(below 50 often means the chain has no fixed shape)</li>" if c["af_mean_plddt"] is not None else "")
            says = (f"<li><b>{verdict}</b></li><li>InterPro (a catalogue of protein families): {ipr}</li>"
                    f"<li>MMseqs2 against Swiss-Prot: {hit}</li>{af}"
                    f"<li>The other {n_tools - 1} tools on the leaderboard: {NYM} on this protein</li>")
        common = f", {esc(s['common_name'])}" if s.get("common_name") else ""
        blocks.append(f"""<article class="prot" id="p-{acc}">
<h3>{esc(info['short'])} <small><i>{esc(s['organism'])}</i>{common}, <a href="https://www.uniprot.org/uniprotkb/{acc}/entry">{acc}</a></small></h3>
<p class="gave">{esc(info['gave'])}</p>
{fig_four_track(s['length'], lmax, c)}
<ul class="says">{says}</ul>
<details class="res"><summary>All {s['length']:,} residues</summary>{residues_block(s['sequence'])}</details>
</article>""")
    four = f"""<figure id="figure-3">
<div class="legend" aria-label="Key to Figure 3">
<span class="item"><span class="swatch backbone" aria-hidden="true"></span>the protein, residue 1 on the left, all four at one scale</span>
<span class="item"><span class="swatch ipr" aria-hidden="true"></span>where an InterPro family matches</span>
<span class="item"><span class="swatch hit" aria-hidden="true"></span>where another Swiss-Prot protein at least 30% identical aligns</span>
<span class="item"><span class="swatch nomark" aria-hidden="true"></span>none found</span></div>
<div class="prots">{''.join(blocks)}</div>
<figcaption><b>Figure 3.</b> Four ocean proteins that changed medicine and biology, and what the tools find on
them today. Residues and features: UniProt 2026_03. Matches: the InterPro web service, and MMseqs2 against Swiss-Prot 2026_03.</figcaption>
</figure>"""

    vmax = 1.0
    srows = []
    for r in st:
        if r["value"] is None:
            srows.append(f'<tr><th scope="row">{esc(r["label"])}</th><td class="barcell"><span class="nomark" aria-hidden="true"></span> '
                         f'{DASH} {NYM}{": " + esc(r["why"]) if r.get("why") else ""}</td><td></td></tr>')
            continue
        pct = 100 * r["value"]
        srows.append(f'<tr><th scope="row">{esc(r["label"])} <small>({r["n"]:,})</small></th>'
                     f'<td class="barcell"><span class="hbar" style="width:{pct / vmax:.1f}%"></span><span class="v">{pct:.0f} in 100</span></td>'
                     f'<td class="muted">{esc(r["tool"])}</td></tr>')
    s_all = next((r for r in st if r.get("axis") == "all"), None)
    s_short = next((r for r in st if r.get("stratum") == "2-15"), None)
    s_rep = next((r for r in st if r.get("stratum") == "REPEAT"), None)
    if s_all and s_short and s_rep and s_all["value"] is not None:
        tie = (f"The proteins that gave us underwater glue, antifreeze, a painkiller and the green light of cell biology "
               f"were all found at the bench, by people, not by a tool. The three that tools label today are found by "
               f"their likeness to proteins people described first; the antifreeze protein still gets no label. "
               f"Three of the four are built from repeats or short peptides. On the leaderboard, the best tool finds {100 * s_short['value']:.0f} in 100 "
               f"known features 2 to 15 residues long, against {100 * s_all['value']:.0f} in 100 features of any length. "
               f"Repeats are found about as often as other features ({100 * s_rep['value']:.0f} in 100), and stretches with no "
               f"fixed shape are {NYM}.")
    else:
        tie = NYM
    strata_fig = f"""<figure id="figure-4">
<div class="legend" aria-label="Key to Figure 4"><span class="item"><span class="swatch bar" aria-hidden="true"></span>share of known features the best tool finds, out of 100, bar to scale from 0 to 100</span>
<span class="item"><span class="swatch nomark" aria-hidden="true"></span>{NYM}</span></div>
<div class="tablewrap"><table class="strata"><thead><tr><th>Known Swiss-Prot features on human proteins (how many)</th><th>Found, out of 100</th><th>Best tool</th></tr></thead>
<tbody>{''.join(srows)}</tbody></table></div>
<figcaption><b>Figure 4.</b> For each group of known features (Swiss-Prot annotations on human proteins), the
share found by the tool that does best on that group. Each tool uses the score cut-off that suits it best. Picking the best tool for each row raises every row. Each tool
labelled human proteins by searching the proteins of zebrafish (a freshwater fish), the relative in the leaderboard's <a href="{up}index.html">default view</a>. Source: the leaderboard tables (footer).</figcaption>
</figure>"""

    screen2 = f"""<section class="screen" id="dark">
<h2>Why ocean animals are dark</h2>
<p class="prose">Science describes proteins where people work on them. Five species hold more reviewed proteins than all ocean animals put together, many times over.</p>
{conc}
{four}
<p class="prose tie">{tie}</p>
{strata_fig}
</section>"""

    # Screen 3
    units = [
        ("Organism pages for ocean animals, each measured on its own proteins", nym(), "built so far"),
        ("Cost per page in CPU-hours (processor hours)", nym(), "to be measured on the first three pages"),
        ("Cost per page in curator hours (a person checking the page)", nym(), "to be measured on the first three pages"),
        ("Tools scored on each page", nym(),
         f"the leaderboard scores {n_tools} tools today, on human proteins only"),
        ("Copy command button on each page", nym(),
         "one click copies the command that runs the best tool on that organism's proteins, so a lab can run it the same day; not built yet"),
        ("Lab results sent in: a protein and what an experiment showed it does", nym(), "counted through a way to send results in that is not built yet"),
        ("Labs that used a page", nym(), "counted by clicks on the Copy command button, with no record of who clicked"),
    ]
    trows = "".join(f'<tr><th scope="row">{esc(a)}</th><td class="num">{b}</td><td>{esc(c)}</td></tr>' for a, b, c in units)
    screen3 = f"""<section class="screen" id="funding">
<h2>What funding buys</h2>
<p class="prose">Each unit below is a count anyone can check on this site. None has been measured yet, so each shows a dash.</p>
<div class="tablewrap"><table class="units"><thead><tr><th>Unit</th><th>Count today</th><th>How it is counted</th></tr></thead>
<tbody>{trows}</tbody></table></div>
<p class="ask"><b>The ask:</b> fund the first 100 ocean animal pages, at the cost per page measured on the first three ({DASH} {NYM}).</p>
<p class="prose">For scale: <a href="{CAFA_URL}">CAFA</a>, the Critical Assessment of protein Function Annotation,
has run as a community challenge since 2010. <a href="{OPEN_PROBLEMS_URL}">Open Problems</a> (Luecken et al., Nature
Biotechnology 2025, co-authored by Olga Botvinnik) does the same for single-cell analysis tools.</p>
</section>"""

    screen4 = f"""<section class="screen" id="who">
<h2>Who is behind it</h2>
<ul class="who">
<li><a href="{SEANOME_URL}">Seanome Inc.</a>, a nonprofit.</li>
<li>The maintainers of <a href="{PROTEINANNOTATOR_URL}">nf-core/proteinannotator</a>, the pipeline that labels a proteome's proteins.</li>
<li>If the proposed fork goes ahead (nothing agreed yet): the LAFA / CAFA group, who run <a href="{LAFA_URL}">functionbench.net</a> (<a href="{LAFA_REPO}">source</a>).</li>
</ul>
<p><b>Licence:</b> the site's own repository has no licence file yet. Each tool keeps its own licence, listed on its <a href="{up}tools/index.html">tool card</a>.</p>
<p><b>Contact:</b> open an issue at <a href="{ISSUES_URL}">github.com/seanome/protein-function-annotation-leaderboard/issues</a>, or write to <span class="mono sel">{CONTACT_EMAIL}</span>.</p>
</section>"""

    used = sorted({k for k in w.n})
    page_src = ("swissprot_dat", "worms_matches", "human_canonical_fasta", "oyster_canonical_fasta",
                "human_interpro_xrefs", "oyster_interpro_xrefs", "sherlock_mmseqs", "afdb_predictions")
    srcs = sorted({k for k in w.sources if k in page_src or k.startswith(("uniprot_", "interpro_matches_", "afdb_confidence_"))}
                  - {"uniprot_relnotes"})
    src_items = "".join(
        f'<li><a href="{esc(w.sources[k]["url"])}">{esc(w.sources[k].get("file", k))}</a>, '
        f'{esc(w.sources[k].get("release", ""))}, fetched {esc(str(w.sources[k].get("fetched_utc", ""))[:10])}</li>'
        for k in srcs if str(w.sources[k].get("url", "")).startswith("http"))
    sources = f"""<section class="screen" id="sources"><details><summary>The files behind these numbers</summary>
<p class="prose">Counts: <span class="mono">results/why/numbers.parquet</span> and <span class="mono">results/why/proteome_labels.parquet</span>, built by <span class="mono">scripts/why/</span> from the files below.
Missing parts: {esc(", ".join(sorted(set(w.missing))) or "none")}.</p><ul class="srclist">{src_items}</ul></details></section>"""

    body = screen1 + screen2 + screen3 + screen4 + sources
    data = {
        "headline": html.unescape(headline.replace("<b>", "").replace("</b>", "")),
        "numbers": {
            "swissprot_entries": total, "ocean_animal_entries": oc["entries"] if oc else None,
            "ocean_animal_species": oc["species"] if oc else None,
            "ocean_animal_share_pct": oc["share_pct"] if oc else None,
            "top5_species": oc["top"] if oc else None, "top5_share_pct": oc["top5_share"] if oc else None,
            "oyster_no_function_label_pct": oy_dark, "human_no_function_label_pct": hu_dark,
            "oyster_proteins": oy_n,
        },
        "maps": None if mp is None else {sp: {k: v for k, v in d.items() if k != "proteins"} for sp, d in mp["species"].items()},
        "four_proteins": {acc: four_calls(w, acc) for acc in FOUR},
        "strata": st,
        "not_yet_measured": ["organism pages built", "cost per page (CPU-hours, curator hours)", "tools scored per page",
                             "wet-lab truth rows contributed", "labs that used a page"] + sorted(set(w.missing)),
    }
    return body, data, mp
