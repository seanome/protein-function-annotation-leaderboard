"""Proteome maps. One tile per protein, the same tile size and 150 tiles per row for both
species.

    python3 -I scripts/why/maps.py
Writes results/why/maps/proteome_map_human_oyster_interpro_or_swissprot_uniprot_2026_03.{svg,png,pdf}:
lit when the protein has an InterPro match or a Swiss-Prot hit above 30% identity, dark
otherwise. And proteome_map_human_oyster_any_label_uniprot_2026_03.{svg,png,pdf}: the same
lit and dark, with the dark tiles split by whether AlphaFold DB has a model.
"""

import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pubfig as pf  # noqa: E402

OUT = pathlib.Path("results/why/maps")
COLS = 150  # tiles per row, both species
# Seanome style guide colours, chosen by Olga on 2026-10-16: Soft Teal for lit, Seanome Blue
# for dark. Lightness carries the meaning, so the maps also read in greyscale.
LIT, DARK = "#89C2BA", "#1E3D6B"
SPECIES = {"human": "Human (Homo sapiens)", "oyster": "Pacific oyster (Crassostrea gigas)"}
LABEL_SETS = {
    "proteome_map_human_oyster_interpro_or_swissprot_uniprot_2026_03": (
        ["interpro_match", "sprot_hit"],
        "an InterPro match or a Swiss-Prot hit above 30% identity (a hit to its own entry counts)",
    ),
}
MODEL_MAP = "proteome_map_human_oyster_any_label_uniprot_2026_03"
# No AlphaFold DB model: black, as in the dark-protein map (dark.py).
NO_MODEL = "#000000"


SOURCES = {
    "interpro_match": "InterPro: the matches UniProt lists from InterPro 110.0",
    "sprot_hit": "Swiss-Prot hit: MMseqs2 18.8cc5c against Swiss-Prot 2026_03, identity > 30%",
    "af_model": "AlphaFold DB: model AF-<accession>-F1 (v6)",
}
PX = 10  # pixels per tile in the image, the last one white so single tiles stay visible


def grid(lit):
    """Lit tiles first, then dark, filled row by row; cells past the last protein are NaN."""
    lit = np.sort(lit.astype(float))[::-1]
    rows = -(-len(lit) // COLS)
    g = np.full(rows * COLS, np.nan)
    g[: len(lit)] = lit
    return g.reshape(rows, COLS)


def tiles_image(g):
    """Each tile becomes PX x PX pixels with a one-pixel white gap on its right and bottom."""
    img = np.repeat(np.repeat(g, PX, axis=0), PX, axis=1)
    img[PX - 1 :: PX, :] = np.nan
    img[:, PX - 1 :: PX] = np.nan
    return img


def draw(table, name, cols, what):
    pf.use_style()
    grids, stats = {}, {}
    for sp in SPECIES:
        t = table.filter(pl.col("species") == sp)
        lit = t.select(pl.any_horizontal(cols)).to_series().to_numpy()
        grids[sp] = grid(lit)
        stats[sp] = (int(lit.sum()), len(lit))
    max_rows = max(g.shape[0] for g in grids.values())
    tile_mm = 85 / COLS  # same tile size in both maps
    top_mm = 22 + 2.3 * len(cols)  # finding, tools, legend, then two lines of text per map, then the tiles
    fig_w_mm, fig_h_mm = pf.TWO_COLUMN_MM, top_mm + max_rows * tile_mm + 2
    fig = plt.figure(figsize=(fig_w_mm / 25.4, fig_h_mm / 25.4))
    cmap = ListedColormap([DARK, LIT])
    cmap.set_bad("white")
    # The one fact the figure rests on, first; then where each label comes from.
    dark = {sp: 100 * (n - k) / n for sp, (k, n) in stats.items()}
    what_dark = "dark: none of these"
    finding = f"{dark['oyster']:.1f}% of Pacific oyster proteins are dark (no label), against {dark['human']:.1f}% of human proteins."
    fig.text(2 / fig_w_mm, 1 - 2.5 / fig_h_mm, finding, fontsize=7, weight="bold", va="center", ha="left")
    tools = ["Proteins: UniProt 2026_03 one-protein-per-gene sets. Dark: none of these labels the protein."] + [SOURCES[c] for c in cols]
    fig.text(2 / fig_w_mm, 1 - 4.6 / fig_h_mm, "\n".join(tools), fontsize=6, va="top", ha="left", linespacing=1.25)
    fig.legend(
        handles=[Patch(facecolor=LIT, label=f"lit: the protein has {what}"), Patch(facecolor=DARK, label=what_dark)],
        loc="upper left", bbox_to_anchor=(2 / fig_w_mm, 1 - (7 + 2.3 * len(cols)) / fig_h_mm), frameon=False, ncol=1, borderaxespad=0,
        title=f"One tile = one protein; {COLS} tiles per row and the same tile size in both maps; lit tiles first, then dark",
        title_fontsize=6, alignment="left",
    )
    for i, sp in enumerate(SPECIES):
        g = grids[sp]
        left_mm = 3 + i * 92
        ax = fig.add_axes([
            left_mm / fig_w_mm, 1 - (top_mm + g.shape[0] * tile_mm) / fig_h_mm,
            COLS * tile_mm / fig_w_mm, g.shape[0] * tile_mm / fig_h_mm,
        ])
        ax.imshow(np.ma.masked_invalid(tiles_image(g)), cmap=cmap, vmin=0, vmax=1, interpolation="nearest", aspect="auto")
        ax.set_axis_off()
        k, n = stats[sp]
        fig.text(left_mm / fig_w_mm, 1 - (top_mm - 0.8) / fig_h_mm,
                 f"{SPECIES[sp]}, {n:,} proteins\n{k:,} lit ({100 * k / n:.1f}%), {n - k:,} dark ({100 * (n - k) / n:.1f}%)",
                 fontsize=6, va="bottom", ha="left", linespacing=1.3)
    OUT.mkdir(parents=True, exist_ok=True)
    pf.save(fig, str(OUT / name), formats=("svg", "png", "pdf"))
    return stats


def draw_model(table, name=MODEL_MAP):
    """The first map again, with the dark tiles split by whether AlphaFold DB has a model.
    A model names no function, so a protein with only a model stays dark: the 19.41% of
    oyster proteins the first map shows dark are all still visible here."""
    pf.use_style()
    states, grids = {}, {}
    for sp in SPECIES:
        t = table.filter(pl.col("species") == sp)
        func = (t["interpro_match"] | t["sprot_hit"]).to_numpy()
        model = t["af_model"].to_numpy()
        state = np.where(func, 2, np.where(model, 1, 0))  # 2 function label, 1 model only, 0 neither
        states[sp] = np.bincount(state, minlength=3)
        grids[sp] = grid(state)
    max_rows = max(g.shape[0] for g in grids.values())
    tile_mm = 85 / COLS
    top_mm = 31
    fig_w_mm, fig_h_mm = pf.TWO_COLUMN_MM, top_mm + max_rows * tile_mm + 2
    fig = plt.figure(figsize=(fig_w_mm / 25.4, fig_h_mm / 25.4))
    cmap = ListedColormap([NO_MODEL, DARK, LIT])
    cmap.set_bad("white")
    o = states["oyster"]
    n_o = o.sum()
    fig.text(2 / fig_w_mm, 1 - 2.5 / fig_h_mm,
             f"A model is not a label: {100 * table.filter(pl.col('species') == 'oyster')['af_model'].mean():.2f}% of oyster proteins have an "
             f"AlphaFold DB model; {100 * (o[0] + o[1]) / n_o:.2f}% have no function label.",
             fontsize=7, weight="bold", va="center", ha="left")
    fig.text(2 / fig_w_mm, 1 - 4.6 / fig_h_mm, "\n".join([
        "Proteins: UniProt 2026_03 one-protein-per-gene sets. The same proteins and order of lit and dark as the first map; "
        "here the dark tiles are split by whether AlphaFold DB has a model.",
        SOURCES["interpro_match"], SOURCES["sprot_hit"], SOURCES["af_model"],
    ]), fontsize=6, va="top", ha="left", linespacing=1.25)
    fig.legend(
        handles=[Patch(facecolor=LIT, label="lit: an InterPro match or a Swiss-Prot hit above 30% identity (a hit to its own entry counts)"),
                 Patch(facecolor=DARK, label="dark, but AlphaFold DB has a model: a predicted structure, which names no function"),
                 Patch(facecolor=NO_MODEL, label="dark, and no AlphaFold DB model")],
        loc="upper left", bbox_to_anchor=(2 / fig_w_mm, 1 - 13.5 / fig_h_mm), frameon=False, ncol=1, borderaxespad=0,
        title=f"One tile = one protein; {COLS} tiles per row, the same tile size as the first map; lit tiles first, then dark with a model, then dark without",
        title_fontsize=6, alignment="left",
    )
    for i, sp in enumerate(SPECIES):
        g = grids[sp]
        left_mm = 3 + i * 92
        ax = fig.add_axes([left_mm / fig_w_mm, 1 - (top_mm + g.shape[0] * tile_mm) / fig_h_mm,
                           COLS * tile_mm / fig_w_mm, g.shape[0] * tile_mm / fig_h_mm])
        ax.imshow(np.ma.masked_invalid(tiles_image(g)), cmap=cmap, vmin=-0.5, vmax=2.5, interpolation="nearest", aspect="auto")
        ax.set_axis_off()
        none_, model_only, lit = states[sp]
        n = none_ + model_only + lit
        fig.text(left_mm / fig_w_mm, 1 - (top_mm - 0.8) / fig_h_mm,
                 f"{SPECIES[sp]}, {n:,} proteins\n{lit:,} lit ({100 * lit / n:.1f}%); dark {model_only + none_:,} ({100 * (model_only + none_) / n:.1f}%): "
                 f"{model_only:,} with a model, {none_:,} without",
                 fontsize=6, va="bottom", ha="left", linespacing=1.3)
    OUT.mkdir(parents=True, exist_ok=True)
    pf.save(fig, str(OUT / name), formats=("svg", "png", "pdf"))
    return {sp: v.tolist() for sp, v in states.items()}


def main():
    table = pl.read_parquet("results/why/proteome_labels.parquet")
    for name, (cols, what) in LABEL_SETS.items():
        print(name, draw(table, name, cols, what))
    print(MODEL_MAP, draw_model(table))


if __name__ == "__main__":
    main()
