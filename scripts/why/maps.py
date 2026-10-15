"""Step 5: proteome maps. One tile per protein, the same tile size and 150 tiles per row for
both species; a tile is lit when the protein has any of the labels, dark otherwise.

    python3 -I scripts/why/maps.py
Writes results/why/maps/proteome_map.{svg,png,pdf} and, for the two function labels alone,
results/why/maps/proteome_map_function_labels.{svg,png,pdf}.
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
LIT, DARK = "#E69F00", "#222222"  # Okabe-Ito orange; near-black
SPECIES = {"human": "Human (Homo sapiens)", "oyster": "Pacific oyster (Crassostrea gigas)"}
LABEL_SETS = {
    "proteome_map": (
        ["af_model", "interpro_match", "sprot_hit"],
        "an AlphaFold DB model, an InterPro match, or a Swiss-Prot hit above 30% identity",
    ),
    "proteome_map_function_labels": (
        ["interpro_match", "sprot_hit"],
        "an InterPro match or a Swiss-Prot hit above 30% identity",
    ),
}


def grid(lit):
    """Lit tiles first, then dark, filled row by row; cells past the last protein are NaN."""
    lit = np.sort(lit.astype(float))[::-1]
    rows = -(-len(lit) // COLS)
    g = np.full(rows * COLS, np.nan)
    g[: len(lit)] = lit
    return g.reshape(rows, COLS)


def draw(table, name, cols, what):
    pf.use_style()
    grids, stats = {}, {}
    for sp in SPECIES:
        t = table.filter(pl.col("species") == sp)
        lit = t.select(pl.any_horizontal(cols)).to_series().to_numpy()
        grids[sp] = grid(lit)
        stats[sp] = (int(lit.sum()), len(lit))
    max_rows = max(g.shape[0] for g in grids.values())
    # Same tile size in both panels: each axis is as tall as its own number of rows.
    tile_mm = 80 / COLS
    fig = plt.figure(figsize=(pf.TWO_COLUMN_MM / 25.4, (max_rows * tile_mm + 34) / 25.4))
    fig_w_mm, fig_h_mm = pf.TWO_COLUMN_MM, max_rows * tile_mm + 34
    cmap = ListedColormap([DARK, LIT])
    cmap.set_bad("white")
    fig.legend(
        handles=[Patch(facecolor=LIT, label=f"lit: the protein has {what}"), Patch(facecolor=DARK, label="dark: none of these")],
        loc="upper left", bbox_to_anchor=(6 / fig_w_mm, 1 - 2 / fig_h_mm), frameon=False, ncol=1,
        title=f"One tile = one protein ({COLS} tiles per row, same tile size in both maps)", title_fontsize=6, alignment="left",
    )
    for i, sp in enumerate(SPECIES):
        g = grids[sp]
        left_mm, top_mm = 6 + i * 91, 22
        ax = fig.add_axes([
            left_mm / fig_w_mm, 1 - (top_mm + g.shape[0] * tile_mm) / fig_h_mm,
            COLS * tile_mm / fig_w_mm, g.shape[0] * tile_mm / fig_h_mm,
        ])
        ax.imshow(np.ma.masked_invalid(g), cmap=cmap, vmin=0, vmax=1, interpolation="nearest", aspect="equal")
        ax.set_axis_off()
        k, n = stats[sp]
        fig.text(left_mm / fig_w_mm, 1 - (top_mm - 1.2) / fig_h_mm,
                 f"{SPECIES[sp]}: {k:,} of {n:,} proteins lit ({100 * k / n:.1f}%), {n - k:,} dark ({100 * (n - k) / n:.1f}%)",
                 fontsize=6, va="bottom", ha="left")
    OUT.mkdir(parents=True, exist_ok=True)
    pf.save(fig, str(OUT / name), formats=("svg", "png", "pdf"))
    return stats


def main():
    table = pl.read_parquet("results/why/proteome_labels.parquet")
    for name, (cols, what) in LABEL_SETS.items():
        print(name, draw(table, name, cols, what))


if __name__ == "__main__":
    main()
