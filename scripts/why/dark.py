"""Step 5: what the dark proteins are. Dark against lit, per species, from the per-protein
table, and a map of the dark proteins alone coloured by AlphaFold DB pLDDT band.

Dark is the page's rule: no InterPro match and no Swiss-Prot hit above 30% identity (a
hit to the protein's own entry counts). The 19.41% of oyster proteins in step 4.

    python3 -I scripts/why/dark.py
Reads results/why/proteome_labels.parquet (and its seg_fraction column when labels.py has
been rerun with --seg). Writes the step 5 rows of results/why/numbers.parquet and
results/why/maps/proteome_map_dark_by_plddt.{svg,png,pdf}.
"""

import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pubfig as pf  # noqa: E402
from four import BANDS  # noqa: E402
from labels import SEG_HIGH, length_bin  # noqa: E402
from maps import COLS, OUT, PX, SPECIES, tiles_image  # noqa: E402
from numbers_table import Numbers  # noqa: E402

STEP = 5
# No AlphaFold DB model for the accession. Black: no other tile or bar on the page is
# black, and it is never folded into "pLDDT 50 or less".
NO_MODEL = "#000000"
REL = "UniProt 2026_03"
TICK = "#555555"
DARK_RULE = "dark: no InterPro match and no Swiss-Prot hit above 30% identity (a hit to its own entry counts)"


def with_dark(table):
    t = table.with_columns(dark=~(pl.col("interpro_match") | pl.col("sprot_hit")))
    if "length_bin" not in t.columns:
        t = t.with_columns(length_bin=length_bin(pl.col("length")))
    return t


def band_index(plddt):
    """0-3 for the four AlphaFold DB bands in four.py's order, 4 for no model."""
    out = np.full(len(plddt), len(BANDS), dtype=float)
    for i, (lo, hi, _, _) in enumerate(BANDS):
        out[(plddt > lo) & (plddt <= hi)] = i
    return out


def leading_no_model(grid):
    """How many of the longest proteins in a row have no model before the first one that has one."""
    flat = grid.ravel()
    flat = flat[~np.isnan(flat)]
    has = np.flatnonzero(flat != len(BANDS))
    return int(has[0]) if has.size else len(flat)


def numbers(t):
    N = Numbers()
    has_seg = "seg_fraction" in t.columns
    for sp in SPECIES:
        for dark, grp in ((True, "dark"), (False, "lit")):
            g = t.filter((pl.col("species") == sp) & (pl.col("dark") == dark))
            n = g.height
            who = f"{sp}, {grp} proteins"
            pid = f"{sp}_{grp}"
            N.add(f"{pid}_n", STEP, f"{who}: count", n, "proteins", "sherlock_mmseqs", REL, note=DARK_RULE)
            N.add(f"{pid}_length_median", STEP, f"{who}: median length", g["length"].median(), "aa", f"{sp}_canonical_fasta", REL)
            N.add(f"{pid}_under_100aa_pct", STEP, f"{who}: share shorter than 100 aa", round(100 * (g["length"] < 100).sum() / n, 2), "% of proteins", f"{sp}_canonical_fasta", REL)
            low = int((g["af_mean_plddt"] <= 50).sum())
            N.add(f"{pid}_plddt_50_or_less_pct", STEP, f"{who}: share whose AlphaFold DB model has mean pLDDT 50 or less", round(100 * low / n, 2), "% of proteins", "afdb_predictions", "AlphaFold DB API",
                  note="denominator is every protein in the group; proteins with no model are not in the numerator")
            N.add(f"{pid}_no_model_pct", STEP, f"{who}: share with no AlphaFold DB model for the accession (AF-accession-F1)", round(100 * (~g["af_model"]).sum() / n, 2), "% of proteins", "afdb_predictions", "AlphaFold DB API")
            frag = int((g["busco_status"] == "Fragmented").sum())
            N.add(f"{pid}_busco_fragmented_pct", STEP, f"{who}: share BUSCO metazoa_odb12 calls Fragmented", round(100 * frag / n, 2), "% of proteins", f"sherlock_busco_{sp}", "BUSCO 6.1.0, metazoa_odb12")
            if has_seg:
                k = int((g["seg_fraction"] > SEG_HIGH).sum())
                N.add(f"{pid}_seg_over_20pct_pct", STEP, f"{who}: share with more than 20% of residues masked by SEG", round(100 * k / n, 2), "% of proteins", f"{sp}_segmasker", "BLAST+ 2.16.0 segmasker, defaults")
    N.write("results/why/numbers.parquet", replace_steps={STEP})
    return N.frame()


def draw(t):
    """Dark proteins only, one tile each, longest first, coloured by pLDDT band."""
    pf.use_style()
    colours = [c for _, _, c, _ in BANDS] + [NO_MODEL]
    cmap = ListedColormap(colours)
    cmap.set_bad("white")
    tile_mm = 85 / COLS  # the same tile size as the step 4 maps
    panels = {}
    for sp in ("oyster", "human"):
        g = t.filter((pl.col("species") == sp) & pl.col("dark")).sort(["length", "accession"], descending=[True, False])
        idx = band_index(g["af_mean_plddt"].fill_null(np.nan).to_numpy())
        rows = -(-len(idx) // COLS)
        grid = np.full(rows * COLS, np.nan)
        grid[: len(idx)] = idx
        panels[sp] = (grid.reshape(rows, COLS), g["length"].to_numpy(), np.bincount(idx.astype(int), minlength=len(colours)))

    o = t.filter((pl.col("species") == "oyster"))
    low_dark = 100 * (o.filter(pl.col("dark"))["af_mean_plddt"] <= 50).sum() / o.filter(pl.col("dark")).height
    low_lit = 100 * (o.filter(~pl.col("dark"))["af_mean_plddt"] <= 50).sum() / o.filter(~pl.col("dark")).height

    gap_mm, label_mm = 6, 32
    top_mm = 19
    heights = [panels[sp][0].shape[0] * tile_mm for sp in panels]
    fig_w_mm = pf.TWO_COLUMN_MM
    fig_h_mm = top_mm + sum(heights) + gap_mm * len(panels) + 2
    fig = plt.figure(figsize=(fig_w_mm / 25.4, fig_h_mm / 25.4))
    fig.text(2 / fig_w_mm, 1 - 2.5 / fig_h_mm,
             f"{low_dark:.1f}% of dark Pacific oyster proteins have an AlphaFold DB model with mean pLDDT 50 or less, "
             f"against {low_lit:.1f}% of lit ones (proteins with no model count in both totals).",
             fontsize=7, weight="bold", va="center", ha="left")
    fig.text(2 / fig_w_mm, 1 - 4.6 / fig_h_mm, "\n".join([
        "Proteins: UniProt 2026_03 one-protein-per-gene sets; only the dark ones are drawn. "
        "Dark: no InterPro match and no Swiss-Prot hit above 30% identity (MMseqs2 18.8cc5c).",
        "Colour: mean pLDDT of the AlphaFold DB model AF-<accession>-F1 (v6), in AlphaFold DB's colours; bands as in the four-protein "
        "figures: above 90, above 70 up to 90, above 50 up to 70, 50 or less.",
    ]), fontsize=6, va="top", ha="left", linespacing=1.25)
    fig.legend(
        handles=[Patch(facecolor=c, label=lab) for (_, _, c, lab) in BANDS] + [Patch(facecolor=NO_MODEL, label="no AlphaFold DB model"),
                 Line2D([], [], marker="<", markersize=3, color=TICK, linestyle="none", label="row holding every 1,000th tile")],
        loc="upper left", bbox_to_anchor=(2 / fig_w_mm, 1 - 9.5 / fig_h_mm), frameon=False, ncol=6, borderaxespad=0,
        title=f"One tile = one dark protein; {COLS} tiles per row, the same tile size as the two proteome maps on this page; "
              "longest protein first, read row by row.",
        title_fontsize=6, alignment="left",
    )
    y_mm = top_mm
    for sp, (grid, lengths, counts) in panels.items():
        n = len(lengths)
        fig.text(3 / fig_w_mm, 1 - (y_mm - 0.8) / fig_h_mm,
                 f"{SPECIES[sp]}: {n:,} dark proteins, {lengths[0]:,} to {lengths[-1]:,} aa. "
                 + ", ".join(f"{lab} {c:,}" for (_, _, _, lab), c in zip(BANDS, counts)) + f", no model {counts[-1]:,}"
                 + (f" (the {lead:,} longest, {lengths[lead - 1]:,} aa and up, all have no model)" if (lead := leading_no_model(grid)) > 1 else ""),
                 fontsize=6, va="bottom", ha="left")
        h_mm = grid.shape[0] * tile_mm
        ax = fig.add_axes([3 / fig_w_mm, 1 - (y_mm + h_mm) / fig_h_mm, COLS * tile_mm / fig_w_mm, h_mm / fig_h_mm])
        ax.imshow(np.ma.masked_invalid(tiles_image(grid)), cmap=cmap, vmin=-0.5, vmax=len(colours) - 0.5,
                  interpolation="nearest", aspect="auto")
        ax.set_xlim(-0.5, COLS * PX - 0.5)
        ax.set_ylim(grid.shape[0] * PX - 0.5, -0.5)
        ax.set_axis_off()
        for k in range(1000, n + 1, 1000):
            r, c = divmod(k - 1, COLS)
            # A triangle in the right margin points at the row holding the tile; inside the
            # grid a tick would be a sixth mark on top of five colours.
            ax.plot([COLS * PX + 6], [r * PX + PX / 2], marker="<", markersize=3, color=TICK, clip_on=False)
            fig.text((3 + COLS * tile_mm + 2.5) / fig_w_mm, 1 - (y_mm + (r + 0.5) * tile_mm) / fig_h_mm,
                     f"tile {k:,} (column {c + 1} of {COLS}): {lengths[k - 1]:,} aa", fontsize=5, va="center", ha="left")
        y_mm += h_mm + gap_mm
    assert 3 + COLS * tile_mm + 1.5 + label_mm < fig_w_mm
    OUT.mkdir(parents=True, exist_ok=True)
    pf.save(fig, str(OUT / "proteome_map_dark_by_plddt"), formats=("svg", "png", "pdf"))
    return {sp: dict(zip([lab for *_, lab in BANDS] + ["no model"], p[2].tolist())) for sp, p in panels.items()}


def main():
    t = with_dark(pl.read_parquet("results/why/proteome_labels.parquet"))
    print(numbers(t).select("id", "value", "unit").to_pandas().to_string())
    print(draw(t))


if __name__ == "__main__":
    main()
