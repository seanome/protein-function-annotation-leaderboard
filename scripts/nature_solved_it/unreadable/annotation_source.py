"""Where each proteome's gene models come from, and the group medians with and without the
five B10K control birds.

The five random-control birds at 300-400 million years have 0.02-0.04% unlabelled proteins
and no protein named "Uncharacterized protein". This script reads, for each of the 42
proteomes of the 2026-10-20 run, UniProt's proteome record (annotation source, assembly,
citation) and, for exactly the proteins the run counted (the accessions in the run's FASTA),
whether UniProt gives a gene symbol and whether it flags the protein as a fragment. It then
recomputes the group medians of PLAN.md's rule with and without the B10K proteomes. The rule
and the verdict are not changed; this is a sensitivity table only.

Inputs, read 2026-10-20 (not committed, ~150 MB):
    <dl>/<upid>.json          https://rest.uniprot.org/proteomes/<upid>?format=json
    <dl>/genes/<upid>.tsv.gz  https://rest.uniprot.org/uniprotkb/stream?query=proteome:<upid>
                              &fields=accession,gene_primary,gene_oln,gene_orf,fragment,
                              protein_existence,protein_name&format=tsv&compressed=true
    <dl>/fasta/<file>         the run's FASTA (URL and sha256 from run_2026-10-20/inputs/
                              <upid>_manifest.json; all 42 matched the run's sha256)

    /Users/olga/anaconda3/envs/2025-kmerseek-analysis/bin/python3 \
        scripts/nature_solved_it/unreadable/annotation_source.py <dl>
"""

import gzip
import json
import pathlib
import sys

import polars as pl

REPO = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts/why"))
import pubfig as pf  # noqa: E402

RUN = REPO / "results/nature-solved-it/unreadable/run_2026-10-20"
OUT = REPO / "results/nature-solved-it/unreadable/annotation_source"
GROUPS = ["model organism", "atlas species", "random control"]
REFERENCE = "leaderboard reference"
BUSCO_MIN = 90.0
B10K_TITLE = "Bird 10,000 Genomes (B10K) Project - Family phase."


def fasta_accessions(path):
    with gzip.open(path, "rt") as fh:
        return {line.split("|")[1] for line in fh if line.startswith(">")}


def proteome_row(dl, upid):
    meta = json.loads((dl / f"{upid}.json").read_text())
    titles = [c.get("title") or "; ".join(c.get("authors") or []) for c in meta.get("citations", [])]
    fasta = json.loads((RUN / "inputs" / f"{upid}_manifest.json").read_text())["fasta"]["file"]
    accs = fasta_accessions(dl / "fasta" / fasta)
    genes = pl.read_csv(gzip.decompress((dl / "genes" / f"{upid}.tsv.gz").read_bytes()), separator="\t",
                        infer_schema_length=0, quote_char=None)
    genes = genes.filter(pl.col("Entry").is_in(accs))
    assert genes.height == len(accs), (upid, genes.height, len(accs))
    n = genes.height
    return {
        "upid": upid,
        "annotation_source": (meta.get("genomeAnnotation") or {}).get("source"),
        "assembly": (meta.get("genomeAssembly") or {}).get("assemblyId"),
        "first_citation": titles[0] if titles else "",
        "b10k": B10K_TITLE in titles,
        "gene_symbol_pct": round(100 * genes["Gene Names (primary)"].is_not_null().sum() / n, 2),
        "fragment_pct": round(100 * (genes["Fragment"] == "fragment").sum() / n, 2),
        "pe_homology_pct": round(100 * (genes["Protein existence"] == "Inferred from homology").sum() / n, 2),
    }


def medians(df):
    return {
        g: {"n": df.filter(pl.col("group") == g).height,
            "median_pct": df.filter(pl.col("group") == g)["unlabelled_pct"].median()}
        for g in GROUPS
    }


def main(dl):
    table = pl.read_csv(RUN / "unreadable_table.tsv", separator="\t")
    src = pl.DataFrame([proteome_row(dl, u) for u in table["upid"]])
    t = table.join(src, on="upid", how="left").with_columns(
        sprot_hit_pct=(100 * pl.col("n_sprot_hit") / pl.col("n_proteins")).round(2)
    )
    OUT.mkdir(exist_ok=True)
    cols = ["group", "species", "upid", "divergence_my", "n_proteins", "unlabelled_pct",
            "uncharacterized_name_pct", "sprot_hit_pct", "busco_completeness_pct",
            "annotation_source", "first_citation", "b10k", "assembly", "gene_symbol_pct",
            "fragment_pct", "pe_homology_pct"]
    t.select(cols).write_csv(OUT / "annotation_source_uniprot_2026_03.tsv", separator="\t")

    groups = t.filter(pl.col("group") != REFERENCE)
    no_b10k = groups.filter(~pl.col("b10k"))
    good = groups.filter(pl.col("busco_completeness_pct") >= BUSCO_MIN)
    sens = {
        "read_utc_date": "2026-10-20",
        "b10k_proteomes": groups.filter(pl.col("b10k")).select("group", "species", "upid").to_dicts(),
        "b10k_with_busco_90": good.filter(pl.col("b10k")).height,
        "medians_all": medians(groups),
        "medians_all_without_b10k": medians(no_b10k),
        "medians_busco_90": medians(good),
        "medians_busco_90_without_b10k": medians(good.filter(~pl.col("b10k"))),
    }
    for k in ("medians_all", "medians_all_without_b10k"):
        m = sens[k]
        sens[k + "_atlas_minus_controls"] = round(m["atlas species"]["median_pct"] - m["random control"]["median_pct"], 2)
    (OUT / "sensitivity_without_b10k.json").write_text(json.dumps(sens, indent=1) + "\n")
    print(json.dumps(sens, indent=1))
    with pl.Config(tbl_rows=50, tbl_cols=20, fmt_str_lengths=40):
        print(t.select(cols[:3] + ["unlabelled_pct", "uncharacterized_name_pct", "sprot_hit_pct",
                                   "busco_completeness_pct", "annotation_source", "b10k",
                                   "gene_symbol_pct", "fragment_pct", "pe_homology_pct"]))
    figure(groups)


def figure(groups):
    pf.use_style()
    colour = dict(zip(GROUPS, [pf.OKABE_ITO["black"], pf.OKABE_ITO["vermillion"], pf.OKABE_ITO["blue"]]))
    fig, axs = pf.figure(pf.TWO_COLUMN_MM, 62, ncols=2, gridspec_kw={"width_ratios": [1, 1.25]})

    # a: unlabelled share per group; B10K proteomes as open squares; medians with and without them
    ax = axs[0]
    for i, g in enumerate(GROUPS):
        d = groups.filter(pl.col("group") == g)
        for is_b10k, marker in ((False, "o"), (True, "s")):
            v = d.filter(pl.col("b10k") == is_b10k)["unlabelled_pct"].to_list()
            if v:
                ax.scatter([i + (0.12 if is_b10k else -0.08)] * len(v), v, s=10, marker=marker, color=colour[g],
                           facecolor="none" if is_b10k else colour[g], linewidths=0.6, zorder=3)
        ax.hlines(d["unlabelled_pct"].median(), i - 0.32, i + 0.32, color="black", lw=1, zorder=4)
        nb = d.filter(~pl.col("b10k"))
        if nb.height < d.height:
            ax.hlines(nb["unlabelled_pct"].median(), i - 0.32, i + 0.32, color="black", lw=1,
                      linestyles=(0, (2, 1.5)), zorder=4)
    ax.annotate("5 B10K birds,\n0.02-0.04%", xy=(2.14, 0.6), xytext=(2.5, 18), ha="center",
                arrowprops={"arrowstyle": "-", "lw": 0.5})
    ax.set_xticks(range(3), [f"{g.split()[0].capitalize()}\n{g.split()[1]}\nn = {groups.filter(pl.col('group') == g).height}" for g in GROUPS])
    ax.set_xlim(-0.6, 2.9)
    ax.set_ylabel("Unlabelled share (% of proteins)")

    # b: share flagged as fragment against unlabelled share, every proteome
    ax = axs[1]
    for g in GROUPS:
        for is_b10k, marker in ((False, "o"), (True, "s")):
            d = groups.filter((pl.col("group") == g) & (pl.col("b10k") == is_b10k))
            if d.height:
                ax.scatter(d["fragment_pct"], d["unlabelled_pct"], s=10, marker=marker, color=colour[g],
                           facecolor="none" if is_b10k else colour[g], linewidths=0.6, zorder=3)
    ax.annotate("5 B10K birds:\nevery protein a fragment", xy=(99, 0.6), xytext=(78, 17), ha="center",
                arrowprops={"arrowstyle": "-", "lw": 0.5})
    ax.set_xlabel("Proteins flagged as fragment by UniProt (%)")
    ax.set_ylabel("Unlabelled share (% of proteins)")
    ax.set_xlim(-3, 103)

    # one legend for both panels: colour = group, shape = gene models from B10K or not, line = median
    from matplotlib.lines import Line2D
    blank = Line2D([], [], ls="", label=" ")
    h = [Line2D([], [], ls="", marker="o", color=colour[g], ms=3.5, label=g.capitalize()) for g in GROUPS]
    h += [Line2D([], [], ls="", marker="o", color="0.3", ms=3.5, label="Gene models from any other project"),
          Line2D([], [], ls="", marker="s", color="0.3", mfc="none", ms=3.5, label="B10K gene models (copied from reference genes)"),
          blank,
          Line2D([], [], color="black", lw=1, label="Group median, all proteomes"),
          Line2D([], [], color="black", lw=1, ls=(0, (2, 1.5)), label="Group median, B10K left out"),
          blank]
    fig.legend(handles=h, loc="outside upper center", ncol=3, frameon=False)
    for ax, letter in zip(axs, "ab"):
        pf.panel_label(ax, letter)
    pf.save(fig, str(OUT / "unlabelled_share_b10k_gene_models_uniprot_2026_03"))


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))
