"""The measured numbers the organism pages read, as one JSON file (docs/data/organism/strata.json).

Everything here comes from data/results.parquet, data/compute.parquet and data/pileup/, the
same tables the leaderboard reads, through build_site.load(). No number is made up: a cell
the tables do not have is left out, and the page shows a dash and "not run" for it.

Strata. A species is placed by two numbers from its cards:
  distance  the TimeTree age to its nearest reviewed relative, matched to the measured
            relative whose age (to human) is closest. The leaderboard searched human proteins
            against 15 relatives, 6 to 1,516 million years (Mya) away by TimeTree 5's
            pairwise ages; E. coli has no TimeTree age and is never matched.
  pLDDT     the mean AlphaFold pLDDT of its proteome, matched to the pLDDT band the tables
            use for the domains (0-50, 50-70, 70-90, 90-100).
"""

import math
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_site as bs  # noqa: E402
import pileup as pu  # noqa: E402
from tools import SPECIES, TOOLS  # noqa: E402

MIN_N = 200
PLDDT_BANDS = ["0-50", "50-70", "70-90", "90-100"]
# NCBI taxon id of each measured relative (the species of the QfO 2020_04 proteome).
RELATIVE_TAXID = {
    "ptroglodytes": 9598, "rnorvegicus": 10116, "clupus": 9615, "btaurus": 9913, "mouse": 10090,
    "mdomestica": 13616, "chicken": 9031, "xtropicalis": 8364, "zebrafish": 7955, "ciona": 7719,
    "fly": 7227, "worm": 6239, "yeast": 4932, "arabidopsis": 3702, "ecoli": 562,
}
SCIENTIFIC = {
    "ptroglodytes": "Pan troglodytes", "rnorvegicus": "Rattus norvegicus", "clupus": "Canis lupus",
    "btaurus": "Bos taurus", "mouse": "Mus musculus", "mdomestica": "Monodelphis domestica",
    "chicken": "Gallus gallus", "xtropicalis": "Xenopus tropicalis", "zebrafish": "Danio rerio",
    "ciona": "Ciona intestinalis", "fly": "Drosophila melanogaster", "worm": "Caenorhabditis elegans",
    "yeast": "Saccharomyces cerevisiae", "arabidopsis": "Arabidopsis thaliana", "ecoli": "Escherichia coli",
}
# The organisms of the dark-protein measurement (dark_set_gain), NCBI species taxon ids.
DARK_TAXID = {"human": 9606, "yeast": 4932, "worm": 6239, "borrelia": 139, "botryllus": 30301}
DARK_NAME = {"human": "human (Homo sapiens)", "yeast": "baker's yeast (Saccharomyces cerevisiae)",
             "worm": "nematode (Caenorhabditis elegans)", "borrelia": "Lyme disease bacterium (Borrelia burgdorferi)",
             "botryllus": "star ascidian (Botryllus schlosseri)"}
PILEUP_SPECIES = "zebrafish"


def f(x):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else float(x)


def band_of(v):
    if v is None:
        return None
    for b in PLDDT_BANDS:
        lo, hi = (float(x) for x in b.split("-"))
        if lo <= v < hi or (hi == 100 and v == 100):
            return b
    return None


def build(data: Path) -> dict:
    res, comp, src, prov, cov = bs.load(data)
    long = pl.read_parquet(data / "results.parquet")
    pf = res.filter((pl.col("truth_set") == "pfam") & (pl.col("calls") == "raw") & (pl.col("tool") != "hmmscan"))
    rel = (pf.filter(pl.col("species") != "all").select("species", "species_mya").unique()
           .sort("species_mya", "species"))
    relatives = [{"id": r["species"], "name": SPECIES[r["species"]], "scientific": SCIENTIFIC[r["species"]],
                  "taxid": RELATIVE_TAXID[r["species"]], "mya": f(r["species_mya"])} for r in rel.iter_rows(named=True)]

    # Best setting of each tool at each relative: highest Fmax over the whole set, ties to the
    # setting name (the rule the tool cards use).
    best = (pf.filter(pl.col("stratum_axis") == "all").sort(["fmax", "setting"], descending=[True, False])
            .group_by("tool", "species", maintain_order=True).first())
    best_setting = {(r["tool"], r["species"]): r["setting"] for r in best.iter_rows(named=True)}

    def cell(q):
        return {"fmax": f(q["fmax"]), "recall": f(q["fmax_recall"]), "precision": f(q["fmax_precision"]),
                "iou": f(q["median_iou_tp"]), "n": int(q["n_truth_instances"])}

    cells = {}
    for r in pf.filter(pl.col("stratum_axis").is_in(["all", "plddt"])).iter_rows(named=True):
        if r["species"] == "all" or best_setting.get((r["tool"], r["species"])) != r["setting"]:
            continue
        band = "all" if r["stratum_axis"] == "all" else r["stratum"]
        c = cells.setdefault(f'{r["species"]}|{band}', {"n": int(r["n_truth_instances"]), "tools": {}})
        if c["n"] != int(r["n_truth_instances"]):
            raise SystemExit(f"two instance counts in {r['species']}|{band}")
        c["tools"][r["tool"]] = dict(cell(r), setting=r["setting"])

    # Distance curve: each tool's Fmax at every relative, for every setting that is some
    # relative's best (so the page can hold the setting fixed along the curve).
    curve = {}
    keep = {(t, s) for (t, _), s in best_setting.items()}
    for r in pf.filter(pl.col("stratum_axis") == "all").iter_rows(named=True):
        if (r["tool"], r["setting"]) in keep and r["species"] != "all":
            curve.setdefault(f'{r["tool"]}|{r["setting"]}', {})[r["species"]] = f(r["fmax"])

    # Short features: recall at the Fmax cut-off on Swiss-Prot features 2-15 residues long.
    sp = res.filter((pl.col("truth_set") == "swissprot") & (pl.col("calls") == "raw")
                    & (pl.col("stratum_axis") == "feature_length_bin") & (pl.col("stratum") == "2-15"))
    short = {}
    for r in sp.iter_rows(named=True):
        if (r["tool"], r["setting"]) in keep:
            short[f'{r["tool"]}|{r["setting"]}|{r["species"]}'] = {"recall": f(r["fmax_recall"]), "precision": f(r["fmax_precision"]),
                                                                      "n": int(r["n_truth_instances"])}

    # Cost per 1,000 query proteins, search step per relative, as the leaderboard computes it.
    cost = {}
    for r in comp.iter_rows(named=True):
        k = f'{r["tool"]}|{r["setting"]}'
        if k in cost:
            continue
        gpu = r["gpu_h_per_1000_queries"] or 0.0
        cost[k] = {"cpu_h_per_1000": f(r["cpu_h_per_1000_queries"]), "gpu_h_per_1000": f(gpu),
                   "cost_usd_est": f(r["cpu_h_per_1000_queries"] * bs.R_CPU + gpu * bs.R_GPU),
                   "db_cpu_h_per_db": None if r["db_cpu_h"] is None else f(r["db_cpu_h"] / r["db_n"]),
                   "n_relatives": r["n_relatives"]}

    # Pile-up (zebrafish only): per tool and setting, per pLDDT band of the domain, the share
    # of domains covered at each of 300 columns; the TP domains (for the overlap figure); and
    # the median distance from a correct call's start to the true start.
    pileup = None
    if (data / "pileup/regions.parquet").exists():
        reg = pl.read_parquet(data / "pileup/regions.parquet").sort("row_id")
        calls = pl.read_parquet(data / "pileup/calls.parquet")
        panels = [{"tool": t, "setting": s, "kind": TOOLS[t]["kind"]}
                  for t, s in calls.select("tool", "setting").unique().sort("tool", "setting").iter_rows()]
        pdata, _ = pu.build(data, panels)
        regions = pdata["regions"]
        # Band from the unrounded pLDDT (pileup.build rounds it to 0.1, which moves domains at
        # 89.95-89.99 into 90-100).
        raw_plddt = dict(zip(reg["row_id"].to_list(), reg["region_plddt"].to_list()))
        band = [band_of(raw_plddt[r["id"]]) for r in regions]
        tp = calls.filter(pl.col("is_tp")).group_by("tool", "setting").agg(pl.col("row_id").unique().sort())
        tp_rows = {f"{t}|{s}": ids for t, s, ids in tp.iter_rows()}
        start = (calls.filter(pl.col("is_tp")).join(reg.select("row_id", "true_start", "region_plddt"), on="row_id")
                 .with_columns((pl.col("call_start") - pl.col("true_start")).abs().alias("d")))
        panels_out = {}
        for p in pdata["panels"]:
            k = f'{p["tool"]}|{p["setting"]}'
            by_band = {}
            for b in PLDDT_BANDS + ["all"]:
                idx = [i for i in range(len(regions)) if b == "all" or band[i] == b]
                if not idx:
                    continue
                m = np.zeros((len(idx), pu.NCOL), dtype=float)
                ok = np.ones((len(idx), pu.NCOL), dtype=bool)
                for j, i in enumerate(idx):
                    r = regions[i]
                    if r["c0"] is None:
                        ok[j, :] = False
                    else:
                        ok[j, :r["c0"]] = False
                        ok[j, r["c1"] + 1:] = False
                    for a, bb, _v in p["rows"][i]["r"]:
                        m[j, a:bb + 1] = 1
                with np.errstate(invalid="ignore", divide="ignore"):
                    curve_b = (m * ok).sum(0) / ok.sum(0)
                ids = {regions[i]["id"] for i in idx}
                d = start.filter((pl.col("tool") == p["tool"]) & (pl.col("setting") == p["setting"])
                                 & pl.col("row_id").is_in(list(ids)))
                by_band[b] = {"n": len(idx),
                              "curve": [None if np.isnan(x) else round(float(x), 4) for x in curve_b],
                              "median_start_offset_aa": f(d["d"].median()) if d.height else None,
                              "n_tp_calls": d.height}
            panels_out[k] = {"bands": by_band, "tp_rows": tp_rows.get(k, [])}
        pileup = {"species": PILEUP_SPECIES, "flank": pu.FLANK, "body": pu.BODY,
                  "median_domain_length": float(np.median([r["len"] for r in regions])),
                  "region_band": {str(r["id"]): b for r, b in zip(regions, band)},
                  "panels": panels_out}

    # Dark proteins: no sequence-tool hit at E < 1 against Swiss-Prot with the organism's own
    # clade removed (dark_set_gain); share = dark / (dark + placed).
    dk = long.filter(pl.col("experiment") == "dark_set_gain")
    dark = []
    for org in DARK_TAXID:
        d = dk.filter(pl.col("query_set").str.starts_with(f"{org} dark") & (pl.col("metric") == "dark_reached"))
        p = dk.filter(pl.col("query_set").str.starts_with(f"{org} placed") & (pl.col("metric") == "placed_reached"))
        if d.height == 0 or p.height == 0:
            continue
        n_dark, n_placed = int(d["n_instances"].max()), int(p["n_instances"].max())
        dark.append({"organism": org, "name": DARK_NAME[org], "taxid": DARK_TAXID[org], "dark": n_dark,
                     "proteins": n_dark + n_placed, "dark_fraction": n_dark / (n_dark + n_placed),
                     "no_relative_fraction": n_dark / (n_dark + n_placed),
                     "disordered_fraction": None, "short_motif_fraction": None,
                     "query_set": d["query_set"][0], "target_db": d["target_db"][0],
                     "source_path": d["source_path"][0], "source_commit": d["source_commit"][0]})

    return {
        "source": src, "min_n": MIN_N, "plddt_bands": PLDDT_BANDS, "relatives": relatives,
        "truth": "Pfam-A domain instances (Pfam 38.2), raw calls, experiment " + bs.EXPERIMENT,
        "n_queries": prov["n_queries"], "cells": cells, "curve": curve, "short": short, "cost": cost,
        "prices": {"r_cpu": bs.R_CPU, "r_gpu": bs.R_GPU, "page": bs.PRICE_PAGE},
        "pileup": pileup, "dark": dark,
        "tools": {t: {"name": TOOLS[t]["name"], "kind": TOOLS[t]["kind"], "shape": TOOLS[t]["shape"],
                      "color_var": f"--k-{TOOLS[t]['kind']}", "command": TOOLS[t]["command"],
                      "container": TOOLS[t]["container"]}
                  for t in sorted(set(pf["tool"].unique()))},
        "setting_label": {f"{t}|{s}": bs.setting_label(t, s) for (t, s) in sorted(keep)},
    }
