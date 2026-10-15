#!/usr/bin/env python3
"""Build the region pile-up data, the way deepTools computeMatrix scale-regions does, on residues.

One row per true region. 300 columns:
    0-99     the 100 residues before the region's start (true_start - 100 ... true_start - 1)
    100-199  the region itself, scaled to 100 bins; bin k takes the residue at its centre,
             true_start + floor((k + 0.5) * L / 100), L = region length
    200-299  the 100 residues after its end (true_end + 1 ... true_end + 100)
A column whose residue lies before residue 1 or past the protein's last residue is null.
Every other cell is the highest score of the panel's calls (same protein, same Pfam family)
covering that residue, or nothing when no call covers it. Coverage is "a score is there".

The page receives each row's runs of equal value, so it rebuilds exactly this matrix; the
heatmap, the column-average curve, the readouts and the caption all come from that rebuild.

Written beside it, for the check: per-residue recall computed from the raw call intervals
without the matrix (fraction of each region's residues covered, averaged over regions).
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import polars as pl

FLANK, BODY = 100, 100
NCOL = FLANK + BODY + FLANK


def column_residues(s: int, e: int) -> np.ndarray:
    L = e - s + 1
    before = s - FLANK + np.arange(FLANK)
    body = s + np.floor((np.arange(BODY) + 0.5) * L / BODY).astype(int)
    after = e + 1 + np.arange(FLANK)
    return np.concatenate([before, body, after])


def runs(values: np.ndarray, covered: np.ndarray) -> list:
    """[[first_col, last_col, value], ...] for consecutive covered columns of equal value."""
    out, j = [], 0
    while j < NCOL:
        if not covered[j]:
            j += 1
            continue
        k = j
        while k + 1 < NCOL and covered[k + 1] and values[k + 1] == values[j]:
            k += 1
        v = float(values[j])
        out.append([j, k, float(f"{v:.4g}")])
        j = k + 1
    return out


def merge(intervals: list[tuple[int, int]]) -> list[list[int]]:
    out: list[list[int]] = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def build(data: Path, panels: list[dict]) -> tuple[dict, dict]:
    reg = pl.read_parquet(data / "pileup/regions.parquet").sort("row_id")
    calls = pl.read_parquet(data / "pileup/calls.parquet")
    rows = reg.to_dicts()
    cols = [column_residues(r["true_start"], r["true_end"]) for r in rows]
    null = [(c < 1) | (c > r["protein_length"]) for c, r in zip(cols, rows)]

    out_panels, check = [], {}
    for p in panels:
        sub = calls.filter((pl.col("tool") == p["tool"]) & (pl.col("setting") == p["setting"]))
        by_row: dict[int, list] = {}
        for r in sub.iter_rows(named=True):
            by_row.setdefault(r["row_id"], []).append((r["call_start"], r["call_end"], r["score"]))
        prow, recall = [], []
        for i, r in enumerate(rows):
            s, e = r["true_start"], r["true_end"]
            cs = by_row.get(r["row_id"], [])
            lo, hi = s - FLANK, e + FLANK
            best = np.full(hi - lo + 1, -np.inf)
            for a, b, sc in cs:
                a2, b2 = max(a, lo), min(b, hi)
                if a2 <= b2:
                    seg = best[a2 - lo:b2 - lo + 1]
                    np.maximum(seg, sc, out=seg)
            res = cols[i]
            val = best[res - lo]
            cov = np.isfinite(val) & ~null[i]
            # Per-residue recall straight from the intervals, no bins: covered residues of
            # the region over its length.
            # Only residues inside the protein count; a region past the protein's end has
            # those residues null in the matrix too.
            inside = np.arange(s, e + 1) <= r["protein_length"]
            if inside.any():
                recall.append(float(np.isfinite(best[s - lo:e - lo + 1])[inside].mean()))
            # The call with the largest overlap with the true region gives the start offset.
            best_ov, d = 0, None
            for a, b, _sc in cs:
                ov = min(b, e) - max(a, s) + 1
                if ov > best_ov:
                    best_ov, d = ov, a - s
            prow.append({
                "r": runs(val, cov),
                "c": merge([(max(a, 1), min(b, r["protein_length"])) for a, b, _ in cs
                            if b >= s and a <= e]),
                "n": len(cs),
                "d": d,
            })
        out_panels.append(dict(p, rows=prow))
        check[f'{p["tool"]}|{p["setting"]}'] = {
            "per_residue_recall_mean_over_regions": float(np.mean(recall)),
            "n_regions_any_body_residue_covered": int(sum(x > 0 for x in recall)),
        }

    regions = []
    for i, r in enumerate(rows):
        nn = np.flatnonzero(~null[i])
        regions.append({
            "id": r["row_id"], "acc": r["accession"], "sym": r["symbol"], "fam": r["pfam_id"],
            "s": r["true_start"], "e": r["true_end"], "len": r["region_length"],
            "plen": r["protein_length"],
            "pid": None if r["best_pident"] is None or math.isnan(r["best_pident"]) else round(r["best_pident"], 1),
            "plddt": None if r["region_plddt"] is None else round(r["region_plddt"], 1),
            # First and last column inside the protein; None when the whole window is past
            # its end (a Pfam region placed beyond the end of this sequence).
            "c0": int(nn[0]) if nn.size else None, "c1": int(nn[-1]) if nn.size else None,
            "past_end": bool(r["true_end"] > r["protein_length"]),
            "seq": r["window_seq"], "off": r["window_offset"],
        })
    return {"regions": regions, "panels": out_panels, "flank": FLANK, "body": BODY}, check


if __name__ == "__main__":
    sys.exit("import build() from build_site.py")
