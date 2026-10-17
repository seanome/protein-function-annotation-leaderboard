#!/usr/bin/env python3
"""Score SCOPe40 per query, from the per-query files Foldseek's scoring script writes.

Writes data/scope40_long.parquet, one row per (tool, setting, level, query), and
data/scope40_summary.json, the AUC, bootstrap interval and coverage per tool and level, overall
and per SCOPe class. The page and its figures read only these two files.

Rules (Foldseek paper Methods and its script, bench.noselfhit.awk):
  family       TP = same family
  superfamily  TP = same superfamily, different family (same family: not counted)
  fold         TP = same fold, different superfamily (same superfamily: not counted)
  FP at every level = a different fold, or a target not in the domain list
Queries: SCOPe 2.08 ASTRAL 40% domains with another member at all three levels (5,737).
Denominators, as in LoL-align's copy of the script (soedinglab/LoLalign-analysis,
scope_benchmark/bench.noselfhit.awk, Reifenrath et al. 2026): the domains that can be found,
family size - 1, superfamily size - family size, fold size - superfamily size (column
`sensitivity`, `n_tp_total`). Foldseek's original script (steineggerlab/foldseek-analysis)
divides the superfamily and fold counts by one more, because it leaves the query itself in
them; that scale is kept as `sensitivity_foldseek_script` and `n_tp_total_foldseek_script`.
The published Foldseek and MMseqs2 files are on Foldseek's original scale: each row's true
positives are its value times that file's denominator (a whole number on every row), divided
again by the LoL-align denominator from the same row's family, superfamily and fold sizes.

A query missing from a tool's file scores 0 with has_hits false and is listed. The script
walks each query's hits in score order and stops counting at the first FP, so every hit
above that FP is a TP at exactly one level. So the 1-based index of the first FP, self-hit
excluded, is the sum of the three TP counts plus 1 when the file reports an FP, and a query
had at least one hit when that sum or the FP count is above 0.

AUC = mean sensitivity over all queries (the area under sensitivity plotted against the
fraction of queries, queries sorted highest first). Interval: 1,000 bootstrap resamples of
the queries, seed 20261017, 2.5th and 97.5th percentiles; every tool is resampled with the
same draws.
"""

import argparse
import collections
import hashlib
import json
from pathlib import Path

import numpy as np
import polars as pl

LEVELS = {"family": "FAM", "superfamily": "SFAM", "fold": "FOLD"}
SEED = 20261017
N_BOOT = 1000
# Published per-query files -> (tool id on this site, setting label).
FILES = {
    "foldseek.rocx": ("foldseek", "published run"),
    "mmseqs.rocx": ("mmseqs2_seqseq", "published run"),
}


def read_fasta(p: Path) -> dict[str, tuple[str, int]]:
    out, cur, n = {}, None, 0
    for line in p.open():
        if line.startswith(">"):
            if cur:
                out[cur[0]] = (cur[1], n)
            f = line[1:].split()
            cur, n = (f[0], f[1]), 0
        else:
            n += len(line.strip())
    if cur:
        out[cur[0]] = (cur[1], n)
    return out


def parts(s: str):
    p = s.split(".")
    return p[0], ".".join(p[:2]), ".".join(p[:3]), s


def queries(dom: dict) -> pl.DataFrame:
    fam = collections.Counter(parts(s)[3] for s, _ in dom.values())
    sf = collections.Counter(parts(s)[2] for s, _ in dom.values())
    fo = collections.Counter(parts(s)[1] for s, _ in dom.values())
    rows = []
    for d, (s, length) in dom.items():
        c, f, s3, s4 = parts(s)
        if fam[s4] > 1 and sf[s3] > fam[s4] and fo[f] > sf[s3]:
            rows.append(dict(query=d, scope_sccs=s, scope_class=c, fold=f, superfamily=s3, family=s4,
                             query_length=length, famcnt=fam[s4], sfamcnt=sf[s3], foldcnt=fo[f]))
    return pl.DataFrame(rows).sort("query")


def denominators(f, s, o):
    """(Foldseek's original script, LoL-align) per level, from sizes that include the query."""
    return {"family": (f - 1, f - 1), "superfamily": (s - (f - 1), s - f), "fold": (o - (s - 1), o - s)}


def score_file(path: Path, q: pl.DataFrame, tool: str, setting: str) -> pl.DataFrame:
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    rec = {r[0]: dict(zip(head, r)) for r in (x.split("\t") for x in lines[1:])}
    rows = []
    for x in q.iter_rows(named=True):
        r = rec.get(x["query"])
        if r is not None:
            f, s, o = int(r["FAMCNT"]), int(r["SFAMCNT"]), int(r["FOLDCNT"])
        else:
            f, s, o = x["famcnt"], x["sfamcnt"], x["foldcnt"]
        den = denominators(f, s, o)
        tp = {lvl: (round(float(r[c]) * den[lvl][0]) if r is not None else 0) for lvl, c in LEVELS.items()}
        n_fp = int(r["FP"]) if r is not None else 0
        n_above = sum(tp.values())
        has_hits = (n_above + n_fp) > 0
        first_fp = n_above + 1 if n_fp > 0 else None
        for lvl in LEVELS:
            d_script, d_lol = den[lvl]
            rows.append(dict(
                tool=tool, setting=setting, level=lvl, query=x["query"], scope_sccs=x["scope_sccs"],
                scope_class=x["scope_class"], fold=x["fold"], superfamily=x["superfamily"], family=x["family"],
                query_length=x["query_length"], n_tp_total=d_lol, n_tp_total_foldseek_script=d_script,
                n_tp_above_cut=tp[lvl], sensitivity=tp[lvl] / d_lol,
                sensitivity_foldseek_script=tp[lvl] / d_script,
                sensitivity_published=float(r[LEVELS[lvl]]) if r else None,
                has_hits=has_hits, first_fp_index=first_fp, in_published_file=r is not None,
                counts_match_208=(r is None) or (f, s, o) == (x["famcnt"], x["sfamcnt"], x["foldcnt"]),
                source_file=path.name))
    return pl.DataFrame(rows)


def read_scored(paths: list[Path], q: pl.DataFrame) -> pl.DataFrame:
    """Per-query rows written by score_hitlist.py (seanome/2024-kmerseek-analysis,
    nextflow-runs/scope40-hitlists) from a tool's own hit list, in this table's columns."""
    if not paths:
        return None
    def norm(df):
        # score_hitlist.py before 2026-10-17 called the LoL-align scale sensitivity_possible
        # and the Foldseek-script scale sensitivity; later versions name both.
        if "sensitivity_possible" in df.columns:
            df = df.rename({"sensitivity": "sensitivity_foldseek_script", "n_tp_total": "n_tp_total_foldseek_script",
                            "sensitivity_possible": "sensitivity", "n_tp_possible": "n_tp_total"})
        return df
    t = pl.concat([norm(pl.read_parquet(p)).with_columns(pl.lit(p.name).alias("source_file")) for p in paths],
                  how="diagonal")
    t = t.join(q.select("query", "query_length", pl.col("famcnt"), pl.col("sfamcnt"), pl.col("foldcnt")),
               on="query", how="inner")
    n = q.height
    bad = t.group_by("tool", "setting", "level").len().filter(pl.col("len") != n)
    if bad.height:
        raise SystemExit(f"scored files do not cover the {n} queries: {bad}")
    return t.select(
        "tool", "setting", "level", "query", "scope_sccs", "scope_class", "fold", "superfamily", "family",
        "query_length", pl.col("n_tp_total").cast(pl.Int64), pl.col("n_tp_total_foldseek_script").cast(pl.Int64),
        pl.col("n_tp_above_cut").cast(pl.Int64), "sensitivity", "sensitivity_foldseek_script",
        pl.lit(None, dtype=pl.Float64).alias("sensitivity_published"), "has_hits",
        pl.col("first_fp_index").cast(pl.Int64), pl.lit(True).alias("in_published_file"),
        pl.lit(True).alias("counts_match_208"), "source_file")


def boot(values: np.ndarray, idx: np.ndarray) -> tuple[float, float]:
    m = values[idx].mean(axis=1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def summarise(long: pl.DataFrame) -> dict:
    out = {"seed": SEED, "n_boot": N_BOOT, "overall": [], "by_class": []}
    rng = np.random.default_rng(SEED)
    qs = long.filter((pl.col("level") == "family")).select("query").unique().sort("query")["query"].to_list()
    n = len(qs)
    idx = rng.integers(0, n, size=(N_BOOT, n))
    classes = sorted(long["scope_class"].unique().to_list())
    cls_idx = {}
    for c in classes:
        rng_c = np.random.default_rng([SEED, ord(c)])
        nc = long.filter((pl.col("level") == "family") & (pl.col("scope_class") == c)
                         & (pl.col("tool") == long["tool"][0]))["query"].len()
        cls_idx[c] = rng_c.integers(0, nc, size=(N_BOOT, nc))
    for (tool, setting, level), g in long.group_by(["tool", "setting", "level"], maintain_order=True):
        g = g.sort("query")
        assert g["query"].to_list() == qs
        v = g["sensitivity"].to_numpy()
        h = g["has_hits"].to_numpy()
        lo, hi = boot(v, idx)
        out["overall"].append(dict(
            tool=tool, setting=setting, level=level, n_queries=n, auc=float(v.mean()), lo=lo, hi=hi,
            coverage=float(h.mean()), n_with_hits=int(h.sum()), auc_hits_only=float(v[h].mean()),
            auc_foldseek_script=float(g["sensitivity_foldseek_script"].to_numpy().mean()),
            n_missing_from_file=int((~g["in_published_file"]).sum())))
        for c in classes:
            gc = g.filter(pl.col("scope_class") == c)
            vc = gc["sensitivity"].to_numpy()
            lo_c, hi_c = boot(vc, cls_idx[c])
            out["by_class"].append(dict(tool=tool, setting=setting, level=level, scope_class=c,
                                        n_queries=len(vc), mean=float(vc.mean()), lo=lo_c, hi=hi_c))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--astral208", type=Path, required=True)
    ap.add_argument("--rocx-dir", type=Path, required=True)
    ap.add_argument("--scored", type=Path, nargs="*", default=[],
                    help="per_query.parquet files from score_hitlist.py (e.g. the kmerseek pair tables)")
    ap.add_argument("--out", type=Path, default=Path("data"))
    a = ap.parse_args()
    dom = read_fasta(a.astral208)
    q = queries(dom)
    long = pl.concat([score_file(a.rocx_dir / f, q, *ts) for f, ts in FILES.items()])
    long = long.with_columns(pl.col("n_tp_total").cast(pl.Int64), pl.col("n_tp_total_foldseek_script").cast(pl.Int64),
                             pl.col("n_tp_above_cut").cast(pl.Int64), pl.col("first_fp_index").cast(pl.Int64))
    extra = read_scored(sorted(a.scored), q)
    if extra is not None:
        long = pl.concat([long, extra.select(long.columns)])
    long.write_parquet(a.out / "scope40_long.parquet")
    s = summarise(long)
    s["n_domains"] = len(dom)
    s["n_queries"] = q.height
    s["missing"] = {t: sorted(long.filter((pl.col("tool") == t) & ~pl.col("in_published_file")
                                          & (pl.col("level") == "family"))["query"].to_list())
                    for t in long["tool"].unique().sort().to_list()}
    s["inputs"] = {"astral208_sha256": hashlib.sha256(a.astral208.read_bytes()).hexdigest(),
                   **{f: hashlib.sha256((a.rocx_dir / f).read_bytes()).hexdigest() for f in FILES},
                   **{p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in a.scored}}
    (a.out / "scope40_summary.json").write_text(json.dumps(s, indent=1) + "\n")
    print(f"{long.height:,} rows; {q.height:,} queries")
    for r in s["overall"]:
        print(f"{r['tool']:16} {r['level']:12} AUC {r['auc']:.3f} [{r['lo']:.3f}, {r['hi']:.3f}] "
              f"coverage {r['coverage']:.3f} hits-only {r['auc_hits_only']:.3f} Foldseek script {r['auc_foldseek_script']:.3f}")


if __name__ == "__main__":
    main()
