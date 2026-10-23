"""The kmerseek settings for the four proteins and the oyster and human proteins: one k per
alphabet, its search flags, and the memory and time each Sherlock task asks for.

- k: the k with the best Pfam Fmax on placed proteins (proteins a sequence tool already
  places), median over the four species of analysis notebook 291, scaled 1, no
  low-complexity mask, extension at the alphabet's own mismatch penalty, transfers not
  deduplicated (the reading that gives notebook 291's hp_pbotc_1st_ed2 k19 at 0.259).
  Ties go to the larger k, which searches with less memory; gbmr7 scores 0 at every k.
- Mismatch penalty C and X-drop: the dark-set pipeline's rule, C = c_opt from
  kappa_by_alphabet.tsv, X = 4 x C, both written as the shortest decimal string.
- Karlin-Altschul fit queries: 500 up to 30 bits per k-mer (k x log2 classes), doubling
  every 3 bits above, capped at 3000, rounded to 50 (the pipeline's kaQueriesFor).
- Search memory and time: the largest peak and longest run of that setting's searches in
  the dark-set 0.4 midi run (2_000-protein chunks against a 572k-entry Swiss-Prot), times
  1.5 for memory and 2 for time.

    python3 -I scripts/why/kmerseek_plan.py --raw ~/data/why-numbers/2026_03
Writes results/why/kmerseek/settings.tsv.
"""

import argparse
import math
import pathlib
import re

import polars as pl

OUT = pathlib.Path("results/why/kmerseek")
SPECIES_N = 4
INDEX_MEM_GB = 64  # pipeline rule 1.4 x max(45 GB build, fit model), floor 48; fits measured 36-47 GB
INDEX_TIME_H = 8
CPUS = 8  # the pipeline's withName: 'kmerseekIndex|kmerseekSearch' { cpus = 8 }


def penalty_string(c):
    s = f"{c:.2f}"
    return s.rstrip("0").rstrip(".") if "." in s else s


def ka_queries(bits):
    if bits <= 30:
        return 500
    n = 500 * 2 ** ((bits - 30) / 3)
    return int(round(min(3000, n) / 50) * 50)


def gb(x):
    if x in (None, "-", ""):
        return None
    v, u = x.split()
    return float(v) * {"B": 1e-9, "KB": 1e-6, "MB": 1e-3, "GB": 1, "TB": 1e3}[u]


def hours(x):
    if x in (None, "-", ""):
        return None
    return sum(float(n) * {"d": 24, "h": 1, "m": 1 / 60, "s": 1 / 3600, "ms": 1 / 3.6e6}[u]
               for n, u in re.findall(r"([\d.]+)(ms|h|m|s|d)", x))


def best_k(metrics):
    d = pl.read_csv(metrics, infer_schema_length=10000).filter(
        pl.col("tool") == "kmerseek", pl.col("truth_set") == "placed", pl.col("dedup_transfers") == False  # noqa: E712
    )
    d = d.with_columns(p=pl.col("variant").str.split(".")).with_columns(
        alphabet=pl.col("p").list.get(0),
        ksize=pl.col("p").list.get(1).str.slice(1).cast(pl.Int64),
        scaled=pl.col("p").list.get(2),
        lowcomp=pl.col("p").list.get(3),
        arm=pl.col("p").list.get(4),
    ).filter(pl.col("scaled") == "s1", pl.col("lowcomp") == "lcfalse", pl.col("arm").str.starts_with("extend"))
    g = d.group_by("alphabet", "ksize").agg(
        fmax_median=pl.col("fmax").median(), n_species=pl.col("species").n_unique()
    ).filter(pl.col("n_species") == SPECIES_N)
    return g.sort(["alphabet", "fmax_median", "ksize"], descending=[False, True, True]).group_by(
        "alphabet", maintain_order=True
    ).first()


def measured_searches(trace):
    t = pl.read_csv(trace, separator="\t", infer_schema_length=0, quote_char=None).filter(
        pl.col("name").str.contains("kmerseekSearch"), pl.col("status").is_in(["COMPLETED", "CACHED"])
    )
    t = t.with_columns(tag=pl.col("name").str.extract(r"\((.*)\)")).with_columns(
        alphabet=pl.col("tag").str.extract(r"chunk_\d+\.([^.]+)\."),
        ksize=pl.col("tag").str.extract(r"\.k(\d+)\.").cast(pl.Int64),
        scaled=pl.col("tag").str.extract(r"\.(s\d+)\."),
        arm=pl.col("tag").str.extract(r"\.(exact|extend-c[\d.]+)$"),
        peak_gb=pl.col("peak_rss").map_elements(gb, return_dtype=pl.Float64),
        run_h=pl.col("realtime").map_elements(hours, return_dtype=pl.Float64),
    ).filter(pl.col("scaled") == "s1", pl.col("arm").str.starts_with("extend"))
    return t.group_by("alphabet", "ksize").agg(
        n_measured=pl.len(), peak_gb_max=pl.col("peak_gb").max(), run_h_max=pl.col("run_h").max()
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=pathlib.Path)
    args = ap.parse_args()
    inp = args.raw / "kmerseek_inputs"
    kappa = pl.read_csv(inp / "kappa_by_alphabet.tsv", separator="\t", comment_prefix="#")
    plan = best_k(inp / "dark_fmax_metrics_2026-10-14.csv").join(
        kappa.select("alphabet", "classes", "c_opt"), on="alphabet", how="left"
    ).join(measured_searches(inp / "trace.ladder-0.4-midi.txt"), on=["alphabet", "ksize"], how="left")
    if plan["c_opt"].null_count() or plan["n_measured"].null_count():
        raise ValueError(plan.filter(pl.col("c_opt").is_null() | pl.col("n_measured").is_null()))
    rows = []
    for r in plan.sort("alphabet").iter_rows(named=True):
        bits = r["ksize"] * math.log2(r["classes"])
        c = penalty_string(r["c_opt"])
        rows.append({
            "alphabet": r["alphabet"], "ksize": r["ksize"], "classes": r["classes"], "bits_per_kmer": round(bits, 2),
            "fmax_median_placed": round(r["fmax_median"], 4), "mismatch_penalty": c,
            "xdrop": penalty_string(4 * float(c)), "ka_queries": ka_queries(bits),
            "measured_searches": r["n_measured"], "measured_peak_gb": round(r["peak_gb_max"], 1),
            "measured_hours": round(r["run_h_max"], 2),
            "search_mem_gb": max(32, 8 * math.ceil(1.5 * r["peak_gb_max"] / 8)),
            "search_time_h": max(2, math.ceil(2 * r["run_h_max"])),
            # More fit queries take longer: the 500-query fits took up to 43 minutes.
            "index_mem_gb": INDEX_MEM_GB, "index_time_h": INDEX_TIME_H if ka_queries(bits) == 500 else 12, "cpus": CPUS,
        })
    OUT.mkdir(parents=True, exist_ok=True)
    out = pl.DataFrame(rows)
    out.write_csv(OUT / "settings.tsv", separator="\t")
    print(out.to_pandas().to_string())


if __name__ == "__main__":
    main()
