"""kmerseek as a fourth label: for each alphabet and query set, how many proteins get a
region to another Swiss-Prot entry at E <= 0.001, and, on the lit oyster control, how
often that region's target shares an InterPro entry with the query.

A dark protein's hit cannot be checked: it has no InterPro entry by definition. The lit
control can, so its share of right hits is printed beside every dark number. The same
check on MMseqs2's best hit to another entry is the reference.

    python3 -I scripts/why/kmerseek_labels.py --raw ~/data/why-numbers/2026_03 \
        --hits ~/data/why-numbers/2026_03/kmerseek_hits
--hits holds the <chunk>.<alphabet>...hits.csv.gz (and .nofit) files the search jobs write.
Writes results/why/kmerseek/best_hits.parquet, results/why/kmerseek/summary.parquet, and
step 7 rows of results/why/numbers.parquet.
"""

import argparse
import collections
import gzip
import math
import pathlib
import re
import sys

import polars as pl

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from numbers_table import Numbers  # noqa: E402

OUT = pathlib.Path("results/why/kmerseek")
EVALUE = 1e-3  # the dark-set pipeline's evalue_call: what counts as "found something"
SETS = ["four", "oyster_dark", "oyster_lit_control", "human_dark"]


def sprot_interpro(raw):
    """Swiss-Prot accession -> its InterPro entries, from the DR InterPro lines of the release."""
    cache = raw / "kmerseek_inputs" / "sprot_interpro.parquet"
    if cache.exists():
        return pl.read_parquet(cache)
    rows, acc = [], None
    with gzip.open(raw / "uniprot_sprot.dat.gz", "rt") as fh:
        for line in fh:
            if line.startswith("AC   ") and acc is None:
                acc = line[5:].split(";")[0].strip()
            elif line.startswith("DR   InterPro;"):
                rows.append((acc, line.split(";")[1].strip()))
            elif line.startswith("//"):
                acc = None
    df = pl.DataFrame(rows, schema=["accession", "ipr"], orient="row")
    df.write_parquet(cache)
    return df


def query_interpro(raw):
    """Oyster and human accession -> InterPro entries, as UniProt cross-references them."""
    out = []
    for tag in ("human", "oyster"):
        x = pl.read_csv(raw / f"{tag}_interpro_xrefs.tsv.gz", separator="\t", infer_schema_length=0)
        out.append(x.select(accession="Entry", ipr=pl.col("InterPro").str.split(";")).explode("ipr")
                   .filter(pl.col("ipr").is_not_null(), pl.col("ipr") != ""))
    return pl.concat(out).unique()


def query_sequences(raw):
    seqs, acc = {}, None
    for f in sorted((raw / "kmerseek_inputs" / "queries").glob("chunk_*.fasta")):
        for line in open(f):
            if line.startswith(">"):
                acc = line[1:].strip()
                seqs[acc] = []
            else:
                seqs[acc].append(line.strip())
    return {a: "".join(v) for a, v in seqs.items()}


def entropy_bits(s):
    """Shannon entropy of the residues in a stretch, in bits: 0 for one residue repeated,
    about 4.2 for a typical protein stretch, at most log2(20) = 4.32."""
    c = collections.Counter(s)
    n = len(s)
    return -sum(v / n * math.log2(v / n) for v in c.values()) if n else None


def read_hits(hits_dir):
    frames, nofit = [], []
    for f in sorted(hits_dir.glob("*.hits.csv.gz")):
        # Only finished searches: a running one may be writing its hits file right now.
        if not f.with_name(f.name.replace(".hits.csv.gz", ".DONE")).exists():
            continue
        m = re.match(r"(chunk_\d+)\.([^.]+)\.k(\d+)\.", f.name)
        df = pl.read_csv(f, schema_overrides={"query_name": pl.String, "target_name": pl.String})
        frames.append(df.with_columns(chunk=pl.lit(m.group(1)), alphabet=pl.lit(m.group(2)), ksize=pl.lit(int(m.group(3)))))
    for f in sorted(hits_dir.glob("*.nofit")):
        m = re.match(r"(chunk_\d+)\.([^.]+)\.k(\d+)\.", f.name)
        nofit.append({"chunk": m.group(1), "alphabet": m.group(2), "ksize": int(m.group(3))})
    return (pl.concat(frames) if frames else pl.DataFrame(schema={"query_name": pl.String, "target_name": pl.String, "chunk": pl.String, "alphabet": pl.String})), pl.DataFrame(nofit, schema={"chunk": pl.String, "alphabet": pl.String, "ksize": pl.Int64})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=pathlib.Path)
    ap.add_argument("--hits", required=True, type=pathlib.Path)
    ap.add_argument("--queries", default="results/why/kmerseek/queries.tsv")
    ap.add_argument("--settings", default="results/why/kmerseek/settings.tsv")
    ap.add_argument("--allow-missing", action="store_true", help="for a test run: not every alphabet and chunk")
    ap.add_argument("--out", type=pathlib.Path, default=OUT)
    args = ap.parse_args()

    q = pl.read_csv(args.queries, separator="\t")
    settings = pl.read_csv(args.settings, separator="\t")
    hits, nofit = read_hits(args.hits)
    done = hits.select("chunk", "alphabet").unique().vstack(nofit.select("chunk", "alphabet"))
    want = settings.select("alphabet").join(q.select("chunk").unique(), how="cross")
    missing = want.join(done, on=["chunk", "alphabet"], how="anti")
    if missing.height and not args.allow_missing:
        raise ValueError(f"{missing.height} alphabet x chunk searches have no output:\n{missing}")
    # A partial run scores only alphabets with every chunk done: a share over some chunks
    # would be a share of a different set of proteins.
    incomplete = missing["alphabet"].unique().to_list()
    settings = settings.filter(~pl.col("alphabet").is_in(incomplete))
    if incomplete:
        print("not scored yet (chunks missing):", ", ".join(sorted(incomplete)))

    # Best region per query and alphabet, to an entry other than the query itself.
    best = (hits.filter(pl.col("query_name") != pl.col("target_name"))
            .sort(["alphabet", "query_name", "region_evalue", "region_ka_bits"], descending=[False, False, False, True])
            .group_by(["alphabet", "query_name"], maintain_order=True).first()
            .rename({"query_name": "accession"}))
    best = best.join(q.select("accession", "query_set"), on="accession", how="left")
    if best["query_set"].null_count():
        raise ValueError("a hit's query is not in queries.tsv")
    t_ipr = sprot_interpro(args.raw).group_by("accession").agg(target_ipr=pl.col("ipr"))
    q_ipr = query_interpro(args.raw).group_by("accession").agg(query_ipr=pl.col("ipr"))
    best = (best.join(t_ipr.rename({"accession": "target_name"}), on="target_name", how="left")
            .join(q_ipr, on="accession", how="left")
            .with_columns(labelled=pl.col("region_evalue") <= EVALUE,
                          shares_interpro=pl.col("query_ipr").list.set_intersection(pl.col("target_ipr")).list.len().fill_null(0) > 0))
    # How varied the matched query stretch is. A hit on a run of one or two residues
    # (EEEEE..., GGLGGGL...) is composition, not shared ancestry.
    seqs = query_sequences(args.raw)
    best = best.with_columns(region_residues=pl.struct("accession", "region_start", "region_end").map_elements(
        lambda r: seqs[r["accession"]][r["region_start"]:r["region_end"]], return_dtype=pl.String))
    best = best.with_columns(region_entropy_bits=pl.col("region_residues").map_elements(entropy_bits, return_dtype=pl.Float64))
    ctrl_floor = best.filter(pl.col("labelled"), pl.col("query_set") == "oyster_lit_control")["region_entropy_bits"].quantile(0.05)
    best = best.with_columns(varied_region=pl.col("region_entropy_bits") >= ctrl_floor)
    print(f"5th percentile of residue entropy in labelled lit-control stretches: {ctrl_floor:.2f} bits")
    best.write_parquet(args.out / "best_hits.parquet")

    # MMseqs2's best hit to another entry, checked the same way, on the same lit control.
    pl_tab = pl.read_parquet("results/why/proteome_labels.parquet").select("accession", "sprot_best_nonself_target")
    ctrl = q.filter(pl.col("query_set") == "oyster_lit_control").join(pl_tab, on="accession", how="left").join(q_ipr, on="accession", how="left")
    ctrl = ctrl.join(t_ipr.rename({"accession": "sprot_best_nonself_target"}), on="sprot_best_nonself_target", how="left").with_columns(
        mm_right=pl.col("query_ipr").list.set_intersection(pl.col("target_ipr")).list.len().fill_null(0) > 0)
    ctrl_checkable = ctrl.filter(pl.col("query_ipr").is_not_null())
    mm_labelled = ctrl_checkable.filter(pl.col("sprot_best_nonself_target").is_not_null())

    N = Numbers()
    rel = "kmerseek 0.4.0 (PR 141 build 3c353d4) against Swiss-Prot 2026_03"
    sizes = dict(q.group_by("query_set").agg(pl.len()).iter_rows())
    n_checkable = ctrl_checkable.height
    if n_checkable:  # a test run has no lit control
        N.add("km_control_checkable", 7, "lit oyster control proteins with at least one InterPro entry, so a hit can be checked", n_checkable, "proteins", "swissprot_dat", rel)
        N.add("km_mmseqs_control_labelled", 7, "of those, proteins with an MMseqs2 hit above 30% identity to another Swiss-Prot entry", mm_labelled.height, "proteins", "sherlock_mmseqs", rel)
        N.add("km_mmseqs_control_right_pct", 7, "share of those MMseqs2 best hits whose target shares an InterPro entry with the query", round(100 * mm_labelled["mm_right"].mean(), 2), "% of labelled proteins", "sherlock_mmseqs", rel)
    rows = []
    for a in settings["alphabet"]:
        b = best.filter(pl.col("alphabet") == a, pl.col("labelled"))
        k = int(settings.filter(pl.col("alphabet") == a)["ksize"][0])
        row = {"alphabet": a, "ksize": k}
        for s in SETS:
            if s not in sizes:  # a test run searches only some sets
                continue
            n = b.filter(pl.col("query_set") == s).height
            row[f"{s}_labelled"] = n
            row[f"{s}_labelled_pct"] = round(100 * n / sizes[s], 2)
        for s in [x for x in ("oyster_dark", "human_dark") if x in sizes]:
            dv = b.filter(pl.col("query_set") == s)
            row[f"{s}_labelled_varied"] = dv["varied_region"].sum()
            row[f"{s}_entropy_median"] = dv["region_entropy_bits"].median()
        row["oyster_lit_control_entropy_median"] = b.filter(pl.col("query_set") == "oyster_lit_control")["region_entropy_bits"].median()
        cb = b.filter(pl.col("query_set") == "oyster_lit_control", pl.col("query_ipr").is_not_null())
        row["control_checkable_labelled"] = cb.height
        row["control_right_pct"] = round(100 * cb["shares_interpro"].mean(), 2) if cb.height else None
        row["nofit_chunks"] = nofit.filter(pl.col("alphabet") == a).height
        rows.append(row)
        for s in [x for x in ("oyster_dark", "oyster_lit_control", "human_dark") if x in sizes]:
            N.add(f"km_{a}_{s}_labelled_pct", 7, f"kmerseek {a} k{k}: share of {s.replace('_', ' ')} proteins with a region to another Swiss-Prot entry at E <= 0.001",
                  row[f"{s}_labelled_pct"], f"% of {sizes[s]} proteins", "sherlock_kmerseek", rel)
        if "oyster_dark" in sizes:
            N.add(f"km_{a}_oyster_dark_labelled_varied", 7, f"kmerseek {a} k{k}: dark oyster proteins labelled on a stretch as varied as 95% of the lit-control stretches (residue entropy >= {ctrl_floor:.2f} bits)",
                  row["oyster_dark_labelled_varied"], "proteins", "sherlock_kmerseek", rel)
        if row["control_right_pct"] is not None:
            N.add(f"km_{a}_control_right_pct", 7, f"kmerseek {a} k{k}: share of its labelled, checkable lit-control proteins whose best target shares an InterPro entry",
                  row["control_right_pct"], "% of labelled proteins", "sherlock_kmerseek", rel)
    summary = pl.DataFrame(rows)
    summary.write_parquet(args.out / "summary.parquet")
    if not args.allow_missing:  # a test run does not touch the committed table
        N.write("results/why/numbers.parquet", replace_steps={7})
    print(summary.to_pandas().to_string())


if __name__ == "__main__":
    main()
