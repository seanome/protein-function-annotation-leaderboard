"""Second, independent recheck of the step-5 numbers (dark vs lit proteins).

Reads only results/why/proteome_labels.parquet and results/why/numbers.parquet.
Writes one verdict row per step-5 id into results/why/verification/verdicts.parquet,
replacing rows with those ids and leaving all other rows as they were.
Also prints checks of the figure claims. Run from the repository root:

    python -I results/why/verification/recompute_dark.py
"""

from pathlib import Path

import polars as pl

LABELS = Path("results/why/proteome_labels.parquet")
NUMBERS = Path("results/why/numbers.parquet")
VERDICTS = Path("results/why/verification/verdicts.parquet")

SCHEMA = {
    "id": pl.String,
    "stated_value": pl.Float64,
    "recomputed_value": pl.Float64,
    "verdict": pl.String,
    "reason": pl.String,
}

d = pl.read_parquet(LABELS)

# Sanity checks on the inputs the rule depends on.
for col in ["species", "accession", "length", "interpro_match", "sprot_hit", "af_model", "busco_status"]:
    nulls = d[col].null_count()
    print(f"nulls in {col}: {nulls}")
print("duplicate (species, accession):", d.height - d.select("species", "accession").unique().height)
print("busco_status values:", sorted(d["busco_status"].unique().to_list()))

d = d.with_columns(
    dark=(~pl.col("interpro_match")) & (~pl.col("sprot_hit")),
)


def pct(mask: pl.Series) -> float:
    return round(100.0 * mask.sum() / mask.len(), 2)


recomputed: dict[str, float] = {}
for species in ["oyster", "human"]:
    for dark in [True, False]:
        g = d.filter((pl.col("species") == species) & (pl.col("dark") == dark))
        tag = f"{species}_{'dark' if dark else 'lit'}"
        recomputed[f"{tag}_n"] = float(g.height)
        recomputed[f"{tag}_length_median"] = float(g["length"].cast(pl.Float64).median())
        recomputed[f"{tag}_under_100aa_pct"] = pct(g["length"] < 100)
        # Denominator: every protein in the group; null pLDDT (no model) is not <= 50.
        recomputed[f"{tag}_plddt_50_or_less_pct"] = pct((g["af_mean_plddt"] <= 50).fill_null(False))
        recomputed[f"{tag}_no_model_pct"] = pct(~g["af_model"])
        recomputed[f"{tag}_busco_fragmented_pct"] = pct(g["busco_status"] == "Fragmented")

stated = pl.read_parquet(NUMBERS).filter(pl.col("step") == 5).select("id", "value")

rows = []
for sid, sval in stated.iter_rows():
    if sid not in recomputed:
        rows.append((sid, sval, None, "UNVERIFIABLE", "id not covered by the step-5 rule"))
        continue
    r = recomputed[sid]
    ok = abs(r - sval) < 1e-9
    reason = (
        "recomputed from proteome_labels.parquet with dark = not interpro_match and not sprot_hit"
        if ok
        else f"recomputed {r} differs from stated {sval}"
    )
    rows.append((sid, float(sval), r, "HOLDS" if ok else "WRONG", reason))

new = pl.DataFrame(rows, schema=SCHEMA, orient="row")
if VERDICTS.exists():
    old = pl.read_parquet(VERDICTS)
    old = old.filter(~pl.col("id").is_in(new["id"].to_list()))
    out = pl.concat([old.select(list(SCHEMA)).cast(SCHEMA), new])
else:
    out = new
out.write_parquet(VERDICTS)

with pl.Config(tbl_rows=50, fmt_str_lengths=40):
    print(new.select("id", "stated_value", "recomputed_value", "verdict"))
print(f"verdicts.parquet now has {out.height} rows")

# ---- Figure claims ----
print("\n== Claim 1: pLDDT bands among dark proteins ==")


def band(p):
    return (
        pl.when(p.is_null()).then(pl.lit("no model"))
        .when(p > 90).then(pl.lit(">90"))
        .when(p > 70).then(pl.lit("(70,90]"))
        .when(p > 50).then(pl.lit("(50,70]"))
        .otherwise(pl.lit("<=50"))
    )


for species in ["oyster", "human"]:
    g = (
        d.filter((pl.col("species") == species) & pl.col("dark"))
        .sort(["length", "accession"], descending=[True, False])
        .with_columns(band=band(pl.col("af_mean_plddt")))
    )
    counts = dict(g.group_by("band").len().iter_rows())
    print(species, {k: counts.get(k, 0) for k in [">90", "(70,90]", "(50,70]", "<=50", "no model"]})
    # Same bands but using af_model for "no model"
    alt = g.with_columns(
        band2=pl.when(~pl.col("af_model")).then(pl.lit("no model")).otherwise(pl.col("band"))
    )
    c2 = dict(alt.group_by("band2").len().iter_rows())
    print(species, "(no model = af_model false)", {k: c2.get(k, 0) for k in [">90", "(70,90]", "(50,70]", "<=50", "no model"]})

    print("\n== Claim 2: length range ==")
    lengths = g["length"].to_list()
    print(species, "longest", lengths[0], "shortest", lengths[-1])
    if species == "oyster":
        print("1000th..5000th longest:", [lengths[i - 1] for i in (1000, 2000, 3000, 4000, 5000)])

print("\n== Claim 3: oyster shares ==")
o = d.filter(pl.col("species") == "oyster")
print("af_model true %:", pct(o["af_model"]))
print("dark %:", pct(o["dark"]))
print("none of af_model, interpro_match, sprot_hit %:",
      pct(~o["af_model"] & ~o["interpro_match"] & ~o["sprot_hit"]))

print("\n== Claim 4: af_mean_plddt null exactly when af_model false ==")
x = d.group_by("species", "af_model", plddt_null=pl.col("af_mean_plddt").is_null()).len().sort("species", "af_model")
print(x)
