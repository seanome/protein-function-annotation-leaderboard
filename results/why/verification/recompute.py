"""Independent recomputation of results/why/numbers.parquet from the raw files.

Own code only (nothing from scripts/why/). Run parse_raw.py first: it caches the
Swiss-Prot .dat and taxdump as parquet in ./cache. Uses taxlib.py (own helper).
Writes verdicts.parquet next to this file.
"""
import gzip
import hashlib
import json
import statistics
import sys
from pathlib import Path

import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from taxlib import Tax  # noqa: E402

R = Path("/Users/olga/data/why-numbers/2026_03")
S = R / "sherlock"
C = HERE / "cache"
W = HERE.parent

# ---------------------------------------------------------------- sha256
src = json.load(open(W / "sources.json"))
path_of = {k: R / v["file"] for k, v in src.items()}
path_of["sherlock_mmseqs"] = S / "mmseqs_vs_sprot.tsv.gz"
path_of["sherlock_busco_human"] = S / "busco_human/run_metazoa_odb12/full_table.tsv"
path_of["sherlock_busco_oyster"] = S / "busco_oyster/run_metazoa_odb12/full_table.tsv"
sha_ok = {}
for k, v in src.items():
    h = hashlib.sha256(open(path_of[k], "rb").read()).hexdigest()
    sha_ok[k] = h == v["sha256"]
print("sha256 mismatches:", [k for k, ok in sha_ok.items() if not ok])

T = Tax()
V = {}  # id -> (recomputed value, extra reason)

# ---------------------------------------------------------------- step 1
sp = pl.read_parquet(C / "sprot.parquet")
oxs = sp["ox"].unique().to_list()
sp = sp.with_columns(
    pl.col("ox").replace_strict({t: T.species_of(t) for t in oxs}, default=None).alias("species"),
    pl.col("ox").replace_strict({t: T.col_group(t) for t in oxs}, default=None).alias("grp"),
)
sp = sp.with_columns(pl.col("ox").replace_strict({t: (T.resolve(t) or t) for t in oxs}, default=None).alias("taxid_ncbi"))
N = sp.height


def top5_share(df, key, denom):
    g = df.drop_nulls(key).group_by(key).len().sort("len", descending=True)
    return 100 * g.head(5)["len"].sum() / denom


V["sp_entries"] = N
V["sp_entries_no_species"] = sp["species"].is_null().sum()
V["sp_organisms"] = sp["taxid_ncbi"].n_unique()
raw_org = sp["ox"].n_unique()
V["sp_species"] = sp["species"].drop_nulls().n_unique()
V["sp_top5_share_by_organisms"] = top5_share(sp, "taxid_ncbi", N)
V["sp_top5_share_by_species"] = top5_share(sp, "species", N)
meta = sp.filter(pl.col("grp") == "Animalia")
V["sp_metazoa_entries"] = meta.height
V["sp_metazoa_species"] = meta["species"].drop_nulls().n_unique()

names = {T.name[s]: s for s in meta["species"].drop_nulls().unique().to_list()}
cnt = dict(meta.drop_nulls("species").group_by("species").len().rows())
matched, marine, marine_only = 0, set(), set()
for line in gzip.open(R / "worms_matches.jsonl.gz", "rt"):
    r = json.loads(line)
    ms = [m for m in r["matches"] if m.get("scientificname") == r["query"]]
    if not ms:
        continue
    m = ([x for x in ms if x.get("status") == "accepted"] or ms)[0]
    matched += 1
    s = names[r["query"]]
    if m.get("isMarine") == 1:
        marine.add(s)
        if m.get("isFreshwater") == 0 and m.get("isTerrestrial") == 0:
            marine_only.add(s)
V["worms_matched_species"] = matched
V["worms_marine_species"] = len(marine)
V["sp_marine_animal_entries"] = sum(cnt[s] for s in marine)
V["sp_marine_animal_share"] = 100 * V["sp_marine_animal_entries"] / N
V["sp_marine_only_animal_share"] = 100 * sum(cnt[s] for s in marine_only) / N

# ---------------------------------------------------------------- step 2
pdb = sp.filter(pl.col("n_pdb") > 0)
V["sp_pdb_entries"] = pdb.height
V["sp_pdb_organisms"] = pdb["taxid_ncbi"].n_unique()
raw_pdb_org = pdb["ox"].n_unique()
V["sp_pdb_species"] = pdb["species"].drop_nulls().n_unique()
V["sp_pdb_top5_share_by_organisms"] = top5_share(pdb, "taxid_ncbi", pdb.height)
V["sp_pdb_top5_share_by_species"] = top5_share(pdb, "species", pdb.height)
alt_sp_pdb_top5 = top5_share(pdb, "species", pdb["species"].drop_nulls().len())

kb = pl.read_csv(R / "uniprotkb_with_pdb.tsv.gz", separator="\t", infer_schema_length=0)
kb = kb.with_columns(pl.col("Organism (ID)").cast(pl.Int64).alias("ox"))
ids = kb["ox"].unique().to_list()
kb = kb.with_columns(
    pl.col("ox").replace_strict({t: T.cellular(t) for t in ids}, default=None).alias("cell"),
    pl.col("ox").replace_strict({t: T.species_of(t) for t in ids}, default=None).alias("species"),
)
kc = kb.filter(pl.col("cell"))
V["kb_pdb_entries_all"] = kb.height
V["kb_pdb_entries"] = kc.height
V["kb_pdb_organism_name"] = kc["Organism"].n_unique()
kc = kc.with_columns(pl.col("ox").replace_strict({t: (T.resolve(t) or t) for t in ids}, default=None).alias("taxid_ncbi"))
V["kb_pdb_taxid_ncbi"] = kc["taxid_ncbi"].n_unique()
raw_kb_tax = kc["ox"].n_unique()
V["kb_pdb_species_taxid"] = kc["species"].drop_nulls().n_unique()
V["kb_pdb_top5_share_by_organism_name"] = top5_share(kc, "Organism", kc.height)
V["kb_pdb_top5_share_by_taxid_ncbi"] = top5_share(kc, "taxid_ncbi", kc.height)
V["kb_pdb_top5_share_by_species_taxid"] = top5_share(kc, "species", kc.height)
alt_kb_sp_top5 = top5_share(kc, "species", kc["species"].drop_nulls().len())

dk = pl.read_parquet(R / "wealth_inequality_uniprot_combined_sprot_trembl_taxonomy_summary.parquet")
g = (dk.filter(pl.col("superdomain") == "Cellular Life").group_by("organism")
     .agg(pl.col("pdb_structures_count").sum()).filter(pl.col("pdb_structures_count") > 0)
     .sort("pdb_structures_count", descending=True))
deck_e = g["pdb_structures_count"].sum()
deck_o = g.height
deck_t = 100 * g.head(5)["pdb_structures_count"].sum() / deck_e
V["deck_entries_recomputed"] = deck_e
V["deck_organisms_recomputed"] = deck_o
V["deck_top5_pct_recomputed"] = deck_t
V["deck_entries_now_minus_then"] = V["kb_pdb_entries"] - deck_e
V["deck_organisms_now_minus_then"] = V["kb_pdb_organism_name"] - deck_o
V["deck_top5_pct_now_minus_then"] = V["kb_pdb_top5_share_by_organism_name"] - deck_t

sifts = {l.split("\t")[0] for l in gzip.open(R / "sifts_uniprot_pdb.tsv.gz", "rt")
         if not l.startswith("#") and not l.startswith("SP_PRIMARY")}
V["sifts_accessions"] = len(sifts)
V["sifts_and_uniprot"] = len(sifts & set(kb["Entry"].to_list()))

# ---------------------------------------------------------------- step 3
col = {c["col_group"]: c["accepted_species"] for c in json.load(open(R / "col_accepted_species.json"))["counts"]}
spg = dict(sp.drop_nulls(["species", "grp"]).group_by("grp").agg(pl.col("species").n_unique()).rows())
for grp in ["Animalia", "Archaea", "Bacteria", "Chromista", "Fungi", "Plantae", "Protozoa"]:
    V[f"col_species_{grp}"] = col[grp]
    V[f"sp_species_{grp}"] = spg[grp]
    V[f"sp_species_share_{grp}"] = 100 * spg[grp] / col[grp]
nb = json.load(open(R / "wealth_inequality_03_make_predicted_n_proteins.ipynb"))
cell6 = "".join(nb["cells"][6]["source"])
larsen = {}
for k in ["Animals", "Plants", "Fungi", "Protists", "Bacteria"]:
    v = cell6.split(f'"{k}":')[1].split(",")[0].strip()
    larsen[k] = float(v)
for k, v in larsen.items():
    V[f"larsen2017_estimated_species_{k}"] = v

# ---------------------------------------------------------------- step 4
def fasta(path):
    seqs, cur = {}, None
    for l in gzip.open(path, "rt"):
        if l.startswith(">"):
            cur = l[1:].split("|")[1]
            seqs[cur] = []
        else:
            seqs[cur].append(l.strip())
    return {k: "".join(v) for k, v in seqs.items()}


af = {}
for l in gzip.open(R / "afdb_predictions.jsonl.gz", "rt"):
    r = json.loads(l)
    af[r["accession"]] = r.get("models") or []
mm = pl.read_csv(S / "mmseqs_vs_sprot.tsv.gz", separator="\t", has_header=False,
                 new_columns=["q", "t", "fident", "alnlen", "qs", "qe", "ts", "te", "evalue", "bits", "qlen", "tlen"],
                 schema_overrides={"q": pl.String, "t": pl.String})
alt4 = {}
for spn, f in [("human", "UP000005640_9606.fasta.gz"), ("oyster", "UP000005408_29159.fasta.gz")]:
    seqs = fasta(R / f)
    A = set(seqs)
    n = len(A)
    xr = pl.read_csv(R / f"{spn}_interpro_xrefs.tsv.gz", separator="\t", infer_schema_length=0)
    cols = [c for c in xr.columns if c != "Entry"]
    xr = xr.filter(pl.col("Entry").is_in(list(A))).filter(
        pl.any_horizontal([pl.col(c).is_not_null() & (pl.col(c).str.strip_chars() != "") for c in cols]))
    ip = set(xr["Entry"].to_list())
    f1 = {a: [x for x in af.get(a, []) if x["entryId"] == f"AF-{a}-F1"] for a in A}
    afm = {a for a in A if f1[a]}
    iso_only = {a for a in A if af.get(a) and not f1[a]}
    f1_seq_differs = sum(1 for a in afm if f1[a][0].get("uniprotSequence") != seqs[a])
    m = mm.filter(pl.col("q").is_in(list(A)))
    hit = set(m.filter(pl.col("fident") > 0.30)["q"].to_list())
    hit_ns = set(m.filter((pl.col("fident") > 0.30) & (pl.col("t") != pl.col("q")))["q"].to_list())
    hit_ge = set(m.filter(pl.col("fident") >= 0.30)["q"].to_list())
    anyl = afm | ip | hit
    pct = lambda x: 100 * x / n  # noqa: E731
    V[f"{spn}_proteins"] = n
    V[f"{spn}_af_model"] = len(afm)
    V[f"{spn}_af_model_pct"] = pct(len(afm))
    V[f"{spn}_af_mean_plddt_median"] = statistics.median(f1[a][0]["globalMetricValue"] for a in afm)
    V[f"{spn}_af_isoform_model_only"] = len(iso_only)
    V[f"{spn}_af_isoform_model_only_pct"] = pct(len(iso_only))
    V[f"{spn}_interpro_match"] = len(ip)
    V[f"{spn}_interpro_match_pct"] = pct(len(ip))
    V[f"{spn}_sprot_hit"] = len(hit)
    V[f"{spn}_sprot_hit_pct"] = pct(len(hit))
    V[f"{spn}_sprot_hit_nonself"] = len(hit_ns)
    V[f"{spn}_sprot_hit_nonself_pct"] = pct(len(hit_ns))
    V[f"{spn}_any_label"] = len(anyl)
    V[f"{spn}_any_label_pct"] = pct(len(anyl))
    V[f"{spn}_no_label"] = n - len(anyl)
    V[f"{spn}_no_label_pct"] = pct(n - len(anyl))
    V[f"{spn}_no_interpro_no_sprot_pct"] = pct(len(A - ip - hit))
    V[f"{spn}_no_interpro_no_sprot_nonself_pct"] = pct(len(A - ip - hit_ns))
    V[f"{spn}_no_label_nonself_pct"] = pct(n - len(afm | ip | hit_ns))
    V[f"{spn}_sprot_hit_at_exactly_30pct"] = len(hit_ge - hit)
    canon = [g_["globalMetricValue"] for a in afm for g_ in af[a][:1]]
    canon_f1 = [x["globalMetricValue"] for a in afm for x in af[a] if x["entryId"] == f"AF-{a}-F1"]
    alt4[spn] = dict(
        no_interpro_no_sprot_nonself_pct=pct(len(A - ip - hit_ns)),
        sprot_hit_ge=len(hit_ge),
        af_isoform_only=sum(1 for a in afm if not any(x.get("uniprotSequence") == seqs[a] for x in af[a])),
        plddt_median_f1_only=statistics.median(canon_f1),
        f1_seq_differs=f1_seq_differs,
        any_label_old_af=len({a for a in A if af.get(a)} | ip | hit),
        no_label_nonself_count=n - len(afm | ip | hit_ns),
        exactly30_nonself=len(set(m.filter((pl.col("fident") >= 0.30) & (pl.col("t") != pl.col("q")))["q"].to_list()) - hit_ns),
        exactly30_rows_fident_eq=m.filter(pl.col("fident") == 0.30)["q"].n_unique(),
    )
    rows = [l.rstrip("\n").split("\t") for l in open(S / f"busco_{spn}/run_metazoa_odb12/full_table.tsv") if not l.startswith("#")]
    for status in ["Complete", "Duplicated", "Fragmented"]:
        prots = {r[2].split("|")[1] for r in rows if r[1] == status and len(r) > 2 and r[2]}
        V[f"{spn}_busco_{status.lower()}"] = len(prots)
print(json.dumps(alt4, indent=1))
print("raw OX counts (no merge):", raw_org, raw_pdb_org, raw_kb_tax)
print("alt top5 (denominator = entries with a species):", alt_sp_pdb_top5, alt_kb_sp_top5)

# ---------------------------------------------------------------- compare
num = pl.read_parquet(W / "numbers.parquet")
bad_meta = [r["id"] for r in num.iter_rows(named=True)
            if r["source_sha256"] != src[r["source"]]["sha256"] or r["source_file"] != src[r["source"]]["file"]]
print("rows whose source_sha256/source_file disagree with sources.json:", bad_meta)
print("duplicate ids:", num["id"].is_duplicated().sum())


def decimals(x):
    s = repr(float(x))
    if "e" in s:
        return 0
    frac = s.split(".")[1]
    return 0 if frac == "0" else len(frac)


REASON = {}
NOTE = {
    "sp_organisms": f"Distinct NCBI taxids after resolving merged ids: 14,903 resolvable + 2 ids deleted from NCBI (3470597 EBV strain, 147272 Paspalum notatum) kept as their own organisms = 14,905. Raw OX ids without merging: {raw_org}.",
    "sp_pdb_organisms": f"Same with or without resolving merged ids ({raw_pdb_org} raw).",
    "kb_pdb_taxid_ncbi": f"Merged ids resolved; {raw_kb_tax} raw OX ids. Includes 0 deleted ids (the deleted EBV id is a virus, so not cellular).",
    "sp_entries_no_species": "56 EBV strain 3470597 + 2 Paspalum notatum 147272 (both deleted from NCBI) + 1 'Corticium sp.' with no species rank. The 2 Paspalum entries drop out of the plant species count.",
    "human_af_model": f"Counts AF-<acc>-F1 models. For {alt4['human']['f1_seq_differs']} of them the F1 model's sequence differs from the current FASTA sequence, so 'of its own sequence' is not exact for those.",
    "oyster_af_model": f"Counts AF-<acc>-F1 models; for {alt4['oyster']['f1_seq_differs']} the model sequence differs from the FASTA.",
    "human_any_label": "Unchanged by the F1 rule: all 104 isoform-only proteins have another label.",
    "human_no_label_nonself_pct": f"{alt4['human']['no_label_nonself_count']} proteins.",
    "oyster_sprot_hit_at_exactly_30pct": f"Self-hits count. Excluding self-hits also gives {alt4['oyster']['exactly30_nonself']}. MMseqs2 rounds fident to 3 decimals, so 'exactly 0.300' means 0.2995 to 0.3005.",
    "human_sprot_hit_at_exactly_30pct": f"Excluding self-hits gives {alt4['human']['exactly30_nonself']}.",
    "larsen2017_estimated_species_Plants": "Matches the notebook; the paper is not on disk. The estimate (340,000) is below COL's described Plantae (389,184).",
    "sp_species_share_Archaea": "NCBI species counted against COL species; the two define Archaea and Protozoa differently, which inflates these two ratios.",
    "sp_species_share_Protozoa": "See sp_species_share_Archaea: COL Protozoa (1,433) is far narrower than 'all other NCBI Eukaryota'.",
}

rows = []
for r in num.iter_rows(named=True):
    i, stated = r["id"], r["value"]
    if not sha_ok.get(r["source"], False):
        rows.append((i, stated, None, "UNVERIFIABLE", "source file sha256 does not match sources.json"))
        continue
    if i not in V:
        rows.append((i, stated, None, "UNVERIFIABLE", "no recomputation written"))
        continue
    mine = float(V[i])
    d = decimals(stated)
    ok = round(mine, d) == round(stated, d) or abs(round(mine, d) - stated) < 10 ** (-d) / 2 + 1e-9
    verdict = "HOLDS" if ok else "WRONG"
    reason = REASON.get(i, "") if not ok else NOTE.get(i, "")
    if not ok and not reason:
        reason = f"recomputed {mine} vs stated {stated}"
    rows.append((i, stated, mine, verdict, reason))

out = pl.DataFrame(rows, schema=["id", "stated_value", "recomputed_value", "verdict", "reason"], orient="row")
out.write_parquet(HERE / "verdicts.parquet")
print(out["verdict"].value_counts())
with pl.Config(tbl_rows=200, fmt_str_lengths=60, tbl_width_chars=250):
    print(out.filter(pl.col("verdict") != "HOLDS"))
    print(out.select("id", "stated_value", "recomputed_value", "verdict"))
