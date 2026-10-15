"""Step 4: one row per protein of the human and oyster one-protein-per-gene sets.

Columns: length; AlphaFold DB model (yes/no, mean pLDDT); any InterPro match (yes/no);
a Swiss-Prot hit above 30% identity by MMseqs2 (yes/no); BUSCO status. A protein is
labelled when any of the three yes/no columns is yes.

    python3 -I scripts/why/labels.py --raw ~/data/why-numbers/2026_03 --sherlock ~/data/why-numbers/2026_03/sherlock
"""

import argparse
import gzip
import json
import pathlib
import sys

import polars as pl

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from numbers_table import Numbers  # noqa: E402

OUT = pathlib.Path("results/why")
PROTEOMES = {"human": "UP000005640_9606", "oyster": "UP000005408_29159"}
MIN_FIDENT = 0.30  # "above 30% identity": MMseqs2 fident > 0.30, E <= 1e-3 (search setting)


def read_fasta(path):
    seqs, acc, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc, buf = line.split("|")[1], []
            else:
                buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def afdb_table(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            r = json.loads(line)
            models = r.get("models") or []
            # A whole-protein model (F1 starting at residue 1) when there is one.
            m = next((m for m in models if (m.get("entryId") or "").endswith("-F1")), models[0] if models else None)
            rows.append(
                {
                    "accession": r["accession"],
                    "af_model": bool(models),
                    "af_mean_plddt": m["globalMetricValue"] if m else None,
                    "af_version": m["latestVersion"] if m else None,
                    "af_model_sequence": m["sequence"] if m else None,
                }
            )
    return pl.DataFrame(rows, schema_overrides={"af_mean_plddt": pl.Float64, "af_version": pl.Int64})


def interpro_set(path):
    """Accessions in an InterPro API proteome listing (oyster only; the human one timed out)."""
    with gzip.open(path, "rt") as fh:
        return {json.loads(line)["accession"] for line in fh}


def interpro_xref_set(path):
    """Accessions UniProt cross-references to any InterPro entry or member-database signature."""
    x = pl.read_csv(path, separator="\t", infer_schema_length=0)
    members = [c for c in x.columns if c != "Entry"]
    return set(x.filter(pl.any_horizontal(pl.col(members).is_not_null() & (pl.col(members) != "")))["Entry"])


def acc_of(col):
    """MMseqs2 writes the accession alone ("P42212"); a full "sp|P42212|GFP_AEQVI" header is
    reduced to its accession too."""
    c = pl.col(col)
    return pl.when(c.str.contains("|", literal=True)).then(c.str.split("|").list.get(1)).otherwise(c)


def mmseqs_table(path):
    cols = ["query", "target", "fident", "alnlen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen"]
    hits = pl.read_csv(path, separator="\t", has_header=False, new_columns=cols)
    hits = hits.with_columns(
        query=acc_of("query"),
        target=acc_of("target"),
    ).filter(pl.col("fident") > MIN_FIDENT)
    hits = hits.with_columns(self_hit=pl.col("query") == pl.col("target"))
    best = hits.sort(["query", "bits"], descending=[False, True]).group_by("query", maintain_order=True).first()
    best_ns = hits.filter(~pl.col("self_hit")).sort(["query", "bits"], descending=[False, True]).group_by("query", maintain_order=True).first()
    return hits, best.select(
        accession="query", sprot_best_target="target", sprot_best_fident="fident", sprot_best_evalue="evalue"
    ).join(
        best_ns.select(accession="query", sprot_best_nonself_target="target", sprot_best_nonself_fident="fident"),
        on="accession", how="left",
    )


def busco_table(path):
    """full_table.tsv: one row per BUSCO gene and the protein(s) it was found in."""
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) > 2 and f[2]:
            rows.append({"busco_id": f[0], "busco_status": f[1], "protein": f[2]})
    df = pl.DataFrame(rows)
    df = df.with_columns(accession=pl.col("protein").str.split("|").list.get(1).fill_null(pl.col("protein")))
    # A protein can hold more than one BUSCO gene: keep the best status (Complete before
    # Duplicated before Fragmented) and how many BUSCO genes point at it.
    order = {"Complete": 0, "Duplicated": 1, "Fragmented": 2}
    df = df.with_columns(rank=pl.col("busco_status").replace_strict(order, default=9))
    return df.sort("rank").group_by("accession").agg(
        busco_status=pl.col("busco_status").first(), busco_genes=pl.len()
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=pathlib.Path)
    ap.add_argument("--sherlock", required=True, type=pathlib.Path, help="the job's out/ folder, copied back")
    args = ap.parse_args()
    raw, sh = args.raw, args.sherlock

    af = afdb_table(raw / "afdb_predictions.jsonl.gz")
    _, best = mmseqs_table(sh / "mmseqs_vs_sprot.tsv.gz")
    N = Numbers()
    frames = []
    for tag, stem in PROTEOMES.items():
        seqs = read_fasta(raw / f"{stem}.fasta.gz")
        fields = pl.read_csv(raw / f"{tag}_proteome_fields.tsv.gz", separator="\t", infer_schema_length=0).select(
            accession="Entry", reviewed="Reviewed", gene="Gene Names (primary)", protein_name="Protein names", interpro_xref="InterPro",
        )
        ipr = interpro_xref_set(raw / f"{tag}_interpro_xrefs.tsv.gz")
        busco = busco_table(next((sh / f"busco_{tag}").glob("run_*/full_table.tsv")))
        df = (
            pl.DataFrame({"accession": list(seqs), "sequence": list(seqs.values())})
            .with_columns(species=pl.lit(tag), length=pl.col("sequence").str.len_chars())
            .join(fields, on="accession", how="left")
            .join(af, on="accession", how="left")
            .join(best, on="accession", how="left")
            .join(busco, on="accession", how="left")
            .with_columns(
                interpro_match=pl.col("accession").is_in(list(ipr)),
                af_model=pl.col("af_model").fill_null(False),
                af_model_matches_sequence=pl.col("af_model_sequence") == pl.col("sequence"),
                sprot_hit=pl.col("sprot_best_target").is_not_null(),
                sprot_hit_nonself=pl.col("sprot_best_nonself_target").is_not_null(),
                busco_status=pl.col("busco_status").fill_null("not a BUSCO gene"),
                busco_genes=pl.col("busco_genes").fill_null(0),
            )
            .with_columns(any_label=pl.col("af_model") | pl.col("interpro_match") | pl.col("sprot_hit"))
            .drop("af_model_sequence")
        )
        if df["reviewed"].null_count():
            raise ValueError(f"{tag}: {df['reviewed'].null_count()} canonical accessions missing from the proteome fields")
        frames.append(df)

        n = df.height
        src = f"{tag}_canonical_fasta"
        rel = "UniProt 2026_03"
        N.add(f"{tag}_proteins", 4, f"{tag}: proteins in the one-protein-per-gene reference set", n, "proteins", src, rel)
        for col, what, source in (
            ("af_model", "with an AlphaFold DB model", "afdb_predictions"),
            ("interpro_match", "with any InterPro match (InterPro entry or member-database signature, as UniProt cross-references them)", f"{tag}_interpro_xrefs"),
            ("sprot_hit", "with a Swiss-Prot hit above 30% identity (MMseqs2, self-hits count)", src),
            ("sprot_hit_nonself", "with a Swiss-Prot hit above 30% identity other than the protein itself", src),
            ("any_label", "with at least one of the three labels", src),
        ):
            k = int(df[col].sum())
            N.add(f"{tag}_{col}", 4, f"{tag}: proteins {what}", k, "proteins", source, rel)
            N.add(f"{tag}_{col}_pct", 4, f"{tag}: share of proteins {what}", round(100 * k / n, 2), "% of proteins", source, rel)
        k = int((~df["any_label"]).sum())
        N.add(f"{tag}_no_label", 4, f"{tag}: proteins with none of the three labels", k, "proteins", src, rel)
        N.add(f"{tag}_no_label_pct", 4, f"{tag}: share of proteins with none of the three labels", round(100 * k / n, 2), "% of proteins", src, rel)
        if tag == "oyster":
            api = interpro_set(raw / "interpro_oyster_proteins_with_match.jsonl.gz")
            accs = set(df["accession"])
            both, x_only, api_only = len(ipr & api & accs), len((ipr - api) & accs), len((api - ipr) & accs)
            N.add("oyster_interpro_api_and_xref", 4, "oyster: proteins with a match by both the InterPro API listing and UniProt's cross-references", both, "proteins", "interpro_oyster", "InterPro API")
            N.add("oyster_interpro_xref_only", 4, "oyster: proteins with a match by UniProt's cross-references only", x_only, "proteins", "interpro_oyster", "InterPro API")
            N.add("oyster_interpro_api_only", 4, "oyster: proteins with a match in the InterPro API listing only", api_only, "proteins", "interpro_oyster", "InterPro API")
        fn = df.filter(~pl.col("interpro_match") & ~pl.col("sprot_hit")).height
        N.add(f"{tag}_no_interpro_no_sprot_pct", 4, f"{tag}: share of proteins with neither an InterPro match nor a Swiss-Prot hit", round(100 * fn / n, 2), "% of proteins", src, rel,
              note="the two labels that name a function; AlphaFold DB models nearly every UniProt sequence")
        N.add(f"{tag}_af_mean_plddt_median", 4, f"{tag}: median of per-protein mean pLDDT, proteins with a model", round(df["af_mean_plddt"].median(), 2), "pLDDT", "afdb_predictions", "AlphaFold DB API")
        for status in ("Complete", "Duplicated", "Fragmented"):
            k = df.filter(pl.col("busco_status") == status).height
            N.add(f"{tag}_busco_{status.lower()}", 4, f"{tag}: proteins BUSCO metazoa_odb12 calls {status}", k, "proteins", src, rel)

    table = pl.concat(frames)
    table.drop("sequence").write_parquet(OUT / "proteome_labels.parquet")
    N.write(OUT / "numbers.parquet", replace_steps={4})
    print(N.frame().select("id", "value", "unit").to_pandas().to_string())


if __name__ == "__main__":
    main()
