#!/usr/bin/env python3
"""Write data/query_sets/<name>/ from a copy of the region benchmark's own inputs.

A query set is everything a tool author needs to run a tool and everything the scorer needs
to grade it, frozen so that every submission is graded against the same files:

    queries.fasta              the query proteins, as the benchmark searched them
    relatives.tsv              each relative proteome a search tool may search: its label,
                               QfO proteome, taxon, time since its line split from ours,
                               the file inside the QfO tarball and that file's sha256
    truth/pfam/                human_truth.parquet (Pfam-A domain instances on the queries)
                               and maps/<relative>.parquet (Pfam-A domains on each relative)
    truth/swissprot/           the same for Swiss-Prot features
    covariates.parquet         per-query values the strata are cut on (pLDDT, disorder)
    identity/<relative>.parquet  per truth instance: best identity to the relative
    SOURCE.json                where every file came from, with its sha256

--src is a local copy of the Sherlock files listed in SOURCE.json (the midi-plus run of
seanome/2024-kmerseek-analysis nextflow-runs/qfo-pfam-region-benchmark, scored 2026-09-11 to
2026-09-13). Nothing is recomputed here: files are copied byte for byte.
"""

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import polars as pl

SHERLOCK = (
    "/scratch/users/olgabot/2024-kmerseek-analysis/nextflow-runs/"
    "qfo-pfam-region-benchmark/data/midi-plus"
)
QFO_TARBALL = (
    "https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/"
    "qfo_release-2020_04_with_updated_UP000008143/"
    "QfO_release_2020_04_with_updated_UP000008143.tar.gz"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--src", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=Path("data/results.parquet"))
    ap.add_argument("--name", default="human_midi_plus")
    ap.add_argument("--out", type=Path, default=Path("data/query_sets"))
    a = ap.parse_args()
    out = a.out / a.name
    out.mkdir(parents=True, exist_ok=True)

    # The relatives the leaderboard already scores this query set against, with the time
    # since each split from the human line exactly as results.parquet carries it.
    res = pl.read_parquet(a.results).filter(
        (pl.col("experiment") == "qfo_region_midi_plus")
        & pl.col("target_db").str.starts_with("QfO ")
    )
    mya = {
        r["target_db"].removeprefix("QfO "): r["target_mya_to_human"]
        for r in res.select("target_db", "target_mya_to_human").unique().iter_rows(named=True)
    }
    species = pl.read_csv(a.src / "qfo_species.tsv", separator="\t", infer_schema=False)
    shas = {}
    for line in (a.src / "qfo_sha256.txt").read_text().splitlines():
        digest, rel = line.split()
        shas[rel] = digest
    rows = []
    for r in species.filter(pl.col("label").is_in(list(mya))).iter_rows(named=True):
        member = f"{r['subdir']}/{r['proteome']}_{r['taxon']}.fasta"
        rows.append(
            {
                "relative": r["label"],
                "scientific_name": r["scientific_name"],
                "taxon": int(r["taxon"]),
                "qfo_proteome": r["proteome"],
                "n_proteins": int(r["n_proteins"]),
                "mya_to_human": mya[r["label"]],
                "file_in_qfo_tarball": member,
                "sha256": shas[member],
            }
        )
    if len(rows) != len(mya):
        raise SystemExit(f"relatives in results.parquet without a QfO row: {set(mya) - {r['relative'] for r in rows}}")
    # E. coli has no TimeTree age (null): it sorts last.
    rel = pl.DataFrame(rows).sort("mya_to_human", "relative", nulls_last=True)
    rel.write_csv(out / "relatives.tsv", separator="\t")

    copies = {
        "queries.fasta": "qfo/Eukaryota/UP000005640_9606.fasta",
        "covariates.parquet": "results/truth/human_query_covariates.parquet",
        "truth/pfam/human_truth.parquet": "results/truth/human_domain_truth.parquet",
        "truth/swissprot/human_truth.parquet": "results/truth_swissprot/human_swissprot_truth.parquet",
        "query_set_summary.json": "midi_testset_summary.json",
    }
    local = {
        "qfo/Eukaryota/UP000005640_9606.fasta": "UP000005640_9606.fasta",
        "results/truth/human_query_covariates.parquet": "truth/human_query_covariates.parquet",
        "results/truth/human_domain_truth.parquet": "truth/human_domain_truth.parquet",
        "results/truth_swissprot/human_swissprot_truth.parquet": "truth_swissprot/human_swissprot_truth.parquet",
        "midi_testset_summary.json": "midi_testset_summary.json",
    }
    for s in rel["relative"]:
        copies[f"truth/pfam/maps/{s}.parquet"] = f"results/truth/{s}_domain_map.parquet"
        local[f"results/truth/{s}_domain_map.parquet"] = f"truth/{s}_domain_map.parquet"
        copies[f"truth/swissprot/maps/{s}.parquet"] = f"results/truth_swissprot/{s}_domain_map.parquet"
        local[f"results/truth_swissprot/{s}_domain_map.parquet"] = f"truth_swissprot/{s}_domain_map.parquet"
        copies[f"identity/{s}.parquet"] = f"results/identity/{s}.domain_identity.parquet"
        local[f"results/identity/{s}.domain_identity.parquet"] = f"identity/{s}.domain_identity.parquet"
    files = {}
    for dst, sherlock_rel in copies.items():
        src = a.src / local[sherlock_rel]
        (out / dst).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, out / dst)
        files[dst] = {"from": f"{SHERLOCK}/{sherlock_rel}", "sha256": sha256(out / dst)}

    n_queries = sum(1 for line in (out / "queries.fasta").open() if line.startswith(">"))
    source = {
        "name": a.name,
        "description": "human proteins from chromosome 6 plus the MHC partners on other chromosomes "
        "(B2M, CD1A-E, MR1, KIR and LILR genes), a test subset of the proteome",
        "n_queries": n_queries,
        "experiment_on_the_leaderboard": "qfo_region_midi_plus",
        "written_by": "seanome/2024-kmerseek-analysis commit a47ea25 (HEAD of the Sherlock clone "
        "when the midi-plus run started)",
        "truth_versions": {"pfam": "Pfam 38.2", "swissprot": "uniprot_sprot.dat (release not recorded)"},
        "relatives_from": QFO_TARBALL,
        "scoring": {"min_overlap": 0.5, "strict_iou": 0.8, "dedup_transfer_modes": "off,on"},
        "files": files,
    }
    (out / "SOURCE.json").write_text(json.dumps(source, indent=1) + "\n")
    print(f"{out}: {n_queries} queries, {rel.height} relatives, {len(files)} files")


if __name__ == "__main__":
    main()
