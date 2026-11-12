"""The proteins kmerseek searches, in chunks of 2_000 like the dark-set pipeline.

Sets (column `query_set` in queries.tsv):
- four: the four proteins of step 6.
- oyster_dark: oyster proteins with neither an InterPro match nor a Swiss-Prot hit above
  30% identity (proteome_labels.parquet).
- oyster_lit_control: the same number of oyster proteins drawn at random (seed 2026) from
  those with an InterPro match or a Swiss-Prot hit. A kmerseek hit on a dark protein cannot
  be checked; on a lit one it can, so this set says how often kmerseek's hits are right.
- human_dark: human proteins with neither label.

    python3 -I scripts/why/kmerseek_queries.py --raw ~/data/why-numbers/2026_03
Writes ~/data/why-numbers/2026_03/kmerseek_inputs/queries/chunk_NNNN.fasta and
results/why/kmerseek/queries.tsv.
"""

import argparse
import gzip
import json
import pathlib
import sys

import polars as pl

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from four import FOUR  # noqa: E402
from labels import PROTEOMES, read_fasta  # noqa: E402

CHUNK = 2000
SEED = 2026


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=pathlib.Path)
    args = ap.parse_args()
    t = pl.read_parquet("results/why/proteome_labels.parquet")
    function_label = pl.col("interpro_match") | pl.col("sprot_hit")
    oyster_dark = t.filter(pl.col("species") == "oyster", ~function_label)["accession"].sort().to_list()
    oyster_lit = t.filter(pl.col("species") == "oyster", function_label).sort("accession").sample(
        n=len(oyster_dark), seed=SEED, shuffle=False)["accession"].sort().to_list()
    human_dark = t.filter(pl.col("species") == "human", ~function_label)["accession"].sort().to_list()
    seqs = {}
    for tag, stem in PROTEOMES.items():
        seqs.update(read_fasta(args.raw / f"{stem}.fasta.gz"))
    for acc in FOUR:
        seqs[acc] = json.loads((args.raw / f"uniprot_{acc}.json").read_text())["sequence"]["value"]
    rows = [(a, "four") for a in FOUR] + [(a, "oyster_dark") for a in oyster_dark] + \
           [(a, "oyster_lit_control") for a in oyster_lit] + [(a, "human_dark") for a in human_dark]
    accs = [a for a, _ in rows]
    if len(set(accs)) != len(accs):
        raise ValueError("an accession is in two query sets")
    out_dir = args.raw / "kmerseek_inputs" / "queries"
    out_dir.mkdir(parents=True, exist_ok=True)
    for f in out_dir.glob("chunk_*.fasta"):
        f.unlink()
    table = []
    for i in range(0, len(rows), CHUNK):
        name = f"chunk_{i // CHUNK:04d}"
        with open(out_dir / f"{name}.fasta", "w") as fh:
            for acc, qset in rows[i : i + CHUNK]:
                fh.write(f">{acc}\n{seqs[acc]}\n")
                table.append({"accession": acc, "query_set": qset, "chunk": name, "length": len(seqs[acc])})
    pathlib.Path("results/why/kmerseek").mkdir(parents=True, exist_ok=True)
    q = pl.DataFrame(table)
    q.write_csv("results/why/kmerseek/queries.tsv", separator="\t")
    print(q.group_by("query_set").agg(n=pl.len(), residues=pl.col("length").sum()).sort("query_set"))
    print(q.group_by("chunk").agg(n=pl.len()).sort("chunk"))


if __name__ == "__main__":
    main()
