"""P24856 against Swiss-Prot eight ways (scripts/why/antifreeze_masking.sbatch): which runs
find the protein's own entry, and with what bit score and E-value.

    python3 -I scripts/why/antifreeze.py --sherlock ~/data/why-numbers/2026_03/sherlock/out_antifreeze
Writes results/why/antifreeze_masking.parquet (one row per run) and the first entry of
results/why/failure_modes.yaml.
"""

import argparse
import datetime
import pathlib

import polars as pl
import yaml

OUT = pathlib.Path("results/why")
ACC = "P24856"
MMSEQS = {
    "mmseqs_as_job_47961379": "masking as job 47961379 ran it (--mask 1 --comp-bias-corr 1, the defaults)",
    "mmseqs_mask0": "low-complexity masking off (--mask 0)",
    "mmseqs_compbias0": "composition-bias correction off (--comp-bias-corr 0)",
    "mmseqs_mask0_compbias0": "both off (--mask 0 --comp-bias-corr 0)",
}
HMMER = {
    "phmmer": ("phmmer", "defaults, bias filter on"),
    "phmmer_nobias": ("phmmer", "bias filter off (--nobias)"),
    "jackhmmer": ("jackhmmer", "defaults, bias filter on"),
    "jackhmmer_nobias": ("jackhmmer", "bias filter off (--nobias)"),
}
COLS = ["query", "target", "fident", "alnlen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen"]


def acc(name):
    return name.split("|")[1] if "|" in name else name


def mmseqs_row(d, name):
    p = d / f"{name}.tsv"
    hits = pl.read_csv(p, separator="\t", has_header=False, new_columns=COLS) if p.stat().st_size else pl.DataFrame(schema={c: pl.Utf8 for c in COLS})
    hits = hits.with_columns(pl.col("target").map_elements(acc, return_dtype=pl.Utf8))
    own = hits.filter(pl.col("target") == ACC)
    return {
        "run": name, "tool": "MMseqs2", "version": (d / "mmseqs_version.txt").read_text().strip(), "setting": MMSEQS[name],
        "hits": hits.height, "self_hit": own.height > 0,
        "self_bits": float(own["bits"][0]) if own.height else None, "self_evalue": float(own["evalue"][0]) if own.height else None,
        "command": (d / f"{name}.cmd").read_text().strip(),
    }


def hmmer_row(d, name):
    """--tblout: one line per target sequence, the full-sequence E-value in column 5 and
    bit score in column 6. jackhmmer writes the hits of every round; the last round's are
    kept (the round number is not in the table, so the last line per target wins)."""
    tool, setting = HMMER[name]
    best = {}
    for line in (d / f"{name}.tblout").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        best[acc(f[0])] = (float(f[4]), float(f[5]))
    own = best.get(ACC)
    version = next(x for x in (d / f"{name}.out").read_text().splitlines() if x.startswith("# HMMER"))
    return {
        "run": name, "tool": tool, "version": version.lstrip("# ").strip(), "setting": setting,
        "hits": len(best), "self_hit": own is not None,
        "self_bits": own[1] if own else None, "self_evalue": own[0] if own else None,
        "command": (d / f"{name}.cmd").read_text().strip(),
    }


def finding(t):
    mm = t.filter(pl.col("tool") == "MMseqs2")
    found = mm.filter(pl.col("self_hit"))
    if found.height == 0:
        return ("MMseqs2 returns no hit for P24856 under any of four masking settings",
                "the k-mer prefilter on a sequence that is mostly the repeat AAT")
    first = found["setting"].to_list()
    base = mm.filter(pl.col("run") == "mmseqs_as_job_47961379")["self_hit"][0]
    if base:
        return ("MMseqs2 finds P24856's own entry when run alone, though job 47961379 did not", "next: compare the two runs' inputs")
    return (f"MMseqs2 finds P24856's own entry only with {'; '.join(first)}", None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sherlock", required=True, type=pathlib.Path, help="out_antifreeze/ from the job")
    args = ap.parse_args()
    d = args.sherlock
    rows = [mmseqs_row(d, n) for n in MMSEQS] + [hmmer_row(d, n) for n in HMMER]
    t = pl.DataFrame(rows, schema_overrides={"self_bits": pl.Float64, "self_evalue": pl.Float64})
    t.write_parquet(OUT / "antifreeze_masking.parquet")
    print(t.drop("command").to_pandas().to_string())
    text, next_step = finding(t)
    mm_found = t.filter((pl.col("tool") == "MMseqs2") & pl.col("self_hit")).height
    entry = {
        "tool": "MMseqs2", "version": t["version"][0], "finding": text,
        "number": f"{mm_found} of 4 MMseqs2 runs find the protein's own entry",
        "command": t.filter(pl.col("run") == "mmseqs_as_job_47961379")["command"][0],
        "source": "results/why/antifreeze_masking.parquet (scripts/why/antifreeze_masking.sbatch)",
        "date": datetime.date.today().isoformat(),
        "status": "open" if next_step else "explained",
    }
    if next_step:
        entry["next"] = next_step
    (OUT / "failure_modes.yaml").write_text(yaml.safe_dump([entry], sort_keys=False, allow_unicode=True, width=120))
    print(yaml.safe_dump([entry], sort_keys=False))


if __name__ == "__main__":
    main()
