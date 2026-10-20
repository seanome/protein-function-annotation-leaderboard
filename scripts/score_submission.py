#!/usr/bin/env python3
"""Score one tool submission with the leaderboard's own scorer and write the rows it adds.

A submission is a folder submissions/<name>/ holding

    predictions.tsv or predictions.tsv.gz   one row per call (CONTRIBUTING.md, step 2)
    timing.tsv                              what each run cost on the author's hardware
    tool.yaml                               the tool card: what it compares, what it needs,
                                            licence, command, container

Steps: check every file and every row, write each (setting, relative) as the table the
scorer reads, run scorer/evaluate_domain_calls.py with the flags the leaderboard used
(scorer/SOURCE.json), turn its metrics into rows in the layout of data/results.parquet with
the leaderboard's own code (scorer/build_leaderboard_long.py), and write

    <out>/results.parquet   new leaderboard rows, experiment = "submission"
    <out>/compute.parquet   the timing log as rows, provenance = "reported", hardware named
    <out>/card.json         the tool card built from tool.yaml
    <out>/comment.md        what the workflow posts on the pull request

A problem with the submission is reported in comment.md and the script exits 1, so the
pull request shows what to fix rather than a stack trace.
"""

import argparse
import gzip
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scorer"))
import build_leaderboard_long as long_tables  # noqa: E402

PRED_COLUMNS = [
    "query_set", "relative", "query", "call_start", "call_end", "target",
    "target_start", "target_end", "label", "score", "tool", "version", "setting",
]
TIMING_COLUMNS = [
    "query_set", "relative", "setting", "step", "wall_s", "cores", "peak_ram_gb",
    "gpu_model", "hardware",
]
CARD_KEYS = {
    "tool": "short name in lower case, letters, digits and _ only",
    "name": "name as people write it",
    "version": "version that made the predictions",
    "kind": "letters, profiles, shapes, kmers or models",
    "compares": "one or two sentences: what the tool compares",
    "needs": "sequence_only, structure_prediction and gpu, each true or false",
    "licence": "licence name",
    "licence_url": "link to the licence",
    "command": "the command that made the predictions",
    "container": "image reference with its digest",
}
KINDS = {"letters", "profiles", "shapes", "kmers", "models"}
STEPS = {"search", "database_build", "query_conversion"}
TRUTH = {
    "pfam": ("Pfam-A domain instances", "Pfam 38.2"),
    "swissprot": ("Swiss-Prot features", None),
    "contributed": ("Wet-lab rows contributed through the site", None),
}
SWISSPROT_TYPES = {
    "ACT_SITE", "BINDING", "COILED", "DNA_BIND", "DOMAIN", "INTRAMEM", "MOTIF",
    "REGION", "REPEAT", "SITE", "TRANSMEM", "ZN_FING",
}
MAX_PROBLEMS = 25


class Problems(list):
    def add(self, msg: str) -> None:
        if len(self) < MAX_PROBLEMS:
            self.append(msg)


def read_fasta_lengths(path: Path) -> dict[str, int]:
    """accession -> length, accession taken the way the scorer takes it (sp|ACC|NAME)."""
    lengths, acc, n = {}, None, 0
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if acc:
                lengths[acc] = n
            name = line[1:].split()[0]
            acc = name.split("|")[1] if "|" in name else name
            n = 0
        else:
            n += len(line.strip())
    if acc:
        lengths[acc] = n
    return lengths


def accession(col: pl.Expr) -> pl.Expr:
    return (
        pl.when(col.str.contains(r"\|"))
        .then(col.str.split("|").list.get(1, null_on_oob=True))
        .otherwise(col)
    )


def read_predictions(sub: Path, problems: Problems) -> pl.DataFrame | None:
    files = [p for p in (sub / "predictions.tsv", sub / "predictions.tsv.gz") if p.exists()]
    if len(files) != 1:
        problems.add("The folder needs one of predictions.tsv or predictions.tsv.gz, not both.")
        return None
    path = files[0]
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if header != PRED_COLUMNS:
        problems.add(
            f"{path.name}: the header must be these {len(PRED_COLUMNS)} tab-separated "
            f"columns, in this order: {', '.join(PRED_COLUMNS)}. Found: {', '.join(header)}."
        )
        return None
    df = pl.read_csv(path, separator="\t", infer_schema=False, missing_utf8_is_empty_string=True)
    return df.with_columns(pl.int_range(2, pl.len() + 2).alias("line"))


def check_predictions(df: pl.DataFrame, qs_dir: Path, card: dict, problems: Problems) -> pl.DataFrame:
    """Parse numbers and check every row. Returns the rows with typed columns."""
    df = df.with_columns(
        *[pl.col(c).cast(pl.Int64, strict=False).alias(f"{c}_i") for c in
          ("call_start", "call_end", "target_start", "target_end")],
        pl.col("score").cast(pl.Float64, strict=False).alias("score_f"),
        accession(pl.col("query")).alias("query_acc"),
        pl.col("label").str.replace(r"^(PF\d{5})\.\d+$", "$1").alias("label"),
    )

    def first_lines(mask: pl.Expr, what: str) -> None:
        bad = df.filter(mask)
        if bad.height:
            lines = ", ".join(str(x) for x in bad["line"].head(5).to_list())
            problems.add(f"Rows that {what}: {bad.height:,} (first lines of the file: {lines}).")

    query_sets = {p.name for p in qs_dir.iterdir() if p.is_dir()}
    first_lines(~pl.col("query_set").is_in(sorted(query_sets)),
                f"name a query set that is not in data/query_sets/ ({', '.join(sorted(query_sets))})")
    for qs in sorted(set(df["query_set"].unique().to_list()) & query_sets):
        lengths = read_fasta_lengths(qs_dir / qs / "queries.fasta")
        rel = set(pl.read_csv(qs_dir / qs / "relatives.tsv", separator="\t")["relative"].to_list())
        here = pl.col("query_set") == qs
        first_lines(here & ~pl.col("relative").is_in(sorted(rel | {"none"})),
                    f"name a relative that is not in {qs}/relatives.tsv and is not 'none'")
        first_lines(here & ~pl.col("query_acc").is_in(list(lengths)),
                    f"name a query protein that is not in {qs}/queries.fasta")
        length = pl.col("query_acc").replace_strict(lengths, default=None, return_dtype=pl.Int64)
        first_lines(here & pl.col("query_acc").is_in(list(lengths))
                    & (pl.col("call_end_i") > length),
                    "end past the end of the query protein")
    first_lines(pl.col("call_start_i").is_null() | pl.col("call_end_i").is_null(),
                "have a call_start or call_end that is not a whole number")
    first_lines(pl.col("call_start_i") < 1, "start before residue 1 (positions count from 1)")
    first_lines(pl.col("call_end_i") < pl.col("call_start_i"), "end before they start")
    first_lines(pl.col("score_f").is_null() | pl.col("score_f").is_nan()
                | pl.col("score_f").is_infinite(), "have a score that is not a finite number")
    search = pl.col("relative") != "none"
    first_lines(search & ((pl.col("target") == "") | pl.col("target_start_i").is_null()
                          | pl.col("target_end_i").is_null()),
                "name a relative but leave target, target_start or target_end empty")
    first_lines(search & (pl.col("target_end_i") < pl.col("target_start_i")),
                "have a target interval that ends before it starts")
    first_lines(search & (pl.col("target_start_i") < 1),
                "have a target_start before residue 1")
    first_lines(search & (pl.col("label") != ""),
                "name a relative and also a label; a search call carries no label")
    first_lines(~search & (pl.col("label") == ""),
                "have relative 'none' and no label; a call with no relative must name its family or feature")
    label_ok = pl.col("label").str.contains(r"^PF\d{5}$") | pl.col("label").is_in(sorted(SWISSPROT_TYPES))
    first_lines(~search & ~label_ok,
                "carry a label that is neither a Pfam accession (PF00001) nor a Swiss-Prot feature type")
    first_lines(pl.col("tool") != str(card.get("tool")), f"name a tool other than tool.yaml's '{card.get('tool')}'")
    first_lines(pl.col("version") != str(card.get("version")),
                f"name a version other than tool.yaml's '{card.get('version')}'")
    first_lines(pl.col("setting") == "", "leave setting empty")
    first_lines(~pl.col("setting").str.contains(r"^[A-Za-z0-9_.+-]+$"),
                "have a setting with characters other than letters, digits and _ . + -")
    return df


def check_card(sub: Path, problems: Problems) -> dict:
    path = sub / "tool.yaml"
    if not path.exists():
        problems.add("tool.yaml is missing.")
        return {}
    try:
        card = yaml.safe_load(path.read_text()) or {}
    except yaml.YAMLError as e:
        problems.add(f"tool.yaml is not valid YAML: {e}")
        return {}
    for key, what in CARD_KEYS.items():
        if card.get(key) in (None, "", {}):
            problems.add(f"tool.yaml needs `{key}`: {what}.")
    if card.get("tool") and not re.fullmatch(r"[a-z0-9_]+", str(card["tool"])):
        problems.add("tool.yaml `tool` may hold only lower-case letters, digits and _.")
    if card.get("kind") and card["kind"] not in KINDS:
        problems.add(f"tool.yaml `kind` must be one of {', '.join(sorted(KINDS))}.")
    needs = card.get("needs") or {}
    for k in ("sequence_only", "structure_prediction", "gpu"):
        if not isinstance(needs.get(k), bool):
            problems.add(f"tool.yaml `needs.{k}` must be true or false.")
    if card.get("container") and "@sha256:" not in str(card["container"]):
        problems.add("tool.yaml `container` must name the image by digest (…@sha256:…), "
                     "so the exact image can be pulled again.")
    return card


def check_timing(sub: Path, preds: pl.DataFrame, problems: Problems) -> pl.DataFrame | None:
    path = sub / "timing.tsv"
    if not path.exists():
        problems.add("timing.tsv is missing.")
        return None
    t = pl.read_csv(path, separator="\t", infer_schema=False, missing_utf8_is_empty_string=True)
    if t.columns != TIMING_COLUMNS:
        problems.add(f"timing.tsv: the header must be {', '.join(TIMING_COLUMNS)}. Found: {', '.join(t.columns)}.")
        return None
    t = t.with_columns(
        pl.col("wall_s").cast(pl.Float64, strict=False),
        pl.col("cores").cast(pl.Int64, strict=False),
        pl.col("peak_ram_gb").cast(pl.Float64, strict=False),
    )
    if t.filter(pl.col("wall_s").is_null() | (pl.col("wall_s") <= 0)).height:
        problems.add("timing.tsv: wall_s must be a number of seconds above 0 on every row.")
    if t.filter(pl.col("cores").is_null() | (pl.col("cores") < 1)).height:
        problems.add("timing.tsv: cores must be a whole number of 1 or more on every row.")
    if t.filter(pl.col("peak_ram_gb").is_null() | (pl.col("peak_ram_gb") <= 0)).height:
        problems.add("timing.tsv: peak_ram_gb must be a number above 0 on every row.")
    if t.filter(pl.col("hardware").str.strip_chars() == "").height:
        problems.add("timing.tsv: hardware must name the machine (CPU model, and GPU model if one was used).")
    if t.filter(pl.col("gpu_model").str.strip_chars() == "").height:
        problems.add("timing.tsv: gpu_model must name the GPU, or say none.")
    if t.filter(~pl.col("step").is_in(sorted(STEPS))).height:
        problems.add(f"timing.tsv: step must be one of {', '.join(sorted(STEPS))}.")
    runs = preds.select("query_set", "relative", "setting").unique()
    timed = t.filter(pl.col("step") == "search").select("query_set", "relative", "setting").unique()
    missing = runs.join(timed, on=["query_set", "relative", "setting"], how="anti")
    if missing.height:
        first = missing.head(3).rows()
        problems.add(f"timing.tsv has no `search` row for {missing.height} (query set, relative, "
                     f"setting) in the predictions, for example {first}.")
    return t


def run_scorer(work: Path, qs_dir: Path, preds: pl.DataFrame, card: dict,
               contributed: Path | None) -> list[Path]:
    """One scorer call per (query set, relative, truth set); returns the metrics files."""
    metrics = []
    rel_tab = {
        qs: pl.read_csv(qs_dir / qs / "relatives.tsv", separator="\t")
        for qs in preds["query_set"].unique().to_list()
    }
    preds = preds.with_columns(
        pl.when(pl.col("call_end_i") == pl.col("call_start_i"))
        .then(pl.col("call_end_i") + 1).otherwise(pl.col("call_end_i")).alias("call_end"))
    for (qs, relative), g in preds.group_by(["query_set", "relative"], maintain_order=True):
        qsd = qs_dir / qs
        species = "all" if relative == "none" else relative
        manifest, direct_files = [], {}
        for (setting,), gs in g.group_by(["setting"], maintain_order=True):
            if relative == "none":
                # A call that names its family is graded only against the truth set whose
                # names it uses: Pfam accessions against Pfam, feature types against Swiss-Prot.
                for ts, ok in (("pfam", pl.col("label").str.contains(r"^PF\d{5}$")),
                               ("swissprot", pl.col("label").is_in(sorted(SWISSPROT_TYPES))),
                               ("contributed", pl.col("label").is_in(sorted(SWISSPROT_TYPES)))):
                    part = gs.filter(ok)
                    if part.height == 0:
                        continue
                    f = work / f"{qs}.{relative}.{setting}.{ts}.direct.tsv"
                    part.select("query", "label", "call_start", "call_end", "score",
                                pl.lit("NA").alias("evalue")).write_csv(f, separator="\t", include_header=False)
                    direct_files.setdefault(ts, {})[setting] = f
            else:
                f = work / f"{qs}.{relative}.{setting}.tsv"
                gs.select("query", "target", "call_start", "call_end", "target_start",
                          "target_end", "score", pl.lit("NA").alias("evalue")).write_csv(
                    f, separator="\t", include_header=False)
                manifest.append(f"{card['tool']}\t{setting}\t{f.resolve()}")
        truths = [("pfam", qsd / "truth/pfam/human_truth.parquet", qsd / f"truth/pfam/maps/{relative}.parquet"),
                  ("swissprot", qsd / "truth/swissprot/human_truth.parquet",
                   qsd / f"truth/swissprot/maps/{relative}.parquet")]
        if contributed is not None and (contributed / f"{qs}.parquet").exists():
            truths.append(("contributed", contributed / f"{qs}.parquet",
                           qsd / f"truth/swissprot/maps/{relative}.parquet"))
        mya = None
        if relative != "none":
            row = rel_tab[qs].filter(pl.col("relative") == relative)
            mya = row["mya_to_human"][0]
        for ts, truth, dmap in truths:
            out_dir = work / f"{qs}.{relative}.{ts}"
            out_dir.mkdir(parents=True, exist_ok=True)
            base = [sys.executable, str(ROOT / "scorer/evaluate_domain_calls.py"),
                    "--species", species, "--truth", str(truth.resolve()),
                    "--covariates", str((qsd / "covariates.parquet").resolve()),
                    "--min-overlap", "0.5", "--strict-iou", "0.8", "--truth-set", ts]
            if relative == "none":
                for setting, f in direct_files.get(ts, {}).items():
                    stem = f"{ts}.{card['tool']}.{setting}.all"
                    cmd = base + ["--regions", str(f.resolve()), "--tool", card["tool"],
                                  "--variant", setting, "--direct-annotation",
                                  "--calls-out", f"{stem}.calls.parquet",
                                  "--metrics-out", f"{stem}.metrics.parquet",
                                  "--curve-out", f"{stem}.curve.parquet"]
                    subprocess.run(cmd, cwd=out_dir, check=True, stdout=subprocess.DEVNULL)
            else:
                (out_dir / "manifest.tsv").write_text("\n".join(manifest) + "\n")
                cmd = base + ["--manifest", "manifest.tsv", "--domain-map", str(dmap.resolve()),
                              "--identity", str((qsd / f"identity/{relative}.parquet").resolve()),
                              "--dedup-transfer-modes", "off,on"]
                # E. coli has no TimeTree age: leave the option out rather than pass "None".
                cmd += ["--species-mya", str(mya)] if mya is not None else []
                subprocess.run(cmd, cwd=out_dir, check=True, stdout=subprocess.DEVNULL)
            metrics += sorted(out_dir.glob("*.metrics.parquet"))
    return metrics


def to_rows(metrics: list[Path], qs_dir: Path, card: dict, source: dict) -> pl.DataFrame:
    parts = []
    for m in metrics:
        qs = m.parent.name.split(".")[0]
        meta = json.loads((qs_dir / qs / "SOURCE.json").read_text())
        df = long_tables.qfo_metrics(
            m, m.name, "qfo_region_midi_plus", meta["description"],
            {card["tool"]: f"{card['name']} {card['version']}"}, truth=TRUTH,
        )
        if m.parent.name.split(".")[1] == "none":
            df = df.with_columns(pl.lit("none: the calls name the family").alias("target_db"),
                                 pl.lit(None, dtype=pl.String).alias("target_db_version"))
        parts.append(df.with_columns(pl.lit(qs).alias("query_set_id")))
    rows = pl.concat(parts, how="vertical")
    return rows.with_columns(
        pl.lit("submission").alias("experiment"),
        pl.lit(source["path"]).alias("source_path"),
        pl.lit(source["commit"]).cast(pl.String).alias("source_commit"),
    ).sort("query_set_id", "target_db", "truth_source", "dedup_or_raw", "setting",
           "stratum_axis", "stratum_bin", "metric")


def to_compute(t: pl.DataFrame, card: dict, n_queries: dict, source: dict) -> pl.DataFrame:
    rows = []
    for r in t.iter_rows(named=True):
        gpu = r["gpu_model"].strip()
        rows.append({
            "tool": card["tool"], "setting": r["setting"], "step": r["step"],
            "target": f"{r['query_set']}|{r['relative']}", "complete": None,
            "wall_s": r["wall_s"], "cores": r["cores"], "peak_ram_gb": r["peak_ram_gb"],
            "gpu_h": (r["wall_s"] / 3600) if gpu.lower() != "none" else None,
            "source": f"{source['path']}/timing.tsv",
        })
    out = []
    for (qs,), g in t.group_by(["query_set"], maintain_order=True):
        sub = [r for r in rows if r["target"].startswith(f"{qs}|")]
        summary = pl.DataFrame(long_tables.summarise_tasks(sub, "submission", n_queries[qs]),
                               schema=long_tables.COMPUTE_SCHEMA, orient="row")
        # A reported row names the submitter's machine, not a cluster default.
        per_step = g.group_by("setting", "step").agg(
            pl.col("gpu_model").unique().sort().str.join("; ").alias("gpu_model"),
            pl.col("hardware").unique().sort().str.join("; ").alias("hardware"),
        )
        summary = (summary.join(per_step, on=["setting", "step"], how="left")
                   .with_columns(
                       pl.when(pl.col("gpu_model").str.to_lowercase() == "none").then(pl.lit("none"))
                       .otherwise(pl.col("gpu_model")).alias("gpu_type"),
                       pl.lit("reported").alias("provenance"),
                       pl.lit(qs).alias("query_set_id"),
                   ).drop("gpu_model"))
        out.append(summary)
    return pl.concat(out).sort("query_set_id", "setting", "step")


def fmt(x, d=3) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    return f"{x:.{d}f}"


def leaderboard_now(site: pl.DataFrame, tool: str) -> tuple[dict, dict]:
    """(same tool and setting already on the leaderboard, best of every tool) per key."""
    base = site.filter((pl.col("experiment") == "qfo_region_midi_plus") & (pl.col("stratum_axis") == "all")
                       & (pl.col("metric") == "fmax") & (pl.col("dedup_or_raw") == "raw"))
    same = {(r["target_db"], r["truth_source"], r["setting"]): r["value"]
            for r in base.filter(pl.col("tool") == tool).iter_rows(named=True)}
    best = {}
    for r in base.sort("value", descending=True).iter_rows(named=True):
        best.setdefault((r["target_db"], r["truth_source"]), (r["value"], r["tool"], r["setting"]))
    return same, best


def comment(rows: pl.DataFrame, comp: pl.DataFrame, card: dict, site: pl.DataFrame,
            source: dict, n_pred: int) -> str:
    head = rows.filter(pl.col("stratum_axis") == "all").with_columns(
        (pl.col("dedup_or_raw") + "_" + pl.col("metric")).alias("col"))
    wide = head.pivot(on="col", index=["query_set_id", "target_db",
                      "target_mya_to_human", "truth_source", "setting"], values="value",
                      aggregate_function="first")
    n_inst = head.group_by("query_set_id", "target_db", "truth_source", "setting").agg(pl.col("n_instances").max())
    wide = wide.join(n_inst, on=["query_set_id", "target_db", "truth_source", "setting"], how="left").sort(
        "query_set_id", "target_mya_to_human", "truth_source", "setting", nulls_last=True)
    same, best = leaderboard_now(site, card["tool"])
    L = []
    L.append(f"## Scored: {card['name']} {card['version']}")
    L.append("")
    L.append(f"The leaderboard's scorer graded {n_pred:,} calls from `{source['path']}` "
             f"at commit `{source['commit'][:8]}`. Every number below is measured by that "
             "scorer. The cost numbers are the submitter's own, marked as reported.")
    L.append("")
    L.append("Table 1. Scores of the submitted calls on the whole query set. Fmax: the best balance "
             "of precision and recall over all score cut-offs (1 is perfect). AUPRC: area under "
             "the precision-recall curve. Median IoU: how well a correct call's start and end "
             "match the true region (1 is exact). Raw keeps every call; dedup merges overlapping "
             "calls of one family on one protein. The last two columns are the leaderboard as it "
             "is now (raw Fmax): the same tool and setting if it is already there, and the best "
             "tool for that relative and truth set.")
    L.append("")
    L.append("| relative (Mya) | truth set | setting | truth regions | Fmax raw | Fmax dedup "
             "| AUPRC raw | median IoU raw | same tool on the leaderboard | best on the leaderboard |")
    L.append("|---|---|---|---:|---:|---:|---:|---:|---:|---|")
    for r in wide.iter_rows(named=True):
        key = (r["target_db"], r["truth_source"])
        s = same.get((r["target_db"], r["truth_source"], r["setting"]))
        b = best.get(key)
        mya = "" if r["target_mya_to_human"] is None else f" ({round(r['target_mya_to_human'])})"
        L.append(
            f"| {r['target_db']}{mya} | {r['truth_source']} | `{r['setting']}` | {r['n_instances']:,} "
            f"| {fmt(r.get('raw_fmax'))} | {fmt(r.get('dedup_fmax'))} | {fmt(r.get('raw_auprc'))} "
            f"| {fmt(r.get('raw_median_iou_tp'))} | {fmt(s)} "
            f"| {fmt(b[0]) + ' ' + b[1] + ' `' + b[2] + '`' if b else '—'} |")
    L.append("")
    L.append("Table 2. Cost, as reported by the submitter (provenance = reported). Core-hours are "
             "cores × wall-clock time; they are not measured by the leaderboard.")
    L.append("")
    L.append("| query set | setting | step | runs | wall-clock s | cores | core-h | peak RAM GB | GPU | hardware |")
    L.append("|---|---|---|---:|---:|---:|---:|---:|---|---|")
    for r in comp.iter_rows(named=True):
        L.append(f"| {r['query_set_id']} | `{r['setting']}` | {r['step']} | {r['n_tasks']} "
                 f"| {fmt(r['wall_s'], 1)} "
                 f"| {r['cores'] if r['cores'] is not None else '—'} | {fmt(r['cpu_h'], 4)} "
                 f"| {fmt(r['peak_ram_gb'], 2)} | {r['gpu_type']} | {r['hardware']} |")
    L.append("")
    needs = card["needs"]
    L.append("Tool card, built from `tool.yaml`:")
    L.append("")
    L.append(f"- **Compares:** {card['compares']}")
    L.append(f"- **Needs:** sequence only: {'yes' if needs['sequence_only'] else 'no'}; "
             f"3D structure prediction: {'yes' if needs['structure_prediction'] else 'no'}; "
             f"GPU: {'yes' if needs['gpu'] else 'no'}")
    L.append(f"- **Licence:** [{card['licence']}]({card['licence_url']})")
    L.append(f"- **Container:** `{card['container']}`")
    L.append("- **Command:**")
    L.append("")
    L.append("```bash")
    L.append(str(card["command"]).rstrip())
    L.append("```")
    L.append("")
    strata = rows.filter((pl.col("stratum_axis") != "all") & (pl.col("metric") == "fmax"))
    L.append(f"<details><summary>Every row this submission adds: {rows.height:,} score rows "
             f"({strata.height:,} of them Fmax by stratum) and {comp.height} cost rows</summary>")
    L.append("")
    L.append("Fmax per stratum (raw calls). A stratum is a slice of the truth regions: by identity "
             "to the nearest same-family region in the relative, by AlphaFold pLDDT of the query, "
             "by Swiss-Prot feature type, or by feature length in residues.")
    L.append("")
    L.append("```text")
    L.append("relative\ttruth\tsetting\taxis\tbin\tfmax\tn_regions")
    for r in strata.filter(pl.col("dedup_or_raw") == "raw").sort(
            "target_db", "truth_source", "setting", "stratum_axis", "stratum_bin").iter_rows(named=True):
        L.append(f"{r['target_db'].removeprefix('QfO ')}\t{r['truth_source']}\t{r['setting']}\t"
                 f"{r['stratum_axis']}\t{r['stratum_bin']}\t{fmt(r['value'])}\t{r['n_instances']}")
    L.append("```")
    L.append("")
    L.append("All rows, in the layout of `data/results.parquet` and `data/compute.parquet`, are in "
             "this run's `scored-rows` artifact.")
    L.append("</details>")
    L.append("")
    L.append("<sub>Posted by the score-submission workflow. Scorer: `scorer/evaluate_domain_calls.py`, "
             "copied from seanome/2024-kmerseek-analysis at commit 6df20c1 (`scorer/SOURCE.json`).</sub>")
    text = "\n".join(L) + "\n"
    if len(text) > 60000:
        cut = text[:58000].rsplit("\n", 1)[0]
        text = cut + "\n```\n\n(The list is cut here to fit a GitHub comment; the artifact has every row.)\n</details>\n"
    return text


def problem_comment(problems: Problems, sub: Path) -> str:
    L = [f"## Not scored yet: {len(problems)} problem(s) in `{sub}`", "",
         "Fix these and push again; the workflow reruns on every push. CONTRIBUTING.md has the layout.", ""]
    L += [f"{i}. {p}" for i, p in enumerate(problems, 1)]
    if len(problems) == MAX_PROBLEMS:
        L.append(f"\nOnly the first {MAX_PROBLEMS} problems are listed.")
    L += ["", "<sub>Posted by the score-submission workflow.</sub>"]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--submission", type=Path, required=True)
    ap.add_argument("--query-sets", type=Path, default=ROOT / "data/query_sets")
    ap.add_argument("--site-data", type=Path, default=ROOT / "data")
    ap.add_argument("--contributed", type=Path, default=ROOT / "data/truth_contributed",
                    help="folder of <query set>.parquet truth tables built by scripts/truth_contributed.py")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--commit", default="local", help="commit of the pull request's head")
    ap.add_argument("--path", default=None, help="the submission's path as the pull request shows it")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    source = {"path": a.path or str(a.submission), "commit": a.commit}

    problems = Problems()
    card = check_card(a.submission, problems)
    raw = read_predictions(a.submission, problems)
    preds = check_predictions(raw, a.query_sets, card, problems) if raw is not None and card else None
    timing = check_timing(a.submission, preds, problems) if preds is not None and not problems else None
    if problems:
        (a.out / "comment.md").write_text(problem_comment(problems, source["path"]))
        print((a.out / "comment.md").read_text())
        return 1

    work = a.out / "work"
    # A run's metrics are found by globbing this folder, so files from an earlier run here
    # would be counted as this run's.
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir()
    contributed = a.contributed if a.contributed.exists() else None
    try:
        metrics = run_scorer(work, a.query_sets, preds, card, contributed)
    except subprocess.CalledProcessError as e:
        (a.out / "comment.md").write_text(
            "## Not scored: the scorer stopped\n\nThe files passed every check, but "
            f"`scorer/evaluate_domain_calls.py` exited with status {e.returncode}. The workflow "
            "log has its message. This is a problem on the leaderboard's side; a maintainer will look.\n\n"
            "<sub>Posted by the score-submission workflow.</sub>\n")
        raise
    rows = to_rows(metrics, a.query_sets, card, source)
    n_queries = {qs: json.loads((a.query_sets / qs / "SOURCE.json").read_text())["n_queries"]
                 for qs in preds["query_set"].unique().to_list()}
    comp = to_compute(timing, card, n_queries, source)
    rows.write_parquet(a.out / "results.parquet")
    comp.write_parquet(a.out / "compute.parquet")
    card_out = {k: card[k] for k in CARD_KEYS} | {
        "numbers_from": {"submission": source["path"], "commit": source["commit"]},
        "cost_provenance": "reported",
    }
    (a.out / "card.json").write_text(json.dumps(card_out, indent=1) + "\n")
    site = pl.read_parquet(a.site_data / "results.parquet")
    (a.out / "comment.md").write_text(comment(rows, comp, card, site, source, preds.height))
    print(f"{rows.height} score rows, {comp.height} cost rows -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
