#!/usr/bin/env python3
"""Regression checks for scripts/claims.py on small made-up tables.

1. table_bins takes its best Fmax only from rows with at least min_instances regions, and is
   untested when no row has that many.
2. table and table_bins kills carry no date; record kills carry the record's date.
3. table_max on a filter with no rows is untested, not a crash.
4. A misspelt key in a test fails the schema.
5. main() reports a broken file and still checks the files after it.
6. A `verdict` template is filled from the resolved test.
7. `--out` twice writes byte-identical JSON, with no wall-clock time in it.
"""

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import jsonschema
import polars as pl
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import claims  # noqa: E402

SCHEMA = json.loads((ROOT / "claims/schema.json").read_text())
failures: list[str] = []


def check(cond: bool, what: str) -> None:
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        failures.append(what)


def rows(*specs) -> pl.DataFrame:
    base = dict(experiment="e", metric="fmax", truth_source="Pfam-N domain instances", dedup_or_raw="raw",
                stratum_axis="identity", tool="kmerseek")
    return pl.DataFrame([dict(base, setting=s, stratum_bin=b, n_instances=n, value=v, target_db=t)
                         for s, b, n, v, t in specs])


def resolver(long: pl.DataFrame, records=None) -> claims.Resolver:
    r = claims.Resolver.__new__(claims.Resolver)
    r.long = long
    r.src = {"repo": "x/y", "commit": "0" * 40, "commit_date": "2026-10-15T00:00:00Z"}
    r.records = records or {}
    r.inv_src = {"repo": "x/y", "path": "data/experiments.json", "commit": "1" * 40, "pr": 1}
    return r


BINS = dict(kind="table_bins", experiment="e", truth_source="Pfam-N domain instances", dedup_or_raw="raw",
            stratum_axis="identity", stratum_bins=["0-20%", "20-30%"], metric="fmax", tool="kmerseek",
            setting_prefix="hp_", min_instances=30, kill_if={"all_below": 0.05})

# 1. A one-region bin with Fmax 1.0 must not keep the claim alive when the 40-region bin is 0.01.
long = rows(("hp_k24", "0-20%", 1, 1.0, "QfO fly"), ("hp_k24", "20-30%", 40, 0.01, "QfO fly"))
r = resolver(long).table_bins(BINS)
check(r["status"] == "retired", f"table_bins ignores bins under min_instances: {r}")
check("date" not in r, "table_bins kill has no date")
r = resolver(rows(("hp_k24", "0-20%", 1, 1.0, "QfO fly"), ("hp_k24", "20-30%", 29, 0.9, "QfO fly"))).table_bins(BINS)
check(r["status"] == "untested" and "at most 29" in r["text"], f"table_bins untested when no row is big enough: {r}")
r = resolver(rows(("hp_k24", "20-30%", 40, 0.2, "QfO fly"), ("hp_k27", "0-20%", 35, 0.04, "QfO fly"))).table_bins(BINS)
check(r["status"] == "alive" and "0.200" in r["text"], f"table_bins alive on a big bin above the cutoff: {r}")

# 2. table kill: no date.
tab = dict(kind="table", experiment="e", truth_source="Pfam-N domain instances", dedup_or_raw="raw",
           stratum_axis="identity", stratum_bin="0-20%", metric="fmax", tool="kmerseek", setting="hp_k24",
           against={"tool": "other", "setting": "d"}, targets=["QfO fly"], kill_if={"abs_difference_above": 0.05})
long = pl.concat([rows(("hp_k24", "0-20%", 40, 0.1, "QfO fly")),
                  rows(("d", "0-20%", 40, 0.5, "QfO fly")).with_columns(pl.lit("other").alias("tool"))])
r = resolver(long).table(tab)
check(r["status"] == "retired" and "date" not in r, f"table kill has no date: {r}")
rec = {"7": {"id": "7", "date": "2026-08-30", "dataset": "benchmark 5,713 queries", "source": {},
             "result": [{"name": "only", "value": "0"}, {"name": "both", "value": "188"}]}}
rc = dict(kind="record_compare", row="7", ours={"name": "only", "quote": "0"}, theirs={"name": "both", "quote": "188"},
          total={"field": "dataset", "quote": "5,713"}, kill_if={"ours_at_most": 0},
          verdict="fired on {a_plus_b:,.0f} of {total:,.0f} ({a_plus_b_share:.1%}); {a:,.0f} missed")
r = resolver(pl.DataFrame(), rec).record_compare(rc)
check(r["status"] == "retired" and r["date"] == "2026-08-30", f"record kill dated by its record: {r}")
# 6. verdict template, with only-ours + both in the fired count.
rec["7"]["result"][0]["value"] = "12"
rc["ours"]["quote"] = "12"
r = resolver(pl.DataFrame(), rec).record_compare(rc)
v = claims.verdict_text(rc, r)
check(v == "fired on 200 of 5,713 (3.5%); 12 missed", f"verdict counts only-ours + both: {v!r}")
bad_total = dict(rc, total={"field": "dataset", "quote": "9,999"})
try:
    resolver(pl.DataFrame(), rec).record_compare(bad_total)
    check(False, "a total that is not in the record fails")
except ValueError:
    check(True, "a total that is not in the record fails")

# 3. table_max on no rows.
r = resolver(rows(("hp_k24", "0-20%", 40, 0.1, "QfO fly"))).table_max(
    dict(kind="table_max", experiment="nope", metric="fmax", tool="kmerseek"))
check(r["status"] == "untested", f"table_max with no rows is untested: {r}")

# 4. a misspelt filter key fails the schema.
doc = {"tool": "t", "title": "t", "updated": "2026-10-16",
       "claims": [{"id": "a", "claim": "c", "kill": "k", "test": dict(tab, stratum_bins_typo=["x"])}]}
try:
    jsonschema.validate(doc, SCHEMA)
    check(False, "misspelt test key fails the schema")
except jsonschema.ValidationError:
    check(True, "misspelt test key fails the schema")
doc["claims"][0]["test"] = tab
jsonschema.validate(doc, SCHEMA)

# 5 and 7. main(): a broken file first, a good one after it; --out twice is identical.
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    (tmp / "claims").mkdir()
    (tmp / "claims/a_broken.yaml").write_text(yaml.safe_dump(
        {"tool": "broken", "title": "b", "updated": "2026-10-16",
         "claims": [{"id": "x", "claim": "c", "kill": "k",
                     "test": {"kind": "record", "row": "no-such-row", "name": "n", "quote": "1"}}]}))
    good = yaml.safe_load((ROOT / "claims/kmerseek.yaml").read_text())
    (tmp / "claims/b_good.yaml").write_text(yaml.safe_dump(good, sort_keys=False))
    cmd = [sys.executable, str(ROOT / "scripts/claims.py"), "--data", str(ROOT / "data"),
           "--claims", str(tmp / "claims")]
    p = subprocess.run(cmd + ["--check"], capture_output=True, text=True)
    check(p.returncode == 1 and "a_broken.yaml: ValueError" in p.stdout and "b_good.yaml: {" in p.stdout,
          f"main reports the broken file and checks the next: rc {p.returncode}, {p.stdout.splitlines()[:2]}")
    outs = []
    for k in (1, 2):
        subprocess.run(cmd + ["--out", str(tmp / f"out{k}")], capture_output=True, text=True)
        outs.append(hashlib.sha256((tmp / f"out{k}/data/claims/kmerseek.json").read_bytes()).hexdigest())
    body = json.loads((tmp / "out1/data/claims/kmerseek.json").read_text())
    check(outs[0] == outs[1] and "resolved_at" not in body, "--out twice gives identical JSON with no clock time")
    check(body["source"]["inventory"].endswith("(PR 97)") and "fe72334b" in body["source"]["inventory"],
          f"inventory source read from data/experiments.SOURCE.json: {body['source']['inventory']}")

# 8. An author's quote altered by one word fails; the quote as fetched passes.
for f in sorted((ROOT / "claims").glob("*.yaml")):
    for c in yaml.safe_load(f.read_text())["claims"]:
        if "author_claim" not in c:
            continue
        claims.check_author_quote(c)  # the committed quote passes
        words = c["author_claim"]["quote"].split()
        altered = dict(c, author_claim=dict(c["author_claim"], quote=" ".join(words[:-1] + ["probably"])))
        try:
            claims.check_author_quote(altered)
            check(False, f"{f.name} {c['id']}: an altered quote fails")
        except ValueError:
            check(True, f"{f.name} {c['id']}: an altered quote fails")
        # Line breaks inside a quote are not a change.
        claims.check_author_quote(dict(c, author_claim=dict(c["author_claim"], quote="\n  ".join(words))))
n_quoted = sum("author_claim" in c for f in (ROOT / "claims").glob("*.yaml")
               for c in yaml.safe_load(f.read_text())["claims"])
check(n_quoted >= 4, f"at least 4 claims quote their authors: {n_quoted}")

# 9. One-sided kill: a tool far ABOVE the other is not disproved; far below is.
one = dict(tab, kill_if={"below_by_more_than": 0.05}, stratum_bin="0-20%")
up = pl.concat([rows(("hp_k24", "0-20%", 40, 0.9, "QfO fly")),
                rows(("d", "0-20%", 40, 0.5, "QfO fly")).with_columns(pl.lit("other").alias("tool"))])
check(resolver(up).table(one)["status"] == "alive", "below_by_more_than ignores a lead")
check(resolver(long).table(one)["status"] == "retired", "below_by_more_than kills a shortfall")

# 10. min_instances in a table test: a small relative is left out and named; all small is untested.
two = pl.concat([rows(("hp_k24", "0-20%", 40, 0.5, "QfO fly"), ("hp_k24", "0-20%", 5, 0.0, "QfO worm")),
                 rows(("d", "0-20%", 40, 0.5, "QfO fly"), ("d", "0-20%", 5, 0.9, "QfO worm"))
                 .with_columns(pl.lit("other").alias("tool"))])
mi = dict(one, targets=["QfO fly", "QfO worm"], min_instances=30)
r = resolver(two).table(mi)
check(r["status"] == "alive" and "left out: worm" in r["text"], f"small relative left out and named: {r}")
r = resolver(two).table(dict(mi, targets=["QfO worm"]))
check(r["status"] == "untested", f"only small relatives: untested: {r}")

if failures:
    print(f"{len(failures)} FAILURES")
    sys.exit(1)
print("all claims checks passed")
