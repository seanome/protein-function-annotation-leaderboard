#!/usr/bin/env python3
"""Resolve claims/*.yaml against the tables and the experiment inventory.

    python scripts/claims.py --data data --check      # validate; exit 1 on any problem
    python scripts/claims.py --data data --out docs   # also write docs/data/claims/<tool>.json
                                                      # and docs/badges/claims-<tool>.svg

A claim's status comes from its test, never from the YAML: alive when the test runs and does
not kill the claim, retired when it does, untested when the test was not run or the data
cannot decide it (too few instances). A `record` quote must equal the inventory's value
character for character, so a re-run notebook that changes a number fails the check here
instead of leaving a stale claim on the site.

A claim about another group's tool quotes the authors (`author_claim`). Its sha256 was written
by scripts/fetch_quote.py when the quote was found on the authors' page; a quote that no longer
matches it fails the check, so a quote cannot be reworded by hand.

A claim retired by a record test is dated with that record's date. A claim retired by a
table test has no date: the tables are recomputed on every build, so the build cannot tell
when the kill first held. The verdict line is the test's optional `verdict` template,
filled from the resolved test; without one it is the test's own text.
"""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import jsonschema
import polars as pl
import yaml

ROOT = Path(__file__).resolve().parents[1]


def num(text: str) -> float:
    m = re.search(r"-?\d+(?:\.\d+)?(?:e-?\d+)?", text)
    if not m:
        raise ValueError(f"no number in {text!r}")
    return float(m.group(0))


def signed(x: float, d: int = 3) -> str:
    # A value that rounds to zero prints as +0.000, not −0.000.
    if round(x, d) == 0:
        x = 0.0
    return f"{x:+.{d}f}".replace("-", "−")


class Resolver:
    def __init__(self, data: Path):
        self.long = pl.read_parquet(data / "results.parquet")
        self.src = json.loads((data / "SOURCE.json").read_text())
        raw = (data / "experiments.json").read_bytes()
        # Where data/experiments.json was copied from; the copy must match its checksum.
        self.inv_src = json.loads((data / "experiments.SOURCE.json").read_text())
        if hashlib.sha256(raw).hexdigest() != self.inv_src["sha256"]:
            raise ValueError("data/experiments.json does not match the sha256 in data/experiments.SOURCE.json")
        self.records = {str(r["id"]): r for r in json.loads(raw)}

    def inventory_words(self) -> str:
        s = self.inv_src
        return f"{s['repo']} {s['path']} @ {s['commit'][:8]} (PR {s['pr']})"

    # ---- tables ----
    def _rows(self, t: dict) -> pl.DataFrame:
        q = self.long.filter(pl.col("experiment") == t["experiment"], pl.col("metric") == t["metric"])
        for k in ("truth_source", "dedup_or_raw", "stratum_axis", "stratum_bin"):
            if k in t:
                q = q.filter(pl.col(k) == t[k])
        return q

    def _diffs(self, t: dict, setting: str, dedup: str | None = None) -> dict:
        q = self._rows(t)
        if dedup:
            q = q.filter(pl.col("dedup_or_raw") == dedup)
        theirs = q.filter(pl.col("tool") == t["against"]["tool"], pl.col("setting") == t["against"]["setting"])
        ours = q.filter(pl.col("tool") == t["tool"])
        out, missing, small = [], [], []
        for target in t["targets"]:
            th = theirs.filter(pl.col("target_db") == target)
            o = ours.filter(pl.col("target_db") == target)
            if setting != "best_per_target":
                o = o.filter(pl.col("setting") == setting)
            if th.height != 1 or o.height == 0:
                missing.append(target)
                continue
            # A relative with too few truth regions in this stratum cannot decide: one call
            # swings Fmax from 0 to 1. It is left out and named, never counted.
            if "min_instances" in t and min(int(th["n_instances"][0]), int(o["n_instances"].min())) < t["min_instances"]:
                small.append(target)
                continue
            out.append((target, float(o["value"].max()) - float(th["value"][0])))
        return {"diffs": out, "missing": missing, "small": small}

    def table(self, t: dict) -> dict:
        r = self._diffs(t, t["setting"])
        small = r["small"]
        if r["missing"] or not r["diffs"]:
            why = f"no rows for {', '.join(r['missing']) or 'any target'}"
            if small and not r["missing"]:
                why = (f"every relative has fewer than {t['min_instances']} truth regions in this stratum; "
                       f"at least {t['min_instances']} are needed to decide")
            return {"status": "untested", "text": why}
        ds = [d for _, d in r["diffs"]]
        k = t["kill_if"]
        if "abs_difference_above" in k:
            killed = any(abs(d) > k["abs_difference_above"] for d in ds)
            worst = max(r["diffs"], key=lambda x: abs(x[1]))
        else:
            # One-sided: only a shortfall against the other tool disproves the claim.
            killed = any(d < -k["below_by_more_than"] for d in ds)
            worst = min(r["diffs"], key=lambda x: x[1])
        text = f"{signed(min(ds))} to {signed(max(ds))} across {len(ds)} of {len(t['targets'])}"
        if small:
            text += (f" ({len(small)} with fewer than {t['min_instances']} truth regions left out: "
                     + ", ".join(x.removeprefix("QfO ") for x in small) + ")")
        if killed:
            text += f"; furthest: {worst[0].removeprefix('QfO ')} {signed(worst[1])}"
        return {"status": "retired" if killed else "alive", "text": text,
                "source": f"{self.src['repo']}@{self.src['commit'][:8]} results.parquet"}

    def table_count_passing(self, t: dict) -> dict:
        q = self._rows(t).filter(pl.col("tool") == t["tool"])
        thr = t["kill_if"]["abs_difference_above"]
        passing, total = [], 0
        for setting in sorted(q["setting"].unique().to_list()):
            for dedup in ("raw", "dedup"):
                r = self._diffs(t, setting, dedup)
                if r["missing"]:
                    continue
                total += 1
                if all(abs(d) <= thr for _, d in r["diffs"]):
                    passing.append(f"{setting} {dedup}")
        n_settings = q["setting"].n_unique()
        return {"text": f"{len(passing)} of {total} setting and dedup choices pass ({n_settings} settings): "
                        + (", ".join(passing) or "none")}

    def table_bins(self, t: dict) -> dict:
        q = self._rows(t).filter(pl.col("tool") == t["tool"], pl.col("stratum_bin").is_in(t["stratum_bins"]),
                                 pl.col("setting").str.starts_with(t["setting_prefix"]))
        if q.height == 0:
            return {"status": "untested", "text": "no rows in those bins"}
        # Only rows with enough regions to decide take part: in a bin of one region, Fmax
        # swings from 0 to 1 on a single call.
        big = q.filter(pl.col("n_instances") >= t["min_instances"])
        if big.height == 0:
            n_max = int(q["n_instances"].max())
            return {"status": "untested",
                    "text": f"the {', '.join(t['stratum_bins'])} bins hold at most {n_max} "
                            f"region{'s' if n_max != 1 else ''} for any setting and relative; "
                            f"at least {t['min_instances']} are needed to decide"}
        best = float(big["value"].max())
        killed = best < t["kill_if"]["all_below"]
        return {"status": "retired" if killed else "alive",
                "text": f"best Fmax {best:.3f} over the {big.height} setting, relative and bin "
                        f"combinations with at least {t['min_instances']} regions"}

    def table_max(self, t: dict) -> dict:
        q = self._rows(t).filter(pl.col("tool") == t["tool"]).sort(["value", "setting"], descending=[True, False])
        if q.height == 0:
            return {"status": "untested", "text": "no rows for this test"}
        r = q.row(0, named=True)
        return {"text": f"{r['value']:.3f} ({r['setting']}, n = {r['n_instances']:,} queries)"}

    # ---- inventory records ----
    def _record(self, row: str, name: str, quote: str) -> tuple[dict, dict]:
        rec = self.records.get(row)
        if rec is None:
            raise ValueError(f"inventory has no row {row!r}")
        res = next((x for x in rec.get("result") or [] if x.get("name") == name), None)
        if res is None:
            raise ValueError(f"row {row} has no result named {name!r}")
        if str(res["value"]).strip() != quote:
            raise ValueError(f"row {row} {name!r}: quote {quote!r} != inventory value {res['value']!r}")
        return rec, res

    def record(self, t: dict) -> dict:
        rec, res = self._record(t["row"], t["name"], t["quote"])
        return {"text": t["quote"], "date": rec.get("date"), "source": self._src(rec, res)}

    def record_compare(self, t: dict) -> dict:
        rec, ro = self._record(t["row"], t["ours"]["name"], t["ours"]["quote"])
        _, rt = self._record(t["row"], t["theirs"]["name"], t["theirs"]["quote"])
        a, b = num(t["ours"]["quote"]), num(t["theirs"]["quote"])
        k = t["kill_if"]
        killed = (k.get("ours_below_theirs") and a < b) or ("ours_at_most" in k and a <= k["ours_at_most"])
        out = {"status": "retired" if killed else "alive", "a": a, "b": b, "a_plus_b": a + b,
               "date": rec.get("date"), "source": self._src(rec, ro)}
        if "total" in t:
            # A denominator quoted from a field of the same record, checked like any quote.
            field, quote = t["total"]["field"], t["total"]["quote"]
            if quote not in str(rec.get(field, "")):
                raise ValueError(f"row {t['row']} field {field!r} does not contain {quote!r}")
            out["total"] = num(quote.replace(",", ""))
            out["a_plus_b_share"] = (a + b) / out["total"]
        return out

    @staticmethod
    def _src(rec: dict, res: dict) -> str:
        s = rec.get("source") or {}
        return f"{s.get('path', '')} @ {s.get('commit', '')}, cell {res.get('source_cell', '?')}"

    def resolve(self, t: dict) -> dict:
        if t["kind"] == "not_run":
            return {"status": "untested", "text": t["reason"]}
        return getattr(self, t["kind"])(t)


def verdict_text(t: dict, r: dict) -> str:
    """The verdict line: the test's `verdict` template filled from the resolved test, or the
    test's own text when it has no template."""
    if "verdict" not in t:
        return r.get("text", "")
    try:
        return t["verdict"].format_map(r)
    except (KeyError, ValueError) as e:
        raise ValueError(f"verdict template {t['verdict']!r} cannot be filled from {sorted(r)}: {e}") from e


def normalise(text: str) -> str:
    """Whitespace collapsed to single spaces, so a quote that wraps across lines in the
    source still matches. Nothing else changes: case, punctuation and words must match."""
    return " ".join(text.split())


def quote_sha256(quote: str) -> str:
    return hashlib.sha256(normalise(quote).encode()).hexdigest()


def check_author_quote(claim: dict) -> None:
    """The author's words must be the words recorded when their page was fetched.
    `sha256` is written by scripts/fetch_quote.py from the quote as found on the page, so a
    quote edited by hand, without fetching the page again, fails here."""
    a = claim.get("author_claim")
    if not a:
        return
    if quote_sha256(a["quote"]) != a["sha256"]:
        raise ValueError(f"claim {claim['id']!r}: the quote does not match the sha256 recorded when "
                         f"{a['text_url']} was fetched on {a['fetched']}; run scripts/fetch_quote.py again")
    n_words = len(normalise(a["quote"]).split())
    if n_words >= 15:
        raise ValueError(f"claim {claim['id']!r}: the quote has {n_words} words; keep it under 15")


def resolve_file(path: Path, res: Resolver, schema: dict) -> dict:
    doc = yaml.safe_load(path.read_text())
    jsonschema.validate(doc, schema)
    for c in doc["claims"]:
        check_author_quote(c)
    out = {"tool": doc["tool"], "title": doc["title"], "updated": str(doc["updated"]),
           "file": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path), "claims": []}
    for c in doc["claims"]:
        r = res.resolve(c["test"])
        ctx = []
        for x in c.get("context", []):
            rx = res.resolve(x["test"])
            ctx.append({"label": x["label"], "text": verdict_text(x["test"], rx), "source": rx.get("source")})
        out["claims"].append({
            "id": c["id"], "claim": c["claim"], "kill": c["kill"], "scope": c.get("scope"),
            "status": r["status"], "verdict": verdict_text(c["test"], r), "date": r.get("date"),
            "source": r.get("source"), "context": ctx, "evidence": c.get("evidence"),
            **({"author_claim": {k: c["author_claim"][k] for k in ("quote", "url", "fetched")}}
               if "author_claim" in c else {}),
        })
    counts = {s: sum(x["status"] == s for x in out["claims"]) for s in ("alive", "retired", "untested")}
    kills = [x["date"] for x in out["claims"] if x["status"] == "retired" and x["date"]]
    out["counts"] = counts
    out["last_kill"] = max(kills) if kills else None
    out["source"] = {"tables": {k: res.src[k] for k in ("repo", "commit", "commit_date")},
                     "inventory": res.inventory_words()}
    return out


def badge(tool: str, counts: dict) -> str:
    """A flat badge, ink on the site's paper colour, with no status colours."""
    left = "claims"
    right = f"{counts['alive']} alive · {counts['retired']} retired · {counts['untested']} untested"
    wl, wr = 7 * len(left) + 12, int(6.4 * len(right)) + 14
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{wl + wr}" height="20" role="img" '
            f'aria-label="{tool} {left}: {right}"><title>{tool} {left}: {right}</title>'
            f'<rect width="{wl}" height="20" fill="#1C211F"/><rect x="{wl}" width="{wr}" height="20" fill="#F6F7F5" '
            f'stroke="#1C211F"/><g font-family="Verdana,DejaVu Sans,sans-serif" font-size="11">'
            f'<text x="{wl / 2}" y="14" fill="#F6F7F5" text-anchor="middle">{left}</text>'
            f'<text x="{wl + wr / 2}" y="14" fill="#1C211F" text-anchor="middle">{right}</text></g></svg>\n')


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--claims", type=Path, default=ROOT / "claims", help="folder of <tool>.yaml files")
    a = ap.parse_args()
    schema = json.loads((ROOT / "claims/schema.json").read_text())
    res = Resolver(a.data)
    ok = True
    for f in sorted(a.claims.glob("*.yaml")):
        try:
            out = resolve_file(f, res, schema)
        except jsonschema.ValidationError as e:
            where = "/".join(map(str, e.absolute_path)) or "top level"
            print(f"{f.name}: does not match claims/schema.json at {where}: {e.message}")
            ok = False
            continue
        except Exception as e:  # one broken file must not stop the other files being checked
            print(f"{f.name}: {type(e).__name__}: {e}")
            ok = False
            continue
        print(f"{f.name}: {out['counts']}, last kill {out['last_kill']}")
        for c in out["claims"]:
            print(f"  {c['status']:9s} {c['claim']}: {c['verdict']}")
        if a.out:
            (a.out / "data/claims").mkdir(parents=True, exist_ok=True)
            (a.out / "badges").mkdir(parents=True, exist_ok=True)
            (a.out / f"data/claims/{out['tool']}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
            (a.out / f"badges/claims-{out['tool']}.svg").write_text(badge(out["tool"], out["counts"]))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
