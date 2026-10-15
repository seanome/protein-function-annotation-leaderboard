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
"""

import argparse
import json
import re
import sys
from datetime import datetime, timezone
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
    return f"{x:+.{d}f}".replace("-", "−")


class Resolver:
    def __init__(self, data: Path):
        self.long = pl.read_parquet(data / "results.parquet")
        self.src = json.loads((data / "SOURCE.json").read_text())
        inv = json.loads((data / "experiments.json").read_text())
        self.records = {str(r["id"]): r for r in inv}

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
        out, missing = [], []
        for target in t["targets"]:
            th = theirs.filter(pl.col("target_db") == target)
            o = ours.filter(pl.col("target_db") == target)
            if setting != "best_per_target":
                o = o.filter(pl.col("setting") == setting)
            if th.height != 1 or o.height == 0:
                missing.append(target)
                continue
            out.append((target, float(o["value"].max()) - float(th["value"][0])))
        return {"diffs": out, "missing": missing}

    def table(self, t: dict) -> dict:
        r = self._diffs(t, t["setting"])
        if r["missing"] or not r["diffs"]:
            return {"status": "untested", "text": f"no rows for {', '.join(r['missing']) or 'any target'}"}
        ds = [d for _, d in r["diffs"]]
        thr = t["kill_if"]["abs_difference_above"]
        killed = any(abs(d) > thr for d in ds)
        worst = max(r["diffs"], key=lambda x: abs(x[1]))
        text = f"{signed(min(ds))} to {signed(max(ds))} across {len(ds)} of {len(t['targets'])}"
        if killed:
            text += f"; furthest: {worst[0].removeprefix('QfO ')} {signed(worst[1])}"
        return {"status": "retired" if killed else "alive", "text": text, "date": self.src["commit_date"][:10],
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
        n = q.group_by("stratum_bin").agg(pl.col("n_instances").max())
        n_max = int(n["n_instances"].max())
        if n_max < t["min_instances"]:
            return {"status": "untested",
                    "text": f"the under-40% bins hold at most {n_max} region{'s' if n_max != 1 else ''} "
                            f"each; at least {t['min_instances']} are needed to decide"}
        best = float(q["value"].max())
        killed = best < t["kill_if"]["all_below"]
        return {"status": "retired" if killed else "alive", "text": f"best Fmax {best:.3f}",
                "date": self.src["commit_date"][:10]}

    def table_max(self, t: dict) -> dict:
        q = self._rows(t).filter(pl.col("tool") == t["tool"]).sort(["value", "setting"], descending=[True, False])
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
        return {"status": "retired" if killed else "alive", "a": a, "b": b, "date": rec.get("date"),
                "source": self._src(rec, ro)}

    @staticmethod
    def _src(rec: dict, res: dict) -> str:
        s = rec.get("source") or {}
        return f"{s.get('path', '')} @ {s.get('commit', '')}, cell {res.get('source_cell', '?')}"

    def resolve(self, t: dict) -> dict:
        if t["kind"] == "not_run":
            return {"status": "untested", "text": t["reason"]}
        return getattr(self, t["kind"])(t)


def verdict_text(c: dict, r: dict) -> str:
    """The number on the verdict line, worded per claim from the resolved test."""
    if c["id"] == "faster-than-mmseqs2":
        return f"{r['a']:.2f} queries per second against MMseqs2's {r['b']:.1f}"
    if c["id"] == "beats-foldseek-scope":
        # 5,713 is FoldSeek's SCOPe40 benchmark size, as nb 076 states in its verdict basis.
        return (f"fired on {int(r['b'])} of 5,713 benchmark queries ({100 * r['b'] / 5713:.1f}%); "
                f"{int(r['a'])} that FoldSeek missed")
    return r.get("text", "")


def resolve_file(path: Path, res: Resolver, schema: dict) -> dict:
    doc = yaml.safe_load(path.read_text())
    jsonschema.validate(doc, schema)
    out = {"tool": doc["tool"], "title": doc["title"], "updated": str(doc["updated"]),
           "file": str(path.relative_to(ROOT)), "claims": []}
    for c in doc["claims"]:
        r = res.resolve(c["test"])
        ctx = []
        for x in c.get("context", []):
            rx = res.resolve(x["test"])
            ctx.append({"label": x["label"], "text": rx.get("text") or verdict_text(c, rx), "source": rx.get("source")})
        out["claims"].append({
            "id": c["id"], "claim": c["claim"], "kill": c["kill"], "scope": c.get("scope"),
            "status": r["status"], "verdict": verdict_text(c, r), "date": r.get("date"),
            "source": r.get("source"), "context": ctx, "evidence": c.get("evidence"),
        })
    counts = {s: sum(x["status"] == s for x in out["claims"]) for s in ("alive", "retired", "untested")}
    kills = [x["date"] for x in out["claims"] if x["status"] == "retired" and x["date"]]
    out["counts"] = counts
    out["last_kill"] = max(kills) if kills else None
    out["source"] = {"tables": {k: res.src[k] for k in ("repo", "commit", "commit_date")},
                     "inventory": "seanome/2024-kmerseek-analysis data/experiments.json @ fe72334b (PR 97)"}
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
    a = ap.parse_args()
    schema = json.loads((ROOT / "claims/schema.json").read_text())
    res = Resolver(a.data)
    ok = True
    for f in sorted((ROOT / "claims").glob("*.yaml")):
        try:
            out = resolve_file(f, res, schema)
        except (ValueError, jsonschema.ValidationError) as e:
            print(f"{f.name}: {e}")
            ok = False
            continue
        print(f"{f.name}: {out['counts']}, last kill {out['last_kill']}")
        for c in out["claims"]:
            print(f"  {c['status']:9s} {c['claim']}: {c['verdict']}")
        if a.out:
            (a.out / "data/claims").mkdir(parents=True, exist_ok=True)
            (a.out / "badges").mkdir(parents=True, exist_ok=True)
            out["resolved_at"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            (a.out / f"data/claims/{out['tool']}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
            (a.out / f"badges/claims-{out['tool']}.svg").write_text(badge(out["tool"], out["counts"]))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
