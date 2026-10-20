#!/usr/bin/env python3
"""Checks of scripts/build_site.py and scripts/copy_tables.py that need no browser.

1. Both MMseqs2 searches carry the MMseqs2 database build, whose rows are stored under
   tool "mmseqs2"; the expected cost is read from compute.parquet here, not through
   build_site.py.
2. A database-build row that no search tool maps to stops the build.
3. strong_weak() never lists one stratum as both strong and weak (3 comparisons).
4. Tool tags come from each tool's licence_tag.
5. Two builds from the same data/, under different Python hash seeds, are byte-identical.
6. copy_tables.py deletes an optional file the new commit lacks, and names a file that
   provenance lists but that was not copied.
7. Every tool card's "Can its weak hits be trusted?" block shows the weak-hit table's numbers,
   recomputed here, or says the tool was not run on the dark set.
"""

import argparse
import filecmp
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import build_site  # noqa: E402
import copy_tables  # noqa: E402
from tools import TOOLS  # noqa: E402

failures: list[str] = []


def check(cond: bool, what: str) -> None:
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        failures.append(what)


def diff_tree(a: Path, b: Path) -> list[str]:
    out = []
    cmp = filecmp.dircmp(a, b)
    out += [str(a / x) for x in cmp.left_only + cmp.right_only + cmp.diff_files + cmp.funny_files]
    # dircmp compares by size and time unless shallow=False; recheck the "same" ones by content.
    _, mismatch, errors = filecmp.cmpfiles(a, b, cmp.same_files, shallow=False)
    out += [str(a / x) for x in mismatch + errors]
    for sub in cmp.common_dirs:
        out += diff_tree(a / sub, b / sub)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    a = ap.parse_args()

    # 1. MMseqs2 build cost on both MMseqs2 cards.
    c = pl.read_parquet(a.data / "compute.parquet").filter(pl.col("experiment") == build_site.EXPERIMENT)
    want = c.filter((pl.col("tool") == "mmseqs2") & (pl.col("setting") == "database build"))["cpu_h"].sum()
    _, comp, *_ = build_site.load(a.data)
    for t in ("mmseqs2_seqseq", "mmseqs2_iterative"):
        got = comp.filter(pl.col("tool") == t)["db_cpu_h"].unique().to_list()
        check(want > 0 and got == [want], f"{t} carries the MMseqs2 build: {got} core-hours, want {want}")

    # 2. An orphan build row stops the build.
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for f in ("results.parquet", "SOURCE.json", "provenance.json", "coverage.md"):
            (tmp / f).symlink_to((a.data / f).resolve())
        full = pl.read_parquet(a.data / "compute.parquet")
        orphan = full.filter((pl.col("tool") == "mmseqs2") & (pl.col("experiment") == build_site.EXPERIMENT)) \
            .with_columns(pl.lit("newtool").alias("tool"))
        pl.concat([full, orphan]).write_parquet(tmp / "compute.parquet")
        try:
            build_site.load(tmp)
            check(False, "a build row for an unknown tool stops the build")
        except SystemExit as e:
            check("newtool" in str(e), f"a build row for an unknown tool stops the build: {e}")

    # 3. strong and weak never share a stratum.
    cs = [{"axis": "identity", "stratum": s, "margin": m} for s, m in (("a", 0.3), ("b", 0.1), ("c", -0.2))]
    strong, weak = build_site.strong_weak(cs)
    check([x["stratum"] for x in strong] == ["a", "b"] and [x["stratum"] for x in weak] == ["c"],
          f"3 comparisons: strong {[x['stratum'] for x in strong]}, weak {[x['stratum'] for x in weak]}")
    cs.append({"axis": "identity", "stratum": "d", "margin": -0.5})
    strong, weak = build_site.strong_weak(cs)
    check([x["stratum"] for x in weak] == ["d", "c"], f"4 comparisons: weak {[x['stratum'] for x in weak]}")
    check(build_site.strong_weak(cs[:2])[1] == [], "2 comparisons: no weak list")

    # 4. Licence tags from licence_tag, never from a substring of the licence text.
    for t, info in TOOLS.items():
        f = build_site.tool_filters(t)
        check(f["licence"] == info["licence_tag"] and f["licence"] in build_site.LICENCE_TAGS,
              f"{t}: licence filter {info['licence_tag']!r} ({info['licence']})")

    # 4b. Figure 1 of the tool cards: bar height is proportional to Fmax on the shared scale,
    # so 0.48 and 0.9 do not draw the same (the old dots stopped growing at about 0.48).
    import re as _re
    rel = ["a", "b", "c"]
    build_site.SPECIES.update({"a": "A", "b": "B", "c": "C"})
    build_site.SPECIES_MYA.update({"a": 1, "b": 2, "c": 3})
    svg = build_site.fmax_strip("hmmer3_phmmer", {"hmmer3_phmmer": {"a": 0.24, "b": 0.48, "c": 0.9}},
                                rel, 1.0)
    hs = [float(h) for h in _re.findall(r'<rect [^>]*height="([0-9.]+)"', svg)]
    check(len(hs) == 3 and abs(hs[1] / hs[0] - 2) < 0.02 and abs(hs[2] / hs[1] - 0.9 / 0.48) < 0.02,
          f"bar heights proportional to Fmax: {hs}")
    svg = build_site.fmax_strip("hmmer3_phmmer", {"hmmer3_phmmer": {"a": 0.3}}, rel, 1.0)
    check(svg.count('class="st-nr"') == 2, "a relative with no row gets the grey x, not a bar")

    # 5. Deterministic builds.
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for seed in ("1", "2"):
            subprocess.run([sys.executable, str(ROOT / "scripts/build_site.py"), "--data", str(a.data),
                            "--out", str(tmp / seed)], check=True, capture_output=True,
                           env=dict(os.environ, PYTHONHASHSEED=seed), cwd=ROOT)
        # Pages copy the assets folder only in the repo's docs/, so compare what the build wrote.
        d = diff_tree(tmp / "1", tmp / "2")
        n = sum(1 for _ in (tmp / "1").rglob("*") if _.is_file())
        check(not d, f"two builds with hash seeds 1 and 2: {n} files, {len(d)} differ {d[:3]}")

    # 6. copy_tables.py with gh replaced by files in a folder.
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        out = tmp / "data"
        out.mkdir()
        (out / "inventory.md").write_text("old copy from an earlier commit\n")
        (out / "pileup").mkdir()
        (out / "pileup/regions.parquet").write_bytes(b"old")
        src = {name: (a.data / name).read_bytes() for name in copy_tables.FILES}

        def fake_gh(*args, raw=False, served=src):
            path = args[0]
            if "/commits/" in path:
                return json.dumps({"sha": "a" * 40, "commit": {"committer": {"date": "2026-10-16T00:00:00Z"}}}).encode()
            name = path.split("results/leaderboard/", 1)[1].split("?", 1)[0]
            if name not in served:
                raise subprocess.CalledProcessError(1, "gh")
            return served[name]

        argv = sys.argv
        copy_tables.gh = fake_gh
        try:
            sys.argv = ["copy_tables.py", "--repo", "x/y", "--commit", "abc", "--out", str(out)]
            copy_tables.main()
            check(not (out / "inventory.md").exists() and not (out / "pileup/regions.parquet").exists(),
                  "optional files the commit lacks are deleted, not left to be credited to it")
            s = json.loads((out / "SOURCE.json").read_text())
            check(s["pileup_commit"] is None, "SOURCE.json has no pileup_commit without pile-up files")
            prov = json.loads(src["provenance.json"])
            prov["outputs"]["extra.parquet"] = "0" * 64
            bad = dict(src, **{"provenance.json": json.dumps(prov).encode()})
            copy_tables.gh = lambda *x, raw=False: fake_gh(*x, raw=raw, served=bad)
            try:
                copy_tables.main()
                check(False, "a provenance file that was not copied gives a clear error")
            except SystemExit as e:
                check("extra.parquet, which was not copied" in str(e), f"clear error: {e}"[:160])
        finally:
            sys.argv = argv

    # A build with pile-up tables but no pileup_commit is refused.
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        shutil.copytree(a.data, tmp / "data", symlinks=False)
        s = json.loads((tmp / "data/SOURCE.json").read_text())
        s["pileup_commit"] = None
        (tmp / "data/SOURCE.json").write_text(json.dumps(s))
        p = subprocess.run([sys.executable, str(ROOT / "scripts/build_site.py"), "--data", str(tmp / "data"),
                            "--out", str(tmp / "out")], capture_output=True, text=True, cwd=ROOT)
        check(p.returncode != 0 and "pileup_commit" in p.stderr, "pile-up page refused without pileup_commit")

    # 7. "Can its weak hits be trusted?" on every card, recomputed from the weak-hit table here.
    wpath = a.data / "weak_hits" / "weak_hits.parquet"
    if wpath.exists():
        import math as _m
        w = pl.read_parquet(wpath).filter(pl.col("split") == "test")
        names = {"human": "human", "worm": "C. elegans", "yeast": "S. cerevisiae",
                 "borrelia": "B. burgdorferi"}
        bands = ["le_0.001", "0.001_1", "1_10"]

        def ci(k, n, z=1.959964):
            p = k / n
            c = (p + z * z / (2 * n)) / (1 + z * z / n)
            h = z * _m.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
            return max(0.0, c - h), min(1.0, c + h)

        def want_cell(r):
            if r is None:
                return '<span class="few">not run</span>'
            k, n = r["right_strict"], r["n_judged"]
            if n == 0:
                return '<span class="few">no judged hits (n = 0)</span>'
            lo, hi = ci(k, n)
            t = f"{100 * k / n:.1f}% ({k:,} of {n:,}; 95% CI {100 * lo:.1f}–{100 * hi:.1f}%)"
            return f'<span class="few">{t}</span>' if n < 30 else t

        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run([sys.executable, str(ROOT / "scripts/build_site.py"), "--data", str(a.data),
                            "--out", tmp], check=True, capture_output=True, cwd=ROOT)
            n_cells = 0
            for t in TOOLS:
                md_path = Path(tmp) / "tools" / f"{t}.md"
                if not md_path.exists():
                    continue
                md = md_path.read_text()
                sec = md.split("## Can its weak hits be trusted?", 1)
                check(len(sec) == 2, f"{t}: card has the weak-hit block")
                if len(sec) != 2:
                    continue
                body = sec[1].split("\n## ", 1)[0]
                wt = w.filter(pl.col("tool") == t)
                if wt.height == 0:
                    check("Not run on the dark set." in body, f"{t}: card says not run on the dark set")
                    continue
                for qset, head in (("dark", "**Dark proteins**"), ("placed", "**Placed proteins**")):
                    part = body.split(head, 1)[1]
                    for sp, name in names.items():
                        line = next((l for l in part.splitlines() if l.startswith(f"| {name} |")), None)
                        cells = [c.strip() for c in line.strip("|").split(" | ")[1:]] if line else []
                        want = []
                        for b in bands:
                            r = wt.filter((pl.col("species") == sp) & (pl.col("query_set") == qset)
                                          & (pl.col("band") == b))
                            want.append(want_cell(r.row(0, named=True) if r.height else None))
                        n_cells += len(want)
                        if cells != want:
                            check(False, f"{t} {qset} {sp}: card {cells} != {want}")
            check(True, f"weak-hit block: {n_cells} card cells recomputed from weak_hits.parquet")

    if failures:
        print(f"{len(failures)} FAILURES")
        return 1
    print("all build checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
