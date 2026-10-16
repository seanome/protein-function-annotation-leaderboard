#!/usr/bin/env python3
"""Check Figures 11 to 14 against results.parquet, read straight from the long table.

Every hover text on Figures 11 (Fmax against Mya), 12 (Fmax by identity) and 13 (Fmax by
feature type) names a tool, a setting, a subset and an Fmax with its n; each is looked up
in data/results.parquet, without going through build_site.py, and must match to the 3
decimals shown. Figure 14 on each tool card is checked the same way for precision and
recall. The count of marks must equal the count of rows the table holds for those filters.
"""

import argparse
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from tools import SPECIES, TOOLS  # noqa: E402

EXPERIMENT = "qfo_region_midi_plus"
NAME2TOOL = {v["name"]: k for k, v in TOOLS.items()}
ORG2SP = {v: k for k, v in SPECIES.items()}


class Tips(HTMLParser):
    """Collect data-tip texts per figure id."""

    def __init__(self):
        super().__init__()
        self.fig, self.tips = None, {}

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "figure" and "id" in a:
            self.fig = a["id"]
        if "data-tip" in a and self.fig:
            self.tips.setdefault(self.fig, []).append(a["data-tip"])

    def handle_endtag(self, tag):
        if tag == "figure":
            self.fig = None


def tips(path: Path) -> dict:
    p = Tips()
    p.feed(path.read_text())
    return p.tips


def split_name(s: str) -> tuple[str, str]:
    """'MMseqs2, three rounds, sensitivity 7, 3 rounds: ...' -> (tool key, rest after 'name, ').
    The longest tool name that starts the text wins, so 'MMseqs2' never takes the
    three-round tool's marks."""
    for name in sorted(NAME2TOOL, key=len, reverse=True):
        if s.startswith(name + ', '):
            return NAME2TOOL[name], s[len(name) + 2:]
    raise SystemExit(f'no tool name starts: {s}')


def label_to_setting(tool: str, label: str, settings: list[str]) -> str:
    if tool == "kmerseek":
        m = re.fullmatch(r"(\S+), k = (\d+)", label)
        return f"{m[1]}_k{m[2]}_lcFalse"
    assert len(settings) == 1, (tool, settings)
    return settings[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=Path("docs"))
    ap.add_argument("--data", type=Path, default=Path("data"))
    a = ap.parse_args()
    long = pl.read_parquet(a.data / "results.parquet").filter(pl.col("experiment") == EXPERIMENT)
    settings = {t: sorted(long.filter(pl.col("tool") == t)["setting"].unique().to_list()) for t in TOOLS
                if long.filter(pl.col("tool") == t).height}
    errors, n_checked = [], 0

    def look(tool, setting, truth, target, axis, bin_, metric):
        q = long.filter((pl.col("tool") == tool) & (pl.col("setting") == setting)
                        & (pl.col("truth_source") == truth) & (pl.col("dedup_or_raw") == "raw")
                        & (pl.col("target_db") == target) & (pl.col("stratum_axis") == axis)
                        & (pl.col("stratum_bin") == bin_) & (pl.col("metric") == metric))
        if q.height != 1:
            raise SystemExit(f"{q.height} rows for {tool} {setting} {target} {axis} {bin_} {metric}")
        return q.row(0, named=True)

    t = tips(a.site / "index.html")
    # Figure 11
    rx = re.compile(r"(.+?): (.+) \(([\d,]+) Mya\), Fmax ([\d.]+), n = ([\d,]+) true domains")
    for s in t.get("fig-distance", []):
        tool, rest = split_name(s)
        g = rx.fullmatch(rest)
        if not g:
            errors.append(f"Figure 11 tip not parsed: {s}")
            continue
        m = [None, tool] + list(g.groups())  # m[2] label, m[3] relative, m[4] Mya, m[5] Fmax, m[6] n
        setting = label_to_setting(tool, m[2], settings[tool])
        r = look(tool, setting, "Pfam-A domain instances", "QfO " + ORG2SP[m[3]], "all", "all", "fmax")
        n_checked += 1
        if f"{r['value']:.3f}" != m[5] or r["n_instances"] != int(m[6].replace(",", "")) \
                or int(r["target_mya_to_human"]) != int(m[4].replace(",", "")):
            errors.append(f"Figure 11: {s} but table has {r['value']:.3f}, n {r['n_instances']}")
    full = long.filter((pl.col("truth_source") == "Pfam-A domain instances") & (pl.col("dedup_or_raw") == "raw")
                       & (pl.col("stratum_axis") == "all") & (pl.col("metric") == "fmax")
                       & pl.col("target_db").str.starts_with("QfO ") & (pl.col("tool") != "kmerseek"))
    n_marks_other = sum(1 for s in t.get("fig-distance", []) if not s.startswith("kmerseek"))
    if n_marks_other != full.height:
        errors.append(f"Figure 11: {n_marks_other} marks for the other tools, table has {full.height} rows")
    # Figure 12
    rx = re.compile(r"(.+?): (.+?), Fmax ([\d.]+), n = ([\d,]+) true domains")
    back = {"0–20%": "0-20%", "20–30%": "20-30%", "30–40%": "30-40%", "40–60%": "40-60%", "60–100%": "60-100%",
            "no same-family domain": "no_homolog"}
    for s in t.get("fig-identity", []):
        tool, rest = split_name(s)
        m = [None, tool] + list(rx.fullmatch(rest).groups())  # m[2] label, m[3] bin, m[4] Fmax, m[5] n
        setting = label_to_setting(tool, m[2], settings[tool])
        r = look(tool, setting, "Pfam-A domain instances", "QfO zebrafish", "identity", back[m[3]], "fmax")
        n_checked += 1
        if f"{r['value']:.3f}" != m[4] or r["n_instances"] != int(m[5].replace(",", "")):
            errors.append(f"Figure 12: {s} but table has {r['value']:.3f}, n {r['n_instances']}")
    # Figure 13: bars carry the tip on the label and on the bar, so each pair is one bar
    rx = re.compile(r"(.+?): Fmax ([\d.]+), n = ([\d,]+) true (.+) features")
    keys = {"membrane-spanning segment": "TRANSMEM", "zinc finger": "ZN_FING", "coiled coil": "COILED",
            "repeat": "REPEAT", "domain": "DOMAIN"}
    for s in t.get("fig-features", [])[1::2]:
        m = rx.fullmatch(s)
        name = m[1]
        if name.startswith("kmerseek, "):
            tool, setting = "kmerseek", label_to_setting("kmerseek", name[len("kmerseek, "):], [])
        else:
            tool = NAME2TOOL[name]
            setting = settings[tool][0]
        r = look(tool, setting, "Swiss-Prot features", "QfO zebrafish", "feature_type", keys[m[4]], "fmax")
        n_checked += 1
        if f"{r['value']:.3f}" != m[2] or r["n_instances"] != int(m[3].replace(",", "")):
            errors.append(f"Figure 13: {s} but table has {r['value']:.3f}")
        if tool == "kmerseek":
            best = long.filter((pl.col("tool") == "kmerseek") & (pl.col("truth_source") == "Swiss-Prot features")
                               & (pl.col("dedup_or_raw") == "raw") & (pl.col("target_db") == "QfO zebrafish")
                               & (pl.col("stratum_axis") == "feature_type") & (pl.col("stratum_bin") == keys[m[4]])
                               & (pl.col("metric") == "fmax"))["value"].max()
            if abs(best - r["value"]) > 1e-12:
                errors.append(f"Figure 13: kmerseek bar for {m[4]} is not its best setting")
    # Figure 14 on each card
    rx = re.compile(r"(.+?): precision ([\d.]+), recall ([\d.]+), Fmax ([\d.]+), n = ([\d,]+) true domains")
    for card in sorted((a.site / "tools").glob("*.html")):
        if card.name == "index.html":
            continue
        ct = tips(card).get("fig-pr", [])
        tool = card.stem
        for s in ct:
            m = rx.fullmatch(s)
            if not m:
                errors.append(f"{card.name} Figure 14 tip not parsed: {s}")
                continue
            q = long.filter((pl.col("tool") == tool) & (pl.col("truth_source") == "Pfam-A domain instances")
                            & (pl.col("dedup_or_raw") == "raw") & (pl.col("metric").is_in(["fmax_precision", "fmax_recall", "fmax"])))
            hits = q.pivot(on="metric", index=["setting", "target_db", "stratum_axis", "stratum_bin", "n_instances"],
                           values="value").filter(
                (pl.col("fmax_precision").round(3) == float(m[2])) & (pl.col("fmax_recall").round(3) == float(m[3]))
                & (pl.col("fmax").round(3) == float(m[4])) & (pl.col("n_instances") == int(m[5].replace(",", ""))))
            n_checked += 1
            if hits.height == 0:
                errors.append(f"{card.name} Figure 14: no table row with {s}")
    # Figures 15 and 17, recomputed here from the pile-up tables with their own code.
    if (a.data / "pileup/calls.parquet").exists():
        c = pl.read_parquet(a.data / "pileup/calls.parquet")
        reg = pl.read_parquet(a.data / "pileup/regions.parquet")
        c = c.join(reg.select("row_id", "accession", "true_start", "true_end"), on="row_id")
        ov = (c["call_end"].clip(upper_bound=None).to_numpy())
        import numpy as np
        s_ = np.maximum(c["call_start"].to_numpy(), c["true_start"].to_numpy())
        e_ = np.minimum(ov, c["true_end"].to_numpy())
        inter = np.clip(e_ - s_ + 1, 0, None)
        union = np.maximum(ov, c["true_end"].to_numpy()) - np.minimum(c["call_start"].to_numpy(), c["true_start"].to_numpy()) + 1
        c = c.with_columns(pl.Series("hit", c["is_tp"].to_numpy() & (inter / union >= 0.5)))
        per_set = c.filter("hit").group_by("tool", "setting").agg(pl.col("row_id").n_unique().alias("n"))
        ks = per_set.filter(pl.col("tool") == "kmerseek").sort(["n", "setting"], descending=[True, False])["setting"][0]
        use = c.filter((pl.col("tool") != "hmmscan") & ((pl.col("tool") != "kmerseek") | (pl.col("setting") == ks)))
        found = use.filter("hit").group_by("row_id").agg(pl.col("tool").unique().sort().alias("tools"))
        combo = reg.select("row_id").join(found, on="row_id", how="left").with_columns(
            pl.col("tools").fill_null([]).list.join("+").alias("k")).group_by("k").agg(pl.len())
        want = sorted(combo["len"].to_list(), reverse=True)
        got = []
        for s in t.get("fig-upset", []):
            m = re.match(r"([\d,]+) regions (found only by|found by [^.]* and no other tool|found by none)", s)
            if m:
                got.append(int(m[1].replace(",", "")))
        n_checked += len(got)
        # every drawn column must be one of the recomputed combination sizes, the largest ones
        if not got or sorted(got, reverse=True)[0] != want[0] or any(g not in want for g in got):
            errors.append(f"Figure 15 columns {sorted(got, reverse=True)[:5]} vs recomputed {want[:5]}")
        n_none = reg.height - found.height
        if f"{n_none:,} regions found by none" not in " ".join(t.get("fig-upset", [])):
            errors.append(f"Figure 15: the none column is not {n_none:,}")
        kinds = {tl: TOOLS[tl]["kind"] for tl in TOOLS}
        pk = use.filter("hit").with_columns(pl.col("tool").replace_strict(kinds).alias("kind")) \
            .group_by("accession").agg(pl.col("kind").n_unique().alias("k"))
        prot = reg.select("accession").unique().join(pk, on="accession", how="left").with_columns(pl.col("k").fill_null(0))
        sent = (a.site / "proteins.html").read_text()
        for label, val in (("none found", prot.filter(pl.col("k") == 0).height),
                           ("found by all four kinds", prot.filter(pl.col("k") == 4).height)):
            pat = (f"{val:,} had none found" if label == "none found" else f"{val:,} had a region found by all four kinds")
            n_checked += 1
            if pat not in sent:
                errors.append(f"Figure 17: expected '{pat}'")
        if len(tips(a.site / "proteins.html").get("fig-map", [])) != prot.height:
            errors.append("Figure 17: tile count differs from the protein count")
    if not t.get("fig-distance") or not t.get("fig-identity") or not t.get("fig-features"):
        errors.append("a figure on index.html has no hover marks")
    for e in errors:
        print("FAIL", e)
    print(f"figures: {n_checked} hover numbers checked against results.parquet")
    if errors:
        sys.exit(1)
    print("figure checks passed")


if __name__ == "__main__":
    main()
