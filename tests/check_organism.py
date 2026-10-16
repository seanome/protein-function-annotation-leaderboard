#!/usr/bin/env python3
"""Checks of the organism pages, for five species across the tree.

    python tests/check_organism.py --site docs --data data --raw ~/data/organism-page [--live]

1. The four cards, recomputed here from the source files by a second route:
   - nearest reviewed relative: every species with >= 100 reviewed entries is scored by the
     depth of its last common ancestor with the species (brute force over the NCBI taxdump),
     not by the build's walk up the lineage;
   - order and class counts: Swiss-Prot entries per taxon, summed over every taxon whose
     lineage holds that order or class (per taxon id, not per species as the build does);
   - AlphaFold share and pLDDT bands: data/organism/alphafold_counts.parquet;
   - Earth BioGenome Project: data/organism/ebp_species.parquet.
2. The verdict's numbers, read straight from data/results.parquet (the long table), not from
   the leaderboard's pivoted copy: Fmax, recall, precision, median IoU and n for each tool,
   that the setting is the tool's best on all domains at that relative, and that no other
   tool's best setting scores higher in the stratum.
3. With --live: the same species built in headless Chrome from live lookups
   (organism/live.html#taxid=…&live=1) must give the same relative, counts, AlphaFold
   numbers and verdict as the pre-built page. The distance may differ in method: the
   pre-built page asks TimeTree for the pair, the browser uses TimeTree's age of their
   shared ancestor (TimeTree's API does not allow calls from a web page).
Writes tests/organism_check.json with every comparison.
"""

import argparse
import functools
import http.server
import json
import socketserver
import sys
import threading
from pathlib import Path

import polars as pl

CHECK = {29159: "bivalve: Pacific oyster", 46731: "coral: cauliflower coral", 2850: "diatom: Phaeodactylum tricornutum",
         5141: "fungus: Neurospora crassa", 59729: "bird: zebra finch"}
BANDS = ["0-50", "50-70", "70-90", "90-100"]
TRUTH = "Pfam-A domain instances"
EXP = "qfo_region_midi_plus"


def close(a, b, tol=1e-9):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


class Report:
    def __init__(self):
        self.rows, self.bad = [], 0

    def add(self, species, what, page, source, ok):
        self.rows.append({"species": species, "check": what, "page": page, "source": source, "ok": bool(ok)})
        if not ok:
            self.bad += 1
            print(f"FAIL {species}: {what}: page {page!r} vs source {source!r}")


def taxonomy(raw):
    n = pl.read_parquet(raw / "taxonomy_nodes.parquet")
    parent = dict(zip(n["taxid"].to_list(), n["parent"].to_list()))
    rank = dict(zip(n["taxid"].to_list(), n["rank"].to_list()))
    m = pl.read_parquet(raw / "taxonomy_merged.parquet")
    merged = dict(zip(m["old"].to_list(), m["new"].to_list()))

    def lineage(t):
        t = merged.get(t, t)
        out = []
        while t in parent:
            out.append(t)
            if parent[t] == t:
                break
            t = parent[t]
        return out
    return lineage, rank


def check_cards(rep, page, raw, data, lineage, rank):
    rec, A = page["data"]["lookups"], page["data"]["answer"]
    t = rec["taxid"]
    lin = lineage(t)
    depth = {a: len(lin) - i for i, a in enumerate(lin)}  # deeper ancestor = larger number
    sp = pl.read_parquet(data / "organism/swissprot_per_species.parquet")
    best = None
    for s, e in sp.filter(pl.col("entries") >= 100).iter_rows():
        if s == t:
            continue
        lca = next((a for a in lineage(s) if a in depth), None)
        key = (depth.get(lca, 0), e, -s)
        if best is None or key > best[0]:
            best = (key, s, e, lca)
    rel = A["cards"]["relative"]
    rep.add(t, "nearest reviewed relative", rel["taxid"], best[1], rel["taxid"] == best[1])
    rep.add(t, "its reviewed entries", rel["entries"], best[2], rel["entries"] == best[2])
    rep.add(t, "shared ancestor", rel["shared_ancestor"]["taxid"], best[3], rel["shared_ancestor"]["taxid"] == best[3])
    per_taxid = pl.read_parquet(data / "organism/swissprot_per_taxid.parquet")
    for r in ("order", "class"):
        node = next((a for a in lin if rank.get(a) == r), None)
        tot = sum(e for x, e in per_taxid.iter_rows() if node in set(lineage(x))) if node else None
        got = A["cards"]["reviewed"][r]
        rep.add(t, f"reviewed entries in the {r}", got and got["entries"], tot, (got and got["entries"]) == tot)
    self_tot = sum(e for x, e in per_taxid.iter_rows() if t in set(lineage(x)))
    rep.add(t, "reviewed entries of the species", A["cards"]["reviewed"]["self"], self_tot, A["cards"]["reviewed"]["self"] == self_tot)
    af = pl.read_parquet(data / "organism/alphafold_counts.parquet").filter(pl.col("upid") == rec["proteome"]["upid"])
    c = A["cards"]["alphafold"]
    if af.height:
        r = af.row(0, named=True)
        rep.add(t, "proteome entries with an AlphaFold model", c["with_model"], r["proteome_entries_with_afdb"], c["with_model"] == r["proteome_entries_with_afdb"])
        rep.add(t, "proteome entries", c["proteome_entries"], r["proteome_entries"], c["proteome_entries"] == r["proteome_entries"])
        bands = {k: v for k, v in r.items() if k.startswith("plddt_")}
        tot = sum(bands.values())
        mean = sum(v * (int(k.split("_")[1]) + int(k.split("_")[2])) / 2 for k, v in bands.items()) / tot if tot else None
        rep.add(t, "mean pLDDT estimate", c["mean_plddt_est"], mean, close(c["mean_plddt_est"], mean))
    else:
        rep.add(t, "AlphaFold counts present", c, None, c is None)
    ebp = pl.read_parquet(data / "organism/ebp_species.parquet").filter(pl.col("taxid") == t)
    rep.add(t, "on the Earth BioGenome list", A["cards"]["ebp"]["on_list"], ebp.height == 1, A["cards"]["ebp"]["on_list"] == (ebp.height == 1))


def check_verdict(rep, page, data):
    A = page["data"]["answer"]
    t = A["taxid"]
    st = A["stratum"]
    if not st.get("relative"):
        rep.add(t, "verdict present", None, None, True)
        return
    rel = st["relative"]["id"]
    long = pl.read_parquet(data / "results.parquet").filter(
        (pl.col("experiment") == EXP) & (pl.col("truth_source") == TRUTH) & (pl.col("dedup_or_raw") == "raw")
        & (pl.col("target_db") == f"QfO {rel}") & (pl.col("tool") != "hmmscan"))
    allrows = long.filter(pl.col("stratum_axis") == "all").pivot(on="metric", index=["tool", "setting"], values="value")
    best = {}
    for r in allrows.sort(["fmax", "setting"], descending=[True, False]).iter_rows(named=True):
        best.setdefault(r["tool"], r["setting"])
    band = st["band"]
    q = long.filter(pl.col("stratum_axis") == "all") if band == "all" else \
        long.filter((pl.col("stratum_axis") == "plddt") & (pl.col("stratum_bin") == band))
    wide = q.pivot(on="metric", index=["tool", "setting", "n_instances"], values="value")
    at = {(r["tool"], r["setting"]): r for r in wide.iter_rows(named=True)}
    ranked = sorted(((at[(tl, s)]["fmax"], tl) for tl, s in best.items() if (tl, s) in at), reverse=True)
    rep.add(t, "the three tools, in order", [v["tool"] for v in A["verdict"]], [tl for _, tl in ranked[:3]],
            [v["tool"] for v in A["verdict"]] == [tl for _, tl in ranked[:3]])
    for v in A["verdict"]:
        rep.add(t, f"{v['tool']} setting is its best at {rel}", v["setting"], best[v["tool"]], v["setting"] == best[v["tool"]])
        r = at[(v["tool"], v["setting"])]
        for k_page, k_src in (("fmax", "fmax"), ("recall", "fmax_recall"), ("precision", "fmax_precision"),
                              ("overlap_median", "median_iou_tp")):
            rep.add(t, f"{v['tool']} {k_src} at {rel}, pLDDT {band}", v[k_page], r[k_src], close(v[k_page], r[k_src]))
        rep.add(t, f"{v['tool']} n at {rel}, pLDDT {band}", v["n"], r["n_instances"], v["n"] == r["n_instances"])
    if st.get("band_asked") and st["band"] != st["band_asked"]:
        n_asked = long.filter((pl.col("stratum_axis") == "plddt") & (pl.col("stratum_bin") == st["band_asked"]))["n_instances"].max()
        rep.add(t, "fallback because the asked stratum has < 200", st.get("n_asked"), n_asked or 0,
                (st.get("n_asked") or 0) == (n_asked or 0) and (n_asked or 0) < 200)


def check_live(rep, site, pages):
    from playwright.sync_api import sync_playwright

    class Quiet(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        for t, page in pages.items():
            # A fresh page each time: going to another #taxid on the same page only changes
            # the hash, and the previous species' answer would still be there.
            pg = b.new_page()
            pg.goto(base + f"organism/live.html#taxid={t}&live=1")
            pg.wait_for_function("window.__organismAnswer !== undefined", timeout=180000)
            got = pg.evaluate("window.__organismAnswer")
            if "error" in got:
                rep.add(t, "live page built", got["error"], None, False)
                continue
            L, A = got["L"], got["A"]
            P = page["data"]["answer"]
            for what, a, b_ in (
                ("live: relative", A["cards"]["relative"]["taxid"], P["cards"]["relative"]["taxid"]),
                ("live: order entries", A["cards"]["reviewed"]["order"], P["cards"]["reviewed"]["order"]),
                ("live: class entries", A["cards"]["reviewed"]["class"], P["cards"]["reviewed"]["class"]),
                ("live: proteome", A["proteome"]["upid"] if A["proteome"] else None, P["proteome"]["upid"]),
                ("live: AlphaFold share", (A["cards"]["alphafold"] or {}).get("share_with_model"), (P["cards"]["alphafold"] or {}).get("share_with_model")),
                ("live: mean pLDDT estimate", (A["cards"]["alphafold"] or {}).get("mean_plddt_est"), (P["cards"]["alphafold"] or {}).get("mean_plddt_est")),
                ("live: EBP list", (A["cards"]["ebp"] or {}).get("on_list"), P["cards"]["ebp"]["on_list"]),
            ) + ((("live: verdict tools (same measured relative)", [v["tool"] for v in A["verdict"]], [v["tool"] for v in P["verdict"]]),)
                 if (A["stratum"].get("relative") or {}).get("id") == (P["stratum"].get("relative") or {}).get("id") else ()):
                ok = close(a, b_, 1e-6) if isinstance(a, float) or isinstance(b_, float) else a == b_
                rep.add(t, what, a, b_, ok)
            rep.rows.append({"species": t, "check": "live: distance (not compared; method differs)",
                             "page": [A["cards"]["relative"]["mya"], A["cards"]["relative"]["mya_method"],
                                      (A["stratum"].get("relative") or {}).get("id"), [v["tool"] for v in A["verdict"]]],
                             "source": [P["cards"]["relative"]["mya"], P["cards"]["relative"]["mya_method"],
                                        (P["stratum"].get("relative") or {}).get("id"), [v["tool"] for v in P["verdict"]]],
                             "ok": True})
        b.close()
    httpd.shutdown()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=Path("docs"))
    ap.add_argument("--data", type=Path, default=Path("data"))
    ap.add_argument("--raw", type=Path, default=Path.home() / "data/organism-page")
    ap.add_argument("--live", action="store_true")
    a = ap.parse_args()
    prebuilt = json.loads((a.site / "data/organism/prebuilt.json").read_text())
    lineage, rank = taxonomy(a.raw)
    rep = Report()
    pages = {}
    for t, what in CHECK.items():
        slug = prebuilt.get(str(t))
        if not slug:
            rep.add(t, f"pre-built page for {what}", None, "expected", False)
            continue
        page = json.loads((a.site / f"organism/{slug}.html.json").read_text())
        pages[t] = page
        check_cards(rep, page, a.raw, a.data, lineage, rank)
        check_verdict(rep, page, a.data)
    if a.live:
        check_live(rep, a.site, pages)
    out = Path(__file__).with_name("organism_check.json")
    out.write_text(json.dumps(rep.rows, indent=1, default=str) + "\n")
    print(f"{len(rep.rows)} comparisons, {rep.bad} failed; written to {out}")
    return 1 if rep.bad else 0


if __name__ == "__main__":
    sys.exit(main())
