#!/usr/bin/env python3
"""Describe the SCOPe40 domain sets and Foldseek's scoring, before any SCOPe40 ranking is scored.

Writes data/scope40_testbed.json, which the Methods page reads. Nothing here runs a search.

Inputs (none are in this repository; their sha256 sums go into the JSON):
  --astral208  ASTRAL SCOPe 2.08 domains at 40% identity, the FASTA every kmerseek SCOPe40
               run in seanome/2024-kmerseek-analysis searched (PR 107, notebook 277).
  --lookup201  Foldseek's own SCOPe40 list, SCOPe 2.01 (domain id, SCOPe id), as copied into
               seanome/2024-kmerseek-analysis data/foldseek-analysis/scop_lookup.fix.tsv.
  --rocx-dir   The per-query result files of Foldseek and MMseqs2 on SCOPe 2.08 that the TEA
               authors published (data/tea_scope40_rocx_files in the same repository).
  --bench-awk  Foldseek's scoring script, scopbenchmark/scripts/bench.noselfhit.awk from
               github.com/steineggerlab/foldseek-analysis.

Levels follow the Foldseek paper and its script: a family true positive is the same family;
a superfamily true positive is the same superfamily but a different family; a fold true
positive is the same fold but a different superfamily. A hit in a different fold is a false
positive at every level. The nested rule counts every same-superfamily hit, family included,
as a superfamily true positive.
"""

import argparse
import collections
import hashlib
import json
from pathlib import Path


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_fasta_ids(p: Path) -> dict[str, str]:
    out = {}
    for line in p.open():
        if line.startswith(">"):
            f = line[1:].split()
            out[f[0]] = f[1]
    return out


def read_lookup(p: Path) -> dict[str, str]:
    return dict(line.split()[:2] for line in p.open() if line.strip())


def parts(sccs: str) -> tuple[str, str, str, str]:
    """(class, fold, superfamily, family) as prefixes: a, a.1, a.1.1, a.1.1.1."""
    p = sccs.split(".")
    return p[0], ".".join(p[:2]), ".".join(p[:3]), sccs


def describe(dom: dict[str, str]) -> dict:
    """Domain and query counts for one domain set."""
    fam = collections.Counter(parts(s)[3] for s in dom.values())
    sf = collections.Counter(parts(s)[2] for s in dom.values())
    fo = collections.Counter(parts(s)[1] for s in dom.values())
    cls = collections.Counter(parts(s)[0] for s in dom.values())
    n = collections.Counter()
    for s in dom.values():
        c, f, s3, s4 = parts(s)
        n_fam = fam[s4] - 1                 # other members of the family
        n_sf_excl = sf[s3] - fam[s4]        # superfamily members outside the family
        n_sf_nested = sf[s3] - 1            # every other superfamily member
        n_fold = fo[f] - sf[s3]             # fold members outside the superfamily
        n["family"] += n_fam > 0
        n["superfamily"] += n_sf_excl > 0
        n["superfamily_nested"] += n_sf_nested > 0
        n["fold"] += n_fold > 0
        n["all_three"] += n_fam > 0 and n_sf_excl > 0 and n_fold > 0
    return {
        "n_domains": len(dom),
        "n_by_class": dict(sorted(cls.items())),
        "n_families": len(fam), "n_superfamilies": len(sf), "n_folds": len(fo),
        "n_queries_with_a_true_positive": {k: n[k] for k in
                                           ("family", "superfamily", "superfamily_nested", "fold", "all_three")},
    }


def all_three_ids(dom: dict[str, str]) -> set[str]:
    fam = collections.Counter(parts(s)[3] for s in dom.values())
    sf = collections.Counter(parts(s)[2] for s in dom.values())
    fo = collections.Counter(parts(s)[1] for s in dom.values())
    return {d for d, s in dom.items()
            if fam[parts(s)[3]] > 1 and sf[parts(s)[2]] > fam[parts(s)[3]] and fo[parts(s)[1]] > sf[parts(s)[2]]}


def read_rocx(p: Path) -> list[dict]:
    lines = p.read_text().splitlines()
    sep = "," if "," in lines[0] else "\t"
    head = lines[0].split(sep)
    return [dict(zip(head, line.split(sep))) for line in lines[1:]]


def denominators(r: dict) -> dict:
    """Foldseek's script denominators and the number of possible true positives, per level."""
    f, s, o = int(r["FAMCNT"]), int(r["SFAMCNT"]), int(r["FOLDCNT"])
    return {"FAM": (f - 1, f - 1), "SFAM": (s - (f - 1), s - f), "FOLD": (o - (s - 1), o - s)}


def rocx_summary(p: Path) -> dict:
    """Mean sensitivity per level, as the file has it and with the possible-TP denominator."""
    rows = read_rocx(p)
    out = {"file": p.name, "n_queries": len(rows), "levels": {}}
    for lvl in ("FAM", "SFAM", "FOLD"):
        whole_script = whole_possible = 0
        a = b = 0.0
        for r in rows:
            v = float(r[lvl])
            d_script, d_possible = denominators(r)[lvl]
            tp = v * d_script
            whole_script += abs(tp - round(tp)) < 1e-3
            whole_possible += abs(v * d_possible - round(v * d_possible)) < 1e-3
            a += v
            b += round(tp) / d_possible
        out["levels"][lvl] = {
            "mean_as_published": a / len(rows),
            "mean_over_possible_true_positives": b / len(rows),
            "rows_whole_with_script_denominator": whole_script,
            "rows_whole_with_possible_denominator": whole_possible,
        }
    return out


def counts_vs(rows: list[dict], dom: dict[str, str]) -> dict:
    """Rows whose family, superfamily or fold size differs from the size in the domain set."""
    fam = collections.Counter(parts(s)[3] for s in dom.values())
    sf = collections.Counter(parts(s)[2] for s in dom.values())
    fo = collections.Counter(parts(s)[1] for s in dom.values())
    out = {"n_rows": len(rows), "FAMCNT": 0, "SFAMCNT": 0, "FOLDCNT": 0, "any": 0, "smaller_only": True}
    for r in rows:
        c = parts(r["SCOP"])
        diff = False
        for col, n in (("FAMCNT", fam[c[3]]), ("SFAMCNT", sf[c[2]]), ("FOLDCNT", fo[c[1]])):
            if int(r[col]) != n:
                out[col] += 1
                diff = True
                out["smaller_only"] &= int(r[col]) < n
        out["any"] += diff
    return out


def example_row(p: Path, name: str) -> dict:
    r = next(x for x in read_rocx(p) if x["NAME"] == name)
    d = denominators(r)
    return {"file": p.name, "name": name, "scop": r["SCOP"], "SFAM": float(r["SFAM"]),
            "FAMCNT": int(r["FAMCNT"]), "SFAMCNT": int(r["SFAMCNT"]),
            "sfam_denominator_script": d["SFAM"][0], "sfam_possible": d["SFAM"][1],
            "sfam_tp": round(float(r["SFAM"]) * d["SFAM"][0])}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--astral208", type=Path, required=True)
    ap.add_argument("--lookup201", type=Path, required=True)
    ap.add_argument("--rocx-dir", type=Path, required=True)
    ap.add_argument("--bench-awk", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("data/scope40_testbed.json"))
    a = ap.parse_args()

    d208 = read_fasta_ids(a.astral208)
    d201 = read_lookup(a.lookup201)
    tea = read_rocx(a.rocx_dir / "eba.rocx")
    fsk = read_rocx(a.rocx_dir / "foldseek.rocx")
    all3 = all_three_ids(d208)
    out = {
        "inputs": {
            "astral208": {"path": "~/data/scope/" + a.astral208.name,
                          "sherlock_path": "$SCRATCH/scope/" + a.astral208.name,
                          "sha256": sha256(a.astral208)},
            "lookup201": {"path": "seanome/2024-kmerseek-analysis data/foldseek-analysis/" + a.lookup201.name,
                          "sha256": sha256(a.lookup201)},
            "rocx_dir": {"path": "seanome/2024-kmerseek-analysis data/tea_scope40_rocx_files",
                         "sha256": {p.name: sha256(p) for p in sorted(a.rocx_dir.glob("*.rocx"))}},
            "bench_awk": {"path": "github.com/steineggerlab/foldseek-analysis scopbenchmark/scripts/"
                          + a.bench_awk.name, "sha256": sha256(a.bench_awk)},
        },
        "scope208_astral40": describe(d208),
        "scope201_foldseek": describe(d201),
        "tea_query_set": {
            "n_eba": len(tea), "eba_equals_all_three": {x["NAME"] for x in tea} == all3,
            "n_foldseek_rocx": len(fsk),
            "foldseek_rocx_subset_of_all_three": {x["NAME"] for x in fsk} <= all3,
            "missing_from_foldseek_rocx": sorted(all3 - {x["NAME"] for x in fsk}),
        },
        "rocx": [rocx_summary(a.rocx_dir / f) for f in ("foldseek.rocx", "mmseqs.rocx")],
        "rocx_counts_vs_208": counts_vs(fsk, d208),
        "example": example_row(a.rocx_dir / "foldseek.rocx", "d6jc0b_"),
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
