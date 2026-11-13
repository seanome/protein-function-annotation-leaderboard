"""Second-agent check: group medians, within-group Spearman, the three-sentence rule.

Reads the shares from check_raw.tsv (recounted from the raw tables by check_shares.py) and
the group, divergence time and BUSCO value of each proteome from verdict.json. Runs the same
numbers from verdict.json's own shares, so the two can be compared.

Usage: python check_rule.py   (from this folder; 2025-kmerseek-analysis env)
"""
import csv
import json
import statistics

from scipy.stats import spearmanr

GROUPS = ["model organism", "atlas species", "random control"]
GAP = 5.0
BUSCO_MIN = 90.0


def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(a, b):
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
    return num / den


def evaluate(rows, share_key):
    out = {}
    for label, keep in [("all", lambda p: True),
                        ("busco_90", lambda p: p["busco_completeness_pct"] >= BUSCO_MIN)]:
        med, n = {}, {}
        for g in GROUPS:
            vals = [p[share_key] for p in rows if p["group"] == g and keep(p)]
            med[g] = statistics.median(vals)
            n[g] = len(vals)
        a, c, m = med["atlas species"], med["random control"], med["model organism"]
        out[label] = {
            "median": med, "n": n,
            "atlas_minus_control": a - c, "atlas_minus_model": a - m,
            "control_minus_model": c - m,
            "s1_atlas_5_above_controls": a - c >= GAP,
            "s2_within_5": abs(a - c) < GAP,
            "s2_both_5_above_models": a - m >= GAP and c - m >= GAP,
        }
    out["spearman"] = {}
    for g in GROUPS:
        sub = [p for p in rows if p["group"] == g]
        x = [p["divergence_my"] for p in sub]
        y = [p[share_key] for p in sub]
        rho_own = pearson(ranks(x), ranks(y))
        rho, pval = spearmanr(x, y)
        out["spearman"][g] = {"n": len(sub), "rho_own": round(rho_own, 3),
                              "rho_scipy": round(float(rho), 3), "p": round(float(pval), 4)}
    return out


def main():
    v = json.load(open("verdict.json"))
    rows = [p for p in v["proteomes"] if p["group"] in GROUPS]
    raw = {r["upid"]: r for r in csv.DictReader(open("check_raw.tsv"), delimiter="\t")}
    for p in rows:
        p["share_recount"] = 100.0 * int(raw[p["upid"]]["n_unlabelled"]) / int(
            raw[p["upid"]]["n_proteins"])
    result = {"from_recount": evaluate(rows, "share_recount"),
              "from_verdict_json": evaluate(rows, "unlabelled_pct"),
              "verdict_json_says": {"medians_all": v["medians_all"],
                                    "medians_busco_90": v["medians_busco_90"],
                                    "n_busco_90": v["n_busco_90"],
                                    "sentence_held": v["sentence_held"]}}
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
