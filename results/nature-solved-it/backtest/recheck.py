"""Independent recomputation of the transferability back-test statistics.

Reads scored.csv (and key.csv + blind_scores.csv + sets.csv as a second route),
recomputes n, median T, Mann-Whitney U, p values and AUC, and compares with
results.json. Written by the checking agent; it does not import analyze.py.
"""
import csv, json, itertools
from pathlib import Path
import numpy as np
from scipy.stats import mannwhitneyu

HERE = Path("/Users/olga/code/pfal-transferability-backtest/results/nature-solved-it/backtest")


def rd(name):
    with (HERE / name).open() as f:
        return list(csv.DictReader(f))


scored = rd("scored.csv")
key = {r["blind_id"]: r for r in rd("key.csv")}
blind = {r["blind_id"]: r for r in rd("blind_scores.csv")}
sets = {r["entry_id"]: r for r in rd("sets.csv")}
res = json.load((HERE / "results.json").open())

# Route 2: rebuild T from blind_scores.csv joined through key.csv, and check scored.csv.
mismatch = []
for r in scored:
    b = blind[r["blind_id"]]
    k = key[r["blind_id"]]
    pts = [int(b["T%d" % i]) for i in range(1, 7)]
    if sum(pts) != int(r["T"]) or int(b["T"]) != int(r["T"]):
        mismatch.append((r["blind_id"], sum(pts), b["T"], r["T"]))
    if k["entry_id"] != r["entry_id"] or k["set"] != r["set"] or sets[r["entry_id"]]["set"] != r["set"]:
        mismatch.append((r["blind_id"], "join", k["entry_id"], r["entry_id"]))
print("rows in scored.csv:", len(scored))
print("join / sum-of-points mismatches:", mismatch if mismatch else "none")

A = [int(r["T"]) for r in scored if r["set"] == "A"]
Bi = [int(r["T"]) for r in scored if r["set"] == "B_i"]
Bii = [int(r["T"]) for r in scored if r["set"] == "B_ii"]
B = Bi + Bii

print("n:", dict(A=len(A), B=len(B), B_i=len(Bi), B_ii=len(Bii)), "| results.json:", res["n"])
med = dict(A=float(np.median(A)), B=float(np.median(B)), B_i=float(np.median(Bi)), B_ii=float(np.median(Bii)))
print("median T:", med, "| results.json:", res["median_T"])
print("mean T:", {k: round(float(np.mean(v)), 3) for k, v in dict(A=A, B=B, B_i=Bi, B_ii=Bii).items()},
      "| results.json:", res["mean_T"])

U, p2 = mannwhitneyu(A, B, alternative="two-sided")
U1, p1 = mannwhitneyu(A, B, alternative="greater")
print("U (A over B):", U, "| results.json:", res["mann_whitney"]["U_A_over_B"])
print("p two-sided:", p2, "| results.json:", res["mann_whitney"]["p_two_sided"])
print("p one-sided A greater:", p1, "| results.json:", res["mann_whitney"]["p_one_sided_A_greater"])

greater = sum(1 for a, b in itertools.product(A, B) if a > b)
tied = sum(1 for a, b in itertools.product(A, B) if a == b)
pairs = len(A) * len(B)
auc_u = U / pairs
auc_pairs = (greater + 0.5 * tied) / pairs
print("pairs A>B / tied / total:", greater, tied, pairs,
      "| results.json:", res["mann_whitney"]["pairs_A_greater"], res["mann_whitney"]["pairs_tied"], res["mann_whitney"]["pairs_total"])
print("AUC U/(nA*nB):", round(auc_u, 4), "| by counting pairs:", round(auc_pairs, 4),
      "| results.json:", res["auc"])

# The three conditions of the pass/fail rule in INCLUSION.md.
c1 = (med["A"] - med["B"]) >= 2
c2 = auc_u >= 0.75
c3 = len(A) >= 15 and len(B) >= 15
print("rule: median gap >= 2:", c1, "(gap %.1f)" % (med["A"] - med["B"]),
      "| AUC >= 0.75:", c2, "| >=15 per set:", c3, "=> T predicts:", c1 and c2 and c3)
print("results.json rule:", res["rule"])

# Split comparisons.
for name, other in (("A_vs_B_i", Bi), ("A_vs_B_ii", Bii)):
    u, pp = mannwhitneyu(A, other, alternative="two-sided")
    u1, pp1 = mannwhitneyu(A, other, alternative="greater")
    print(name, "U", u, "auc", round(u / (len(A) * len(other)), 4), "p2", pp, "p1", pp1,
          "| results.json:", res["split_comparisons"][name])

# Source-class share of Set A, under the definition in analyze.py.
external = {"A01", "A02", "A03", "A04", "A05", "A06"}
internal = {"A10"}
a_ids = [r["entry_id"] for r in scored if r["set"] == "A"]
print("external numerator/denominator:", len(external & set(a_ids)), "/", len(a_ids),
      "| incl internal:", len((external | internal) & set(a_ids)), "/", len(a_ids),
      "| results.json:", res["venom_saliva_secretion_share_of_A"]["numerator"],
      res["venom_saliva_secretion_share_of_A"]["denominator"])
t_ext = [int(r["T"]) for r in scored if r["entry_id"] in external]
t_rest = [int(r["T"]) for r in scored if r["set"] == "A" and r["entry_id"] not in external]
print("median T external:", float(np.median(t_ext)), "rest of A:", float(np.median(t_rest)),
      "| results.json:", res["venom_saliva_secretion_share_of_A"]["median_T_external"],
      res["venom_saliva_secretion_share_of_A"]["median_T_rest_of_A"])

# Sensitivity: no record failed verification, so nothing is dropped. These two rows are
# the ones whose molecule is approved by an agency other than the FDA or the EMA.
for drop in (["B111"], ["B114"], ["B111", "B114"], ["B115"]):
    keep_B = [int(r["T"]) for r in scored if r["set"].startswith("B") and r["entry_id"] not in drop]
    keep_A = [int(r["T"]) for r in scored if r["set"] == "A" and r["entry_id"] not in drop]
    u, pp = mannwhitneyu(keep_A, keep_B, alternative="two-sided")
    gap = float(np.median(keep_A)) - float(np.median(keep_B))
    a_ = u / (len(keep_A) * len(keep_B))
    print("drop", drop, "nA", len(keep_A), "nB", len(keep_B), "gap %.1f" % gap,
          "AUC %.4f" % a_, "p2 %.4f" % pp,
          "=> T predicts:", gap >= 2 and a_ >= 0.75 and len(keep_A) >= 15 and len(keep_B) >= 15)
