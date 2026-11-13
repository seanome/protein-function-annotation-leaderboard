"""Rejoin the blind T scores with the outcome, run the three statistics, draw the figure.

Inputs in this directory:
  sets.csv          the two sets, collected before any score existed
  key.csv           blind id to entry id (held back from the scoring agent)
  blind_scores.csv  T1..T6 with reasons and T, from an agent that read only
                    blind_input.csv

Outputs in this directory:
  scored.csv        one row per entry: set, entry, T1..T6, T, year
  results.json      every number this script reports
  figures/transferability_t_distribution_backtest.{png,pdf,svg}

No score is changed here. The pass/fail rule is the one in INCLUSION.md and is applied
as written.

Run: /Users/olga/anaconda3/envs/2025-kmerseek-analysis/bin/python analyze.py
"""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

import pubfig as pf

HERE = Path(__file__).resolve().parent
FIGDIR = HERE / "figures"

# Which Set A entries come from a secretion the animal puts outside its body (venom,
# saliva, other oral or skin secretion). Written as a rule, then applied by entry id:
# "external" = venom, saliva or another secretion released out of the animal;
# "internal" = a molecule the animal secretes into its own blood;
# "tissue" = taken from an organ or body tissue;
# "small_molecule" = a non-protein chemical the animal makes.
SOURCE_CLASS = {
    "A01": "external", "A02": "external", "A03": "external", "A04": "external",
    "A05": "external", "A06": "external",
    "A07": "tissue", "A08": "tissue", "A09": "tissue",
    "A10": "internal",
    "A11": "tissue", "A12": "tissue", "A13": "tissue",
    "A14": "tissue", "A15": "tissue", "A16": "tissue",
    "A17": "small_molecule", "A18": "small_molecule", "A19": "tissue",
}
# Whether the molecule that carries the activity is a peptide or protein.
IS_PEPTIDE = {
    "A01": False,  # the marketed molecule is a small synthetic copy of a venom peptide
    "A02": True, "A03": True, "A04": True, "A05": True, "A06": True,
    "A07": False,  # DNA fragments
    "A08": False,  # sugar chain
    "A09": True, "A10": True, "A11": True, "A12": True,
    "A13": True,   # a mixture of enzymes
    "A14": False, "A15": False, "A16": False,  # mostly lipid
    "A17": False, "A18": False, "A19": False,
}


def read_csv(name):
    with (HERE / name).open() as f:
        return list(csv.DictReader(f))


def main():
    sets = {r["entry_id"]: r for r in read_csv("sets.csv")}
    key = {r["blind_id"]: r for r in read_csv("key.csv")}
    scores = {r["blind_id"]: r for r in read_csv("blind_scores.csv")}

    assert set(key) == set(scores), "blind ids do not match"
    assert len(sets) == len(key) == 43, (len(sets), len(key))

    rows = []
    for blind_id, k in key.items():
        s = scores[blind_id]
        e = sets[k["entry_id"]]
        parts = [int(s[f"T{i}"]) for i in range(1, 7)]
        assert sum(parts) == int(s["T"]), blind_id
        rows.append(dict(
            blind_id=blind_id, entry_id=k["entry_id"], set=e["set"],
            drug_or_programme=e["drug_or_programme"], species=e["species"],
            molecule=e["molecule"], year=int(e["year"]), year_type=e["year_type"],
            T=int(s["T"]), **{f"T{i}": parts[i - 1] for i in range(1, 7)},
            **{f"T{i}_reason": s[f"T{i}_reason"] for i in range(1, 7)}))
    rows.sort(key=lambda r: (r["set"], r["entry_id"]))

    with (HERE / "scored.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    A = np.array([r["T"] for r in rows if r["set"] == "A"])
    Bi = np.array([r["T"] for r in rows if r["set"] == "B_i"])
    Bii = np.array([r["T"] for r in rows if r["set"] == "B_ii"])
    B = np.concatenate([Bi, Bii])

    med = {k: float(np.median(v)) for k, v in
           dict(A=A, B=B, B_i=Bi, B_ii=Bii).items()}

    u_two = mannwhitneyu(A, B, alternative="two-sided")
    u_greater = mannwhitneyu(A, B, alternative="greater")
    n_a, n_b = len(A), len(B)
    # U from scipy is the count of (A, B) pairs where A > B, plus half the ties.
    auc_hand = float(u_greater.statistic) / (n_a * n_b)
    pairs_greater = int(sum(1 for a in A for b in B if a > b))
    pairs_tied = int(sum(1 for a in A for b in B if a == b))
    auc_by_pairs = (pairs_greater + 0.5 * pairs_tied) / (n_a * n_b)
    labels = np.concatenate([np.ones(n_a), np.zeros(n_b)])
    auc_sklearn = float(roc_auc_score(labels, np.concatenate([A, B])))

    # Descriptive, not part of the rule: how Set A compares with each kind of Set B
    # entry on its own. This says where the separation in the pooled comparison comes
    # from.
    def auc_against(other):
        u = mannwhitneyu(A, other, alternative="greater")
        return dict(n=len(other), median=float(np.median(other)),
                    U=float(u.statistic), auc=round(float(u.statistic)
                                                    / (n_a * len(other)), 4),
                    p_one_sided=float(u.pvalue),
                    p_two_sided=float(mannwhitneyu(A, other,
                                                   alternative="two-sided").pvalue))

    split = dict(A_vs_B_i=auc_against(Bi), A_vs_B_ii=auc_against(Bii))

    # The venom, saliva and secretion count that the narrow claim rests on.
    a_rows = [r for r in rows if r["set"] == "A"]
    n_external = sum(1 for r in a_rows if SOURCE_CLASS[r["entry_id"]] == "external")
    n_external_or_internal = sum(1 for r in a_rows
                                 if SOURCE_CLASS[r["entry_id"]] in ("external", "internal"))
    n_peptide = sum(1 for r in a_rows if IS_PEPTIDE[r["entry_id"]])
    external_ids = [r["entry_id"] for r in a_rows
                    if SOURCE_CLASS[r["entry_id"]] == "external"]
    t_external = [r["T"] for r in a_rows if SOURCE_CLASS[r["entry_id"]] == "external"]
    t_other_a = [r["T"] for r in a_rows if SOURCE_CLASS[r["entry_id"]] != "external"]

    years_a = sorted(r["year"] for r in a_rows)
    years_bi = sorted(r["year"] for r in rows if r["set"] == "B_i")
    years_bii = sorted(r["year"] for r in rows if r["set"] == "B_ii")

    # The rule from INCLUSION.md, applied as written.
    rule = dict(
        median_gap_at_least_2=bool(med["A"] - med["B"] >= 2),
        auc_at_least_0_75=bool(auc_hand >= 0.75),
        at_least_15_per_set=bool(n_a >= 15 and n_b >= 15),
    )
    rule["t_predicts"] = bool(all(rule.values()))

    results = dict(
        generated="2026-10-19",
        n=dict(A=n_a, B=n_b, B_i=len(Bi), B_ii=len(Bii)),
        median_T=med,
        mean_T={k: round(float(np.mean(v)), 3) for k, v in
                dict(A=A, B=B, B_i=Bi, B_ii=Bii).items()},
        T_counts={k: {str(t): int((v == t).sum()) for t in range(7)} for k, v in
                  dict(A=A, B=B, B_i=Bi, B_ii=Bii).items()},
        mann_whitney=dict(
            U_A_over_B=float(u_greater.statistic),
            p_two_sided=float(u_two.pvalue),
            p_one_sided_A_greater=float(u_greater.pvalue),
            pairs_A_greater=pairs_greater, pairs_tied=pairs_tied,
            pairs_total=n_a * n_b),
        split_comparisons=split,
        auc=dict(by_hand_U_over_nA_nB=round(auc_hand, 4),
                 by_counting_pairs=round(auc_by_pairs, 4),
                 sklearn_roc_auc_score=round(auc_sklearn, 4)),
        venom_saliva_secretion_share_of_A=dict(
            numerator=n_external, denominator=n_a,
            share=round(n_external / n_a, 4), entry_ids=external_ids,
            definition="the source material is a secretion the animal releases outside "
                       "its body: venom, saliva or another oral or skin secretion",
            including_internal_secretion=dict(numerator=n_external_or_internal,
                                              denominator=n_a),
            peptide_or_protein_molecules=dict(numerator=n_peptide, denominator=n_a),
            median_T_external=float(np.median(t_external)),
            median_T_rest_of_A=float(np.median(t_other_a))),
        years=dict(
            A_approval=dict(min=years_a[0], median=int(np.median(years_a)),
                            max=years_a[-1], all=years_a),
            B_i_trial_start=dict(min=years_bi[0], median=int(np.median(years_bi)),
                                 max=years_bi[-1], all=years_bi),
            B_ii_paper=dict(min=years_bii[0], median=int(np.median(years_bii)),
                            max=years_bii[-1], all=years_bii)),
        rule=rule,
        rule_text=("T predicts if the median T of Set A is at least 2 points above the "
                   "median of Set B, the AUC is at least 0.75, and each set has at least "
                   "15 entries (INCLUSION.md, written before collection)."),
    )

    with (HERE / "results.json").open("w") as f:
        json.dump(results, f, indent=2)
        f.write("\n")

    draw(results, A, B)

    print(json.dumps({k: results[k] for k in
                      ("n", "median_T", "mann_whitney", "auc", "rule")}, indent=2))


def draw(results, A, B):
    FIGDIR.mkdir(exist_ok=True)
    pf.use_style()
    fig, ax = pf.figure(pf.ONE_COLUMN_MM, 66)

    ts = np.arange(7)
    count_a = np.array([int((A == t).sum()) for t in ts])
    count_b = np.array([int((B == t).sum()) for t in ts])
    assert count_a.sum() == len(A) and count_b.sum() == len(B)

    blue, vermillion = pf.OKABE_ITO["blue"], pf.OKABE_ITO["vermillion"]
    dodge = 0.12
    # Stems first, so the markers sit on top where the two series meet.
    for t, c in zip(ts, count_b):
        if c:
            ax.plot([t + dodge, t + dodge], [0, c], color=vermillion, lw=0.8, zorder=1)
    for t, c in zip(ts, count_a):
        if c:
            ax.plot([t - dodge, t - dodge], [0, c], color=blue, lw=0.8, zorder=1)
    ax.plot(ts - dodge, count_a, "o", color=blue, markersize=4, zorder=3,
            label=f"became an approved drug (n = {len(A)})")
    ax.plot(ts + dodge, count_b, "s", color="white", markeredgecolor=vermillion,
            markeredgewidth=0.8, markersize=4, zorder=3,
            label=f"did not become a drug (n = {len(B)})")

    med_a = results["median_T"]["A"]
    med_b = results["median_T"]["B"]
    ymax = max(count_a.max(), count_b.max())
    ax.set_ylim(0, ymax + 3.4)
    ax.set_xlim(-0.6, 6.6)
    ax.set_xticks(ts)
    ax.set_xlabel("Transferability score T (points out of 6; higher means easier to move)")
    ax.set_ylabel("Entries (n)")
    ax.set_yticks(range(0, ymax + 2, 2))
    for t in ts:
        ax.axvline(t, color="#DDDDDD", lw=0.3, zorder=0)

    # The two medians, each as a dashed line at its own T with its label on that line.
    ax.plot([med_a, med_a], [0, ymax + 1.9], color=blue, lw=0.6, ls=(0, (3, 2)),
            zorder=2)
    ax.text(med_a, ymax + 2.1, f"median T = {med_a:g}", color=blue, ha="center",
            va="bottom")
    ax.plot([med_b, med_b], [0, ymax + 0.8], color=vermillion, lw=0.6, ls=(0, (3, 2)),
            zorder=2)
    ax.text(med_b, ymax + 1.0, f"median T = {med_b:g}", color=vermillion, ha="center",
            va="bottom")

    pf.shared_legend(fig, ncol=2)
    pf.save(fig, str(FIGDIR / "transferability_t_distribution_backtest"))


if __name__ == "__main__":
    main()
