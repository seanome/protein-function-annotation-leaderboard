"""Build the tolerance atlas table, the query log, the counts and the funnel figure.

Input: results/nature-solved-it/atlas/rows/block-*.jsonl, one JSON object per disease,
plus results/nature-solved-it/atlas/scope.csv, the 60 in-scope diseases fixed before any
species search.

Output, all under results/nature-solved-it/atlas/:
  atlas.csv        eleven columns, one row per disease and species
  queries.csv      every PubMed query, verbatim, with the disease it was run for
  counts.json      the five funnel counts, the negative-row count, and the same
                   restricted to rows whose transferability score is 5 or 6
  funnel_tolerance_atlas_pubmed_2026-10-19.{pdf,svg,png}

Run: python scripts/build_tolerance_atlas.py
"""

import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import pubfig as pf  # noqa: E402
from transferability import DRUG_OVERRIDE, HUMAN_OVERRIDE, POINTS  # noqa: E402

ATLAS = Path(__file__).resolve().parents[1] / "results" / "nature-solved-it" / "atlas"
ROWS = ATLAS / "rows"

COLUMNS = [
    "human_disease_and_insult_class",
    "species",
    "measured_ability",
    "citation_for_the_number",
    "independent_origins",
    "molecule_named",
    "moved_to_a_human_target",
    "approved_drug",
    "transferability_T",
    "nature_solved_it_fit_N",
    "verified",
]

STEP_LABELS = [
    "in scope",
    "named species\nwith a number",
    "molecule\nnamed",
    "reported human-\ntarget claim",
    "approved\ndrug",
]


def points(text, letter):
    """Count the scored points in a string like 'T1=1 one gene; T2=0 intracellular'."""
    return len(re.findall(letter + r"[1-6]\s*=\s*1", str(text)))


def says_yes(text):
    return str(text).strip().lower().startswith("yes")


def has_number(species, measured):
    s = str(species).strip().lower()
    m = str(measured).strip().lower()
    if s.startswith("none found") or s.startswith("not found") or s == "":
        return False
    if m.startswith("not found") or m.startswith("no number") or m.startswith("none found"):
        return False
    return True


def load():
    scope = pd.read_csv(ATLAS / "scope.csv")
    label = {
        r.disease_id: f"{r.human_disease} [{r.insult_class}]" for r in scope.itertuples()
    }
    records, queries = [], []
    for path in sorted(ROWS.glob("block-*.jsonl")):
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            d = json.loads(line)
            did = d["id"]
            if did not in label:
                raise SystemExit(f"{path.name}: {did} is not in scope.csv")
            for q in d["queries"]:
                queries.append({"disease_id": did, "human_disease": label[did], "pubmed_query": q})
            for i, row in enumerate(d["rows"]):
                key = (did, i)
                if key not in POINTS:
                    raise SystemExit(f"{did} row {i} has no transferability score in "
                                     "scripts/transferability.py")
                t_text = POINTS[key]
                t_points = points(t_text, "T")
                n_points = points(row["N"], "N")
                drug = DRUG_OVERRIDE.get(key, row["drug"])
                human = HUMAN_OVERRIDE.get(key, row["human"])
                records.append(
                    {
                        "disease_id": did,
                        "human_disease_and_insult_class": label[did],
                        "species": row["species"],
                        "measured_ability": row["measured"],
                        "citation_for_the_number": row["cite"],
                        "independent_origins": row["origins"],
                        "molecule_named": row["molecule"],
                        "moved_to_a_human_target": human,
                        "approved_drug": drug,
                        "transferability_T": f"T={t_points}. {t_text}",
                        "nature_solved_it_fit_N": f"N={n_points}. {row['N']}",
                        "verified": row["verified"],
                        "T_points": t_points,
                        "N_points": n_points,
                        "clears_step2": has_number(row["species"], row["measured"]),
                        "molecule_yes": says_yes(row["molecule"]),
                        "human_yes": says_yes(human),
                        "drug_yes": says_yes(drug),
                    }
                )
    atlas = pd.DataFrame(records)
    missing = sorted(set(label) - set(atlas.disease_id))
    if missing:
        raise SystemExit("no rows written for: " + ", ".join(missing))
    return scope, atlas, pd.DataFrame(queries)


def funnel(atlas, n_scope, subset=None):
    """Five nested counts of diseases, plus the number of diseases with no number."""
    df = atlas if subset is None else atlas[subset]
    step2 = set(df.loc[df.clears_step2, "disease_id"])
    step3 = set(df.loc[df.clears_step2 & df.molecule_yes, "disease_id"])
    step4 = set(df.loc[df.clears_step2 & df.molecule_yes & df.human_yes, "disease_id"])
    step5 = set(
        df.loc[df.clears_step2 & df.molecule_yes & df.human_yes & df.drug_yes, "disease_id"]
    )
    return [n_scope, len(step2), len(step3), len(step4), len(step5)]


def draw(all_counts, t56_counts, negatives, negatives_t56, out_stem):
    pf.use_style()
    fig, axs = pf.figure(pf.ONE_COLUMN_MM, 110, nrows=2, sharex=True)
    blue = pf.OKABE_ITO["blue"]
    vermillion = pf.OKABE_ITO["vermillion"]
    x = range(5)
    for ax, counts, neg, title in (
        (axs[0], all_counts, negatives, "every row"),
        (axs[1], t56_counts, negatives_t56, "rows with transferability 5 or 6"),
    ):
        ax.bar(x, counts, width=0.68, color=blue, label="diseases reported to clear the step (filled bar)")
        ax.bar(
            [5],
            [neg],
            width=0.42,
            facecolor="none",
            edgecolor=vermillion,
            linewidth=1,
            linestyle="--",
            label="diseases with no species and no number in this set (hollow dashed bar)",
        )
        for xi, c in zip(list(x) + [5], list(counts) + [neg]):
            ax.text(xi, c + 1.2, str(c), ha="center", va="bottom")
        ax.set_ylim(0, max(all_counts + [negatives]) * 1.22)
        ax.set_ylabel("Diseases (n)")
        ax.text(0.98, 0.95, title, transform=ax.transAxes, va="top", ha="right")
    axs[1].set_xticks(list(range(6)))
    axs[1].set_xticklabels(STEP_LABELS + ["no species\nwith a number"])
    axs[1].set_xlabel("Step of the five-step test, left to right")
    pf.shared_legend(fig, ncol=1)
    for ax, letter in zip(axs, "ab"):
        pf.panel_label(ax, letter)
    pf.save(fig, str(out_stem))


IDENT = re.compile(r"(PMID\s*\d+|https?://doi\.org/\S+|UniProt\s+[A-Z0-9]{6,10})")


def write_verify_todo(atlas, counts):
    """List, per row, every identifier a second agent has to open."""
    lines = [
        "# What the checker must open",
        "",
        "Written 2026-10-19 by the agent that built the atlas. A second agent that has not",
        "seen README.md does this list. Expect some rows to fail; report them, do not fix",
        "them silently.",
        "",
        "## What to check in every row of atlas.csv",
        "",
        "1. Open the citation in column 4 and confirm the number in column 3 is in that",
        "   paper, under that condition. A number that is in the body and not the abstract",
        "   counts only if you opened the body.",
        "2. Open the phylogeny paper in column 5 and confirm the count of separate origins.",
        "3. Read the UniProt accession in column 6 today and confirm the protein and the",
        "   species match the row.",
        "4. Open the citation in column 7 and confirm it shows the molecule acting in a",
        "   human cell, on a human receptor, or in a human assay.",
        "5. Open the agency page in column 8 and confirm the approval and the agency.",
        "6. Recompute the five funnel counts and the negative-row count from atlas.csv with",
        "   scripts/build_tolerance_atlas.py, and compare them with counts.json and the",
        "   numbers in README.md.",
        "7. Check that every row marked 'none found' really has at least three queries in",
        "   queries.csv. The atlas names the rows where it has fewer.",
        "",
        "## Identifiers, row by row",
        "",
    ]
    for r in atlas.itertuples():
        ids = []
        for col in ("citation_for_the_number", "independent_origins", "molecule_named",
                    "moved_to_a_human_target", "approved_drug"):
            ids += IDENT.findall(str(getattr(r, col)))
        seen, uniq = set(), []
        for i in ids:
            if i not in seen:
                seen.add(i)
                uniq.append(i)
        lines.append(f"- **{r.disease_id} {r.human_disease_and_insult_class}** "
                     f"({r.species}): " + ("; ".join(uniq) if uniq else "no identifier to open"))
    lines += [
        "",
        "## Counts to recompute",
        "",
        "```json",
        json.dumps(counts["funnel_all_rows"], indent=2),
        "```",
        "",
        f"Negative rows (diseases with no species and no number): "
        f"{counts['diseases_with_no_species_and_no_number']}.",
        "",
    ]
    (ATLAS / "VERIFY_TODO.md").write_text("\n".join(lines) + "\n")


def main():
    scope, atlas, queries = load()
    n_scope = len(scope)
    all_counts = funnel(atlas, n_scope)
    t56 = atlas.T_points >= 5
    t56_counts = funnel(atlas, n_scope, subset=t56)
    with_number = set(atlas.loc[atlas.clears_step2, "disease_id"])
    negatives = n_scope - len(with_number)
    with_number_t56 = set(atlas.loc[t56 & atlas.clears_step2, "disease_id"])
    negatives_t56 = n_scope - len(with_number_t56)

    atlas[COLUMNS].to_csv(ATLAS / "atlas.csv", index=False)
    queries.to_csv(ATLAS / "queries.csv", index=False)
    counts = {
        "diseases_in_scope": n_scope,
        "funnel_all_rows": dict(zip(
            ["in_scope", "named_species_with_a_number", "molecule_named",
             "shown_on_a_human_target", "approved_drug"], all_counts)),
        "funnel_rows_with_T_5_or_6": dict(zip(
            ["in_scope", "named_species_with_a_number", "molecule_named",
             "shown_on_a_human_target", "approved_drug"], t56_counts)),
        "diseases_with_no_species_and_no_number": negatives,
        "diseases_with_no_T5_or_T6_row_with_a_number": negatives_t56,
        "rows_total": len(atlas),
        "rows_with_a_number": int(atlas.clears_step2.sum()),
        "queries_total": len(queries),
        "drug_rows": [
            {"disease": r.human_disease_and_insult_class, "species": r.species,
             "T": int(r.T_points), "drug": r.approved_drug}
            for r in atlas[atlas.drug_yes].itertuples()
        ],
        "rows_with_T_2_or_below_and_a_drug": int(
            ((atlas.T_points <= 2) & atlas.drug_yes).sum()
        ),
    }
    (ATLAS / "counts.json").write_text(json.dumps(counts, indent=2) + "\n")

    write_verify_todo(atlas, counts)
    draw(all_counts, t56_counts, negatives, negatives_t56,
         ATLAS / "funnel_tolerance_atlas_pubmed_2026-10-19")
    print(json.dumps(counts, indent=2)[:4000])


if __name__ == "__main__":
    main()
