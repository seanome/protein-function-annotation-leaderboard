# Transferability back-test: the two sets, collected blind of the scoring

Does the transferability score T, scored without knowing the outcome, separate
animal-derived molecules that became approved drugs from animal-derived programmes that
did not? `INCLUSION.md` holds the inclusion rules and the pass/fail rule, written and
committed before any entry was collected. This page says what is in the files.

This half of the work is collection only. T has not been scored here, and no statistic
has been computed. A separate agent that has not seen these sets scores T from
`blind_input.csv`.

## Files

| File | What it holds |
|---|---|
| `INCLUSION.md` | The inclusion rules and the pass/fail rule, fixed before collection |
| `build_sets.py` | The records read on 2026-10-19, and the code that writes the CSVs |
| `sets.csv` | One row per entry, with the set, species, molecule, human target, year, agency or registry record, and citation |
| `excluded.csv` | Candidates checked and left out, with the reason and the record checked |
| `blind_input.csv` | Opaque id and mechanism description only, row order shuffled with seed 20261019 |
| `key.csv` | Blind id to entry, held back from the scoring agent |

## Counts

| Set | Entries |
|---|---|
| A, approved | 19 |
| B(i), reached a registered trial, no approval | 15 |
| B(ii), named gene from a tolerant animal, no candidate in trials | 9 |
| Excluded | 30 |

Set B has 24 entries in all. The pass/fail rule asks for at least 15 per set; both sets
clear that, and the two kinds of Set B entry are kept separate so the analysis can report
them apart as well as together.

## Where the facts come from

Every year, application number, trial identifier, DOI and PMID in `sets.csv` was read on
2026-10-19 from one of three places, named per row in the `citation` and `source_type`
columns:

* Drugs@FDA and FDA product labels, through the openFDA `drugsfda` and `label` APIs.
* ClinicalTrials.gov, through its v2 API.
* PubMed, through the PubMed tool, which returned the DOI and PMID as written.

Where a source does not give a species binomial, the species column says so in the words
the source used (for example "pig (the label says porcine intestinal mucosa; no binomial
given)"). Rows whose origin claim rests on a review rather than the primary paper are
marked "review" in `source_type`.

## How the blind file was made

`blind_input.csv` carries the mechanism description and nothing else: species, molecule,
where it acts, what it binds, size, and whether it needs partners. It holds no drug name,
no brand name, no agency, no year, no trial identifier, and no word that says how the
entry ended up. Row order is shuffled with seed 20261019, recorded in `key.csv` and in
`build_sets.py`. The descriptions are the same text as the `mechanism_description` column
of `sets.csv`, so the scoring agent and the rejoined table cannot drift apart.

## Known weakness, stated before the run

Set A is survivors; Set B is whatever could be found. The sets differ in age and in how
well they are documented. `sets.csv` carries the approval year for Set A and the start
year for Set B so the gap can be read off directly.

## Still to do

Blind scoring of T, the rejoin, the figure, the three statistics (median difference,
Mann-Whitney U with its p value, AUC), the share of Set A that is a venom, saliva or
secretion peptide, and the second-agent check.
