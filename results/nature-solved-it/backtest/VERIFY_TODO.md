# What the checking agent does

For a second agent that has not read `README.md`, `results.json` or `scored.csv`. Report
what held and what failed. Expect some rows to fail; report them, do not fix them quietly.

## 1. Reopen every approval and trial record

For each of the 19 Set A rows in `sets.csv`, open the record named in `record_id` on
Drugs@FDA (openFDA `drugsfda`, or the Drugs@FDA website) and confirm three things:

* the application number exists and belongs to that drug,
* the original approval date matches the year in the `year` column,
* the label or paper named in `citation` says the source species or tissue the `species`
  column claims, in the words the column quotes.

For each of the 15 Set B(i) rows, open the trial identifier in `record_id` on
ClinicalTrials.gov (v2 API) and confirm the start date, phase, status, and that the
intervention is the molecule in the `molecule` column. For the 9 Set B(ii) rows, confirm
the DOI or PMID resolves to the paper cited, that the paper names that molecule in that
species, and that a ClinicalTrials.gov search for the molecule still returns nothing.

Nine rows say the source gives no species binomial. Check that each of those really does
lack one in the source, rather than the binomial having been missed.

## 2. Re-score ten entries blind

Take ten `blind_id` rows from `blind_input.csv` at random, with the seed recorded, and
score T1 to T6 from the mechanism description alone, without opening `blind_scores.csv`,
`key.csv` or `scored.csv`. Then compare with `blind_scores.csv` and report how many of the
ten total T values match exactly, and which single points differ and why. The six points
are in `INCLUSION.md`.

## 3. Recompute the statistics from the CSVs

From `scored.csv` (or from `key.csv` joined to `blind_scores.csv` and `sets.csv`),
recompute and compare with `results.json`:

* n for Set A, Set B, Set B(i), Set B(ii),
* the median T of Set A and of Set B,
* Mann-Whitney U for Set A against Set B, both p values,
* the AUC, computed as U divided by the number of approved and not-approved pairs, and
  again by counting pairs directly,
* the venom, saliva and secretion numerator and denominator for Set A, under the
  definition written in `analyze.py`,
* the pass or fail of each of the three conditions in the rule in `INCLUSION.md`.

## 4. Check the claims in README.md against those numbers

Read each sentence of `README.md` that carries a number and check it against what you
recomputed, including the sentence that says what the result means for the Idea Scoring
Board. Report any sentence that claims more than its number supports.
