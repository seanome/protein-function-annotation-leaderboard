# Inclusion rules for the transferability back-test

Written 2026-10-19, before any entry was collected. These rules are fixed. If a rule turns
out to be wrong, the result is still reported under it, with the problem named.

The test: does the transferability score T, scored blind to the outcome, separate
animal-derived molecules that became approved drugs from animal-derived programmes that did
not? T is the six-point score from `twenty-diseases.md` (2026-10-19):

| | Criterion | 1 point if |
|---|---|---|
| T1 | Genetic basis | one gene, not polygenic or life-history |
| T2 | Location | secreted, not intracellular |
| T3 | Function | catalytic or binding, not structural or network |
| T4 | Target | conserved across phyla, not a species-specific partner |
| T5 | Size | peptide or single domain |
| T6 | Dependencies | works alone, no pathway or cofactors needed |

## Set A, approved

An entry is in Set A if all of these hold.

1. A drug with marketing approval from the FDA or the EMA, confirmed on that agency's own
   record (Drugs@FDA or the EMA medicine page) opened today, or in a PubMed-indexed paper
   opened today that names the agency and the year.
2. The active molecule comes from an animal: venom, saliva, other secretion, or tissue. The
   molecule may be made synthetically or recombinantly, and may be a modified version of
   the animal molecule, as long as the animal molecule is the starting point and the source
   species is named.
3. Excluded: drugs whose active molecule is a human protein, even when it is produced in an
   animal or in animal cells. Also excluded: molecules from plants, fungi, bacteria,
   archaea and viruses, which are not animals.
4. Excluded: whole-animal antisera and polyclonal antibody products, where the active
   molecule is the animal's immune response to an antigen rather than an evolved animal
   molecule with its own function.
5. One row per drug. The source species is named by binomial where the source paper or
   agency record names the species.

Recorded per entry: drug, source species, molecule, human target, approval year, agency,
citation (Chicago, with the DOI or PMID as the PubMed tool returned it).

## Set B, did not become a drug

Two kinds, recorded separately.

### B(i), reached a registered clinical trial, no approval

1. An animal-derived candidate molecule, meeting Set A rules 2, 3 and 4, with at least one
   study registered on ClinicalTrials.gov, read today through the v2 API, where the
   candidate is the intervention.
2. No marketing approval from the FDA or the EMA for that molecule in any indication, as of
   today.
3. Recorded: the trial identifier, the trial's start year, the phase and the status.
4. One row per candidate molecule, not per trial. Where several trials exist, the earliest
   start year that the search returned is recorded, and that trial's identifier is the one
   given.

### B(ii), a named gene from a tolerant animal, no candidate in trials

1. A molecule or gene named in a PubMed-indexed paper opened today, from a named animal
   species, proposed there as a route to treating a human disease.
2. No candidate in a registered clinical trial, checked today on ClinicalTrials.gov.
3. These are the rows of `twenty-diseases.md` that reach a named molecule but no approval,
   plus the naked mole rat hyaluronan route and tardigrade Dsup.
4. Recorded: the year of the paper that names the molecule as a route, as the start year, so
   Set A approval years and Set B start years can be compared.

Set B is the weak set: "did not become a drug" has no registry, so Set B is whatever could
be found under these rules. The number of entries excluded, and why, is in `excluded.csv`.

## Blind scoring

The entries go to a scoring agent that has not seen the sets, as `blind_input.csv`: an
opaque id and the mechanism description only (species, molecule, where it acts, what it
binds, size, whether it needs partners). No drug name, no outcome, no word that gives the
outcome away. The row order is shuffled with a recorded seed. The id-to-entry key is held
back in `key.csv` and is not given to the scoring agent. The outcome is rejoined only after
the scores come back.

## Pass/fail rule, fixed now

T predicts if all three hold:

* the median T of Set A is at least 2 points above the median T of Set B,
* the AUC is at least 0.75,
* each set has at least 15 entries.

If T predicts, it can rank programmes. If it does not, T stays a checklist of six questions
and is never used as a ranking, and the Idea Scoring Board says so.

If a set has fewer than 15 entries, the numbers are reported and the test is called
underpowered, not called either way.

## Known weakness, stated before the run

Set A is survivors and Set B is whatever could be found. The two sets differ in age and in
how well they are documented. The approval years in Set A and the start years in Set B are
both reported so the reader can see the gap.
