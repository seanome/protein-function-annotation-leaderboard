# The tolerance atlas: how many human diseases clear each of five steps

Written 2026-10-19 for Olga Botvinnik (Seanome). The belief under test: for any insult a
lineage met repeatedly over deep time, some lineage evolved a tolerance. Each disease is
checked against five steps: (1) it is in scope, (2) a named species tolerates the insult
with a measured number, (3) the molecule behind the tolerance is named, (4) that molecule
has been shown to act in a human cell, on a human receptor or in a human assay, (5) it
became an approved drug.

This first section is committed on its own, before any species search, so that git
history shows the scope list and the counting rules were fixed first.

## Scope list, fixed before any species search

`scope.csv` lists 60 diseases. A disease is in scope if it is a physical or chemical
insult that many lineages have met: a toxin or metabolite that accumulates, loss of
oxygen, freezing, drying, radiation, long immobility, injury that scars instead of
regrowing, or cancer in a large or long-lived body. No disease is added after a species
for it turns up.

| insult class | diseases |
|---|---|
| toxin | 16 (D01-D16) |
| metabolite that accumulates | 14 (D17-D30) |
| loss of oxygen | 10 (D31-D40) |
| freezing | 3 (D41-D43) |
| drying | 2 (D44-D45) |
| radiation | 2 (D46-D47) |
| long immobility | 3 (D48-D50) |
| injury that scars instead of regrowing | 9 (D51-D59) |
| cancer in a large or long-lived body | 1 (D60) |

D28 (protein aggregation) is placed under "metabolite that accumulates" because the
damaging thing is a misfolded protein that builds up. That placement is a judgement and
may be argued.

### Rows of `twenty-diseases.md` that are out of scope

These eight starting rows fail the scope rule. Their claims are still checked (they feed
the transferability back-test), in `starting-twenty.csv`, but they are not counted in
the funnel.

| row | disease | why out of scope |
|---|---|---|
| 2 | severe chronic pain | a symptom, not an insult from the list |
| 3 | hypertension | not an insult from the list |
| 4 | thrombosis (in general) | not an insult from the list; clotting during immobility is D50 |
| 5 | antibiotic-resistant infection | a pathogen that co-evolves against us |
| 7 | vaccine and biologic stability | not a human disease; drying is covered by D44 and D45 |
| 10 | autoimmune and inflammatory disease | autoimmunity |
| 17 | transplant rejection | an immune response to another person's tissue, not an insult from the list |
| 20 | ageing and healthspan | not a disease and not an insult from the list; cancer in a long-lived body is D60 |

The other twelve starting rows map to D01, D19, D22, D28, D31 and D32, D36, D42, D46,
D48, D49, D51 to D54, D60.

## Counting rules, fixed before any species search

- A row is one disease and one species. A disease with no species found gets one
  negative row: species "none found", with its PubMed queries in `queries.csv`.
- Step 2 needs a tolerance: the species withstands the insult, with a number, a unit and a
  condition from one primary paper. A species that is only the source of a molecule used
  against the disease (for example a venom peptide) does not clear step 2 by being the
  source. Such a row is kept, marked "source species, no tolerance measured", and counted
  separately.
- The funnel is nested. A disease counts at step k if at least one of its rows clears
  step 2 and every step up to k. The same counts without nesting are reported beside it.
- Transferability T (0 to 6) is scored from columns 6 and 7 only, before column 8 is
  read for that row. One point each: T1 one gene; T2 secreted; T3 catalytic or binding;
  T4 target conserved across phyla; T5 peptide or single domain; T6 works alone, no
  partner or cofactor needed. A row with no named molecule gets T = 0 with the points
  written as not scorable.
- Nature-solved-it fit N (0 to 3): N1 a named species with a measured number; N2 the
  output is a molecule; N3 a lab exists that holds the animal and runs the assay (named
  from the cited paper's affiliation).

## Rule, fixed before any species search

- The atlas supports the belief as stated if at least half of the in-scope diseases (30
  of 60) have a named species with a measured number and a citation.
- It supports the narrow version ("venoms and secretions are the only class with a drug
  record") if every row with an approved drug has T of 5 or 6, and no row with T of 2 or
  below does.
- If fewer than a quarter (15 of 60) clear step 2, the belief is a weak foundation.
