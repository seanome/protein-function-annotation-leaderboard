# Independent check of the tolerance atlas

Checked 2026-10-19 by a second agent that did not read the atlas README or the PR body. It
worked from `atlas.csv`, `queries.csv` and `VERIFY_TODO.md`, and from the papers
themselves. Script: `scripts/check_tolerance_atlas.py`. Row-by-row result:
`check_rows.csv`.

What was opened today: every PMID in columns 4, 5, 6, 7 and 8 (96 of them, all 96
returned by the NCBI E-utilities), every UniProt accession in column 6 (53 of them, all 53
returned by `rest.uniprot.org`), the PubMed Central full text of the 12 papers whose
number the atlas says is in the body, and the two FDA Drugs@FDA pages in column 8.

## Totals

| | count |
| --- | --- |
| checks HELD | 182 |
| checks FAILED | 5 |
| checks COULD NOT OPEN | 0 |
| nothing cited for that column, so nothing to open | 185 |
| total, 62 rows x 6 columns | 372 |

Every identifier resolved, and every one resolved to the paper the row cites. No PMID
pointed at a different paper and no DOI disagreed with its PMID. Two apparent DOI
mismatches turned out to be the checker's own parsing: PMID 10942883 in row 59 carries the
DOI `10.1002/1097-4695(20000905)44:3<289::aid-neu1>3.0.co;2-h`, whose semicolon the first
version of the script read as the end of the DOI, and PMID 28629163 in row 13 and PMID
34681069 in row 49 are each cited by PMID alone, with no DOI beside them.

## The five failures

Four of the five are the same mistake in column 7, the human target: the paper opened
shows a human protein, or a synthetic drug, acting in human cells, not the tolerant
animal's own molecule. The method spec asks for "a paper shows the molecule acting in a
human cell, on a human receptor, or in a human assay", and "the molecule" is the one named
in column 6.

| row | disease | column | what the paper actually says |
| --- | --- | --- | --- |
| 3 | D03 cardiac glycoside poisoning | 7 human target | The pump put into human embryonic kidney cells is the *Drosophila* Na,K-ATPase carrying the monarch's N122H change, not the *Danaus plexippus* protein Q2XN70 of column 6 (Holzinger and Wink, PMID 24227116). |
| 23 | D23 iron overload | 7 human target | The paper measures bat hepcidin mRNA after an iron challenge and compares bats with healthy humans (Stasiak, PMID 29594459). No bat molecule acts in a human cell or assay. The accession read is human hepcidin P81172, and the row itself says no *Desmodus rotundus* entry exists. |
| 27 | D27 hypercapnia | 7 human target | The paper tests aprepitant, a synthetic NK1 receptor antagonist, on human pre-B leukaemia cells (Bayati, PMID 27609608). Naked mole-rat substance P was not put into a human cell or assay. |
| 33 | D32 myocardial infarction | 7 human target | The paper raises human SIRT3 in human AC16 cardiomyocytes through HOTAIR and FUS (Liu, PMID 36991356). The *Urocitellus parryii* protein was not tested in a human cell. |
| 30 | D30 hypercholesterolaemia | 3 the number | The 1.5- to 2.0-fold syndecan-4 difference, the faster proliferation of the susceptible cells and their lower cell-surface heparan sulfate all appear in the abstract (Bortoff and Wagner, PMID 15596055). The added sentence that the two breeds show no consistent difference in plasma cholesterol appears in neither cited paper. |

Row 29, hyperoxaluria, is the fifth case of the same pattern but was left as held, because
the row says plainly what it is: the human evidence is a randomised trial of 80 patients
given the whole bacterium *Oxalobacter formigenes* (82.5% of patients with hyperoxaluria
down to 15% at one month, PMID 25645091), not the purified enzymes of column 6, and not
anything belonging to the woodrat *Neotoma albigula* in column 2. A reader who wants step 4
to mean the tolerant animal's own molecule should drop this row too.

## Smaller mismatches that do not change a count

* Row 3: the paper says untransfected controls showed "substantial mortality"; the row
  says they died.
* Row 28: the 56% and 11% alive at 72 hours are survival relative to the polyQ19 control,
  not absolute shares of cells. The median survival of 70 against 20 hours and the hazards
  of 1.9 against 7.2 are as the row gives them.
* Row 33: 6.7 divided by 0.9 is 7.4-fold, and the row says about 8-fold.
* Row 35: -3.8% against -12.1% is a 69% cut, and the row says 68%.
* Row 51: the claim that the bear and human IGF-1 peptide sequences are identical is not
  in the abstract of PMID 35077510 and was not checked.
* Row 56: the full text of PMID 23018966 does name *Acomys kempi* and *Acomys percivali*,
  so the two binomials the species column says were unavailable could have been read.
* Row 60: PMID 23671090 does have a PubMed Central full text, PMC3670376, although the row
  says the tool returned none.
* Row 18: the comparison with human plasma urea in kidney failure is the atlas's own, not
  the dogfish paper's.

## Recomputed counts

Counted from `atlas.csv` with this checker's own code, before `counts.json` was opened. A
disease clears a step only if one row clears that step and every earlier one.

| step | recomputed | counts.json | match |
| --- | --- | --- | --- |
| diseases in scope | 60 | 60 | yes |
| with a named species and a measured number | 42 | 42 | yes |
| with a named molecule | 29 | 29 | yes |
| moved to a human target | 7 | 7 | yes |
| with an approved drug | 0 | 0 | yes |
| negative rows (diseases with no number) | 18 | 18 | yes |

Restricted to rows scoring 5 or 6 on transferability:

| step | recomputed | counts.json | match |
| --- | --- | --- | --- |
| with a named species and a measured number | 7 | 7 | yes |
| with a named molecule | 7 | 7 | yes |
| moved to a human target | 3 | 3 | yes |
| with an approved drug | 0 | 0 | yes |

Every count matches. Two further numbers in `counts.json` also match: 62 rows in all, 44
of them with a number, and 235 queries in `queries.csv`. Every one of the 60 diseases has
at least three queries recorded, and the fewest for any one disease is three.

Two labels are loose rather than wrong. "Negative rows" is 18 when counted as diseases
with no measured number, but only 17 rows carry the species "none found": the eighteenth
disease is D02, scorpion sting, which names *Onychomys torridus* and has no number. The
key `diseases_with_no_species_and_no_number` therefore counts diseases without a number,
some of which do name a species. Separately, the atlas has 62 rows for 60 diseases because
D31 (ischaemic stroke) and D42 (freezing injury in cryopreservation) each carry two
species.

## How much the counts depend on the counting rule

This checker read column 6, 7 and 8 as a yes only when the cell begins with "yes", and
read a row as having a number unless the species is "none found" or the measured column
begins with "no number", "not found" or "none found". Three places where a stricter
reading of the method spec gives a different count:

* Column 6 as written asks for "the gene or protein name and its UniProt accession, read
  from UniProt today". Two diseases clear step 3 on a cell that begins with "yes" and
  names no accession: D05, the softshell clam sodium channel, and D42, the bacterial
  ice-binding protein FfIBP, whose cell says "yes by name, no by accession". Requiring an
  accession takes step 3 from 29 diseases to 27. Both cells say plainly that no accession
  was read.
* Three rows clear step 2 on a weak number. D37 gives the habitat altitude of 3600 m and
  no outcome value, and its own fit score says so. D49 gives a zero change in muscle fibre
  cross-sectional area in seven bears with no percentage. D54 names a genus, *Acomys* spp.,
  where the spec asks for a binomial, although the paper's full text does carry the two
  binomials. Dropping all three takes step 2 from 42 to 39.
* The build script carries two overrides, in `scripts/transferability.py`, that rewrite a
  "yes" into a "no" for D22's human target and D27's drug. Both are already written into
  `atlas.csv`, so this checker read the overridden text, and the matching counts are not
  an accident of reading the raw collection instead. The same reasoning applied to D27's
  drug column is the reasoning missing from the four human-target failures above.

## Which funnel counts would change

Dropping the four rows whose human-target check failed takes step 4 from 7 diseases to 3:
D29 hyperoxaluria, D43 severe dehydration and D46 acute radiation syndrome. Of those,
only two rest on the tolerant animal's own molecule in human cells: the brine shrimp
protein p26 in human 293H cells (PMID 15963489) and the tardigrade damage suppressor
protein in human HEK293 cells (PMID 27649274). If row 29 goes as well, step 4 is 2.

In the transferability 5-6 restriction, step 4 falls from 3 to 1, leaving only D46, the
tardigrade protein. The two that go are D23 and D27, both failures above.

No other count moves. Steps 1, 2 and 3 and the approved-drug count of 0 are untouched,
because the five failures are all in columns 3 and 7 and none of them removes a species, a
number or a molecule.

One thing the drug count of 0 hides: the atlas does hold one approved drug, rasburicase
(ELITEK, US FDA BLA 103946, Sanofi Synthelabo, original approval 12 July 2002, page
fetched today). It sits in row 21, hyperuricaemia and gout, whose species column is "none
found", so the nested funnel never reaches step 5 for it. The drug is a fungal and
mammalian enzyme resurrected in the lab, not a molecule taken from a tolerant animal, so
the zero is the right answer to the question the funnel asks. The approval itself held on
checking, as did the aprepitant page in row 27 (EMEND, NDA 021549, original approval 26
March 2003), which that row correctly does not claim as a step-5 pass.

## What this checker did not do

* Did not read the atlas README, the PR body, or `scripts/build_tolerance_atlas.py`'s
  counting code before recomputing the counts.
* Did not change `atlas.csv`, `counts.json`, `scope.csv` or `queries.csv`.
* Did not re-run the PubMed queries in `queries.csv`, so the claim that a given search
  returned no articles rests on the atlas's own record. Rows 10, 20, 22, 24 and 57 say so.
* Did not check the transferability T or fit N scores in columns 9 and 10, or the lab
  affiliations behind N3.
* Rows 28 and 29 needed their numbers from the paper body, and the E-utilities reply for
  both held only the abstract. Both were read from the PubMed Central web page instead.
