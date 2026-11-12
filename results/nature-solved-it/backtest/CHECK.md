# Independent check of the transferability back-test

Written 2026-10-19 by a second agent, following `VERIFY_TODO.md`. The order was fixed in
advance: the ten blind re-scores were written and saved before this agent opened
`sets.csv`, `key.csv`, `blind_scores.csv`, `scored.csv`, `results.json` or `README.md`.
Nothing in the back-test folder was edited. The files added by this check are `CHECK.md`,
`check_rows.csv`, `rescore.csv` (the blind re-score), `recheck.py` (the statistics),
`compare_scores.py` (the score comparison), `check_records.py` (the record queries),
`build_check_rows.py` (the table) and the raw API answers those queries returned,
`check_fda_raw.json` and `check_ctg_raw.json`.

What this check opened today: the openFDA `drugsfda` and `label` endpoints for the 19
Set A applications, the ClinicalTrials.gov v2 API for the 15 Set B(i) trials and for
searches of the 9 Set B(ii) molecules, openFDA `drugsfda` searches for each Set B(i)
molecule and its code name, and PubMed (through the PubMed MCP tools) for all 25 distinct
DOI/PMID pairs cited in `sets.csv`.

## 1. Every record, held or failed

Per-record detail, including what each record showed, is in `check_rows.csv`.

| Verdict | Count | Which |
|---|---|---|
| HELD | 43 | all 19 Set A, all 15 Set B(i), all 9 Set B(ii) |
| FAILED | 0 | — |
| COULD NOT OPEN | 0 records, 1 check | the EMA half of the "no approval" test for Set B(i) (below) |

Set A, 19 of 19 held. Each application number in `record_id` exists in Drugs@FDA, belongs
to the drug named in the row, and carries an original approval date whose year is the year
in the `year` column, to the day in every case: NDA018343 CAPOTEN 1981-04-06, NDA021773
BYETTA 2005-04-28, NDA021060 PRIALT 2004-12-28, NDA020718 INTEGRILIN 1998-05-18,
NDA020873 ANGIOMAX 2000-12-15, BLA021271 IPRIVASK 2003-04-04, NDA208114 DEFITELIO
2016-03-30, NDA017029 HEPARIN SODIUM 1972-02-22, ANDA089454 PROTAMINE SULFATE 1987-04-07,
NDA017808 MIACALCIN 1986-07-03, BLA021665 AMPHADASE 2004-10-26, BLA021640 VITRASE
2004-05-05, BLA020725 CREON 2009-04-30, BLA020032 SURVANTA 1991-07-01, BLA020521 INFASURF
1998-07-01, BLA020744 CUROSURF 1999-11-18, NDA201532 HALAVEN 2010-11-15, NDA207953
YONDELIS 2015-10-23, NDA021654 LOVAZA 2004-11-10.

Thirteen Set A rows cite an FDA label rather than a paper. The species or tissue wording
quoted in the `citation` column was found in the current openFDA label for all 13, in the
words the column quotes, including
"originally identified in the lizard Heloderma suspectum" (BYETTA), "piscivorous marine
snail, Conus magus" (PRIALT), "porcine intestinal mucosa" (heparin), "the sperm of salmon
and certain other species of fish" (protamine), "calcitonin of salmon origin",
"purified bovine testicular hyaluronidase" (AMPHADASE), "natural bovine lung extract"
(SURVANTA), "natural bovine lung" (INFASURF), "natural porcine lung" (CUROSURF),
"marine sponge Halichondria okadai" (HALAVEN) and "sourced from fish oils" (LOVAZA).
A brand-name query returned no label for six brands: CAPOTEN, INTEGRILIN, ANGIOMAX,
IPRIVASK, VITRASE and LOVAZA. LOVAZA's label came back under its generic name, with the
quoted wording. The other five rows cite a paper rather than a label, and each of those
papers opened and said what the row claims, so those rows hold on the paper.

Set B(i), 15 of 15 held. Every NCT id resolves, and the start date, phase, status and
intervention match the row: NCT02446340 dalazatide 2013-11 phase 1 completed, NCT00040573
131I-TM-601 2002-06 phase 1/2 completed, NCT00141011 Ancrod (Viprinex) 2005-09 phase 3
terminated, NCT00073554 alfimeprase 2003-06 phase 2 completed, NCT00638781 desmoteplase
2001-03 phase 2 completed, NCT00116012 rNAPc2 2005-06 phase 2 completed, NCT00563394
MSI-78 1994-08 phase 3 completed, NCT01481532 crotoxin 2026-10 phase 1 not yet recruiting,
NCT01578564 SOR-C13 2012-07 phase 1 completed, NCT00482937 CD-NP 2007-01 phase 1
completed, NCT04269954 batroxobin 2020-03-01 phase 4 unknown, NCT00027248 MBI 226 2000-09
phase 3 completed, NCT00022373 iseganan 2000-12 phase 3 unknown, NCT00317512 HBOC-201
2003-12 phase 2 completed, NCT00980083 exendin(9-39)amide 2004-10 phase 1 completed.
The crotoxin start date of 2026-10 is what the registry says today, which is what the row
records.

A Drugs@FDA search of substance, brand and generic name for each of the 15 molecules and
for its code name returned nothing, so no FDA approval was found for any of them. The same
query form returns records for captopril, exenatide and ziconotide, so the empty answers
are real and not a broken query. This is a name search: a molecule approved under a name
none of the sets records would be missed.

Set B(ii), 9 of 9 held. All nine papers opened in PubMed and name the molecule in that
species. Three rows claim the source gives no binomial and that is what the source shows:
the 2023 CAHS D paper says only "tardigrade" in its abstract, the 2012 opossum review
names Didelphidae and no binomial, and the 2004 antifreeze protein paper says "arctic
fish". The Dsup row says the binomial sits in the paper's keywords rather than the
abstract, and it does: Ramazzottius varieornatus. A ClinicalTrials.gov search for each of
the eight searchable molecules returned no study.

### What this check could not open, and what it flags

* The EMA side of the Set B(i) no-approval test. `INCLUSION.md` asks for no FDA or EMA
  approval. The EMA website returned HTTP 403 to every request from here, and the EMA has
  no open API matching openFDA, so only the FDA half was opened today. Two molecules make
  this matter: batroxobin is on the market through agencies outside the FDA and the EMA,
  and so is HBOC-201. Both rows hold under the rule as written, because the rule looks
  only at the FDA and the EMA, but "did not become a drug" is weaker for them than for the
  other 13.
* Row B115 is exendin(9-39), a shortened form of the same lizard peptide whose full-length
  form is approved as exenatide in row A02. Under the rule the two are different
  molecules, so both rows are placed correctly, but Set A and Set B(i) each hold a form of
  one animal peptide.
* Row B110 cenderitide is a chimera of a human peptide and a snake venom peptide. Set A
  rule 2 asks for the animal molecule as the starting point. Here only the tail comes from
  the animal, which is why the blind scorer gave it T1 = 0.
* The `year` column for Set A is the year of the application in `record_id`, not the year
  the molecule first reached the market. Heparin (1972) and protamine (an abbreviated
  application in 1987) were sold long before those dates. The column name and `year_type`
  say this, and nothing downstream treats these years as first market entry.
* Row A12 names Ovis aries in the `species` column while the cited paper says "sheep
  testicular hyaluronidase" and gives no binomial. The same cell says so in its
  parenthesis, so the claim is hedged, but the binomial is the agent's, not the paper's.

## 2. The ten blind re-scores against `blind_scores.csv`

Scored from the mechanism descriptions in `rescore_input.csv` only, saved before any other
file was opened, and compared with `compare_scores.py`.

**3 of the 10 T totals match exactly. Points agree on 52 of 60.**

| blind id | entry | programme | this check | agent | points that differ |
|---|---|---|---|---|---|
| E032 | A01 | captopril | 6 | 6 | none |
| E042 | B102 | chlorotoxin | 6 | 6 | none |
| E004 | B110 | cenderitide | 5 | 5 | none |
| E030 | B202 | tardigrade Dsup | 4 | 3 | T5: 1 here, 0 there |
| E021 | A13 | pancrelipase | 3 | 2 | T6: 1 here, 0 there |
| E014 | B206 | salamander nAG | 4 | 2 | T3 and T5: 1 here, 0 there |
| E031 | B105 | desmoteplase | 3 | 4 | T4: 0 here, 1 there |
| E006 | A12 | hyaluronidase, ovine | 6 | 5 | T5: 1 here, 0 there |
| E010 | B108 | crotoxin | 3 | 2 | T4: 1 here, 0 there |
| E003 | B106 | rNAPc2 | 5 | 6 | T6: 0 here, 1 there |

Agreement per point: T1 10/10, T2 10/10, T3 9/10, T4 8/10, T5 7/10, T6 8/10.

Why the eight single points differ:

* T5, size, four of the eight. The agent scores 0 whenever the description does not state
  a residue count or a domain count, so Dsup, nAG and ovine hyaluronidase all lose the
  point for saying "one chain" without a number. This check read "one chain" and "one
  catalytic chain" as a single domain and gave the point. Neither reading is wrong under
  the one-line rule in `INCLUSION.md`; the rule does not say whether an unstated size
  scores 0 or whether "one chain" is enough. T5 is the least reproducible of the six.
* T4, target conserved across phyla, two. For desmoteplase this check read "acts through a
  human partner protein, so it needs that protein to be there" as a species-specific
  partner and scored 0; the agent read plasminogen as a conserved target and scored 1. For
  crotoxin this check gave the point because lipid membranes are universal; the agent
  scored 0 because no defined target is named. The rule mixes two tests, how conserved the
  target is and whether it is a partner of one species, and the two tests point opposite
  ways on these rows.
* T6, works alone, two. For pancrelipase this check gave the point because each enzyme
  needs no cofactor; the agent scored 0 because many enzymes act together. For rNAPc2 this
  check scored 0 because the molecule needs an assembled two-protein complex; the agent
  gave the point because the description says "no cofactor".
* T3, catalysis or binding, one. For nAG the description says it is released from gland
  cells and is enough to let a stump regrow digits. This check counted a secreted ligand
  as binding; the agent scored 0 because the description names no binding event.

T1 and T2 reproduced perfectly. T4, T5 and T6 did not, and a one-point move on three of
six points swings a total by up to 2. Under the pass rule, which asks for a median gap of
2 points, T values that vary by this much between two scorers cannot support a ranking.

## 3. Recomputed statistics

`recheck.py` reads `scored.csv`, rebuilds T from `blind_scores.csv` joined through
`key.csv`, and recomputes everything with scipy 1.16.3. It does not import `analyze.py`.
Every number matches `results.json`.

| Quantity | Recomputed | `results.json` |
|---|---|---|
| n Set A / Set B / B(i) / B(ii) | 19 / 24 / 15 / 9 | same |
| median T, Set A | 5.0 | 5.0 |
| median T, Set B | 3.5 | 3.5 |
| median T, B(i) / B(ii) | 5.0 / 2.0 | same |
| mean T, A / B / B(i) / B(ii) | 3.895 / 3.375 / 4.267 / 1.889 | same |
| Mann-Whitney U, A over B | 259.5 | 259.5 |
| p, two-sided | 0.4396 | 0.4396 |
| p, one-sided A greater | 0.2198 | 0.2198 |
| pairs A > B / tied / total | 222 / 75 / 456 | same |
| AUC, U/(nA·nB) | 0.5691 | 0.5691 |
| AUC, counting pairs | 0.5691 | 0.5691 |
| A against B(i): U, AUC, p two-sided | 134.0, 0.4702, 0.7754 | same |
| A against B(ii): U, AUC, p two-sided | 125.5, 0.7339, 0.0476 | same |
| venom/saliva/secretion share of A | 6/19, with the internal one 7/19 | same |
| median T, external-secretion rows / rest of A | 6.0 / 3.0 | same |

The 43 rows of `scored.csv` also match `blind_scores.csv` and `key.csv` row for row: the
six points sum to the stored T in every row, and every `entry_id` and `set` agrees with
both `key.csv` and `sets.csv`.

## 4. The rule in `INCLUSION.md`, applied here

| Condition | Required | Found | Pass |
|---|---|---|---|
| median T of Set A above median of Set B | at least 2 points | 1.5 points (5.0 against 3.5) | no |
| AUC | at least 0.75 | 0.5691 | no |
| entries per set | at least 15 | 19 and 24 | yes |

Two of the three conditions fail, so T does not predict. The verdict in `results.json`
holds: `t_predicts` is false. Under the rule's own words, T stays a checklist of six
questions and is not used as a ranking.

**Would dropping failed rows change it?** No row failed, so nothing is dropped and the
verdict is unchanged. The nearest thing to a failure is the two molecules approved outside
the FDA and the EMA. Dropping them one at a time or together still fails the rule: drop
batroxobin and the gap is 2.0 with AUC 0.5744 and p 0.41; drop HBOC-201 and the gap is 1.0
with AUC 0.5515; drop both and the gap is 1.5 with AUC 0.5562; drop exendin(9-39) and the
gap is 2.0 with AUC 0.5858. The AUC condition fails in all four, so T does not predict in
any of them.

One split is worth naming because `results.json` records it: Set A against Set B(ii), the
nine named genes with no candidate in trials, gives AUC 0.7339 and two-sided p 0.0476,
while Set A against Set B(i), the 15 candidates that reached a trial, gives AUC 0.4702 and
p 0.7754. T separates approved drugs from gene-level proposals and does not separate them
from candidates that reached a trial and stopped. Both splits have fewer than 15 entries
on the B side, so under the rule each is underpowered and neither is called either way.
