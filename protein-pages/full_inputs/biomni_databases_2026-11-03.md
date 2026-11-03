# Which Swiss-Prot each search used: Biomni's note, 2026-11-03

Biomni's answer, pasted by Olga on 2026-11-03, to the question of which Swiss-Prot release each
program searched. Kept as Biomni wrote it, apart from the layout. The 2025_01 match was checked
on Sherlock the same day by comparing accessions (`../test_inputs/PROVENANCE.md`); the other
statements are Biomni's.

| Program | Database searched | Part searched | Proteins |
|---|---|---|---|
| kmerseek | Swiss-Prot release 2026_03 | all of it | 575,748 |
| MMseqs2 three rounds, phmmer, jackhmmer | Swiss-Prot release 2025_01 (released 05-Feb-2025) | without its 120 sea-squirt (Ascidiacea) proteins | 572,850 |
| MMseqs2 (Biomni run) | Swiss-Prot release 2026_03 (released 02-Sep-2026; downloaded 2026-10-31 from `current_release`) | all of it, Ascidiacea included | 575,748 |
| DIAMOND (Biomni run) | Swiss-Prot release 2026_03, the same file as the MMseqs2 run | all of it (`--sensitive`, E ≤ 1e-3, top 25 per query) | 575,748 in the database; 533,805 distinct targets in the hit table |
| ProstT5 + Foldseek | Foldseek's prebuilt AlphaFold/Swiss-Prot database, downloaded 2026-10-31 | sea-squirt proteins kept | 590,183 |

- **2025_01:** its release notes give 572,970 Swiss-Prot entries; 572,970 − 120 = 572,850. Nearby
  releases have other counts: 2025_02 573,230, 2025_03 573,661, 2025_04 573,661, 2026_01 574,627,
  2026_02 575,503.
- **Biomni's runs used 2026_03:** `reldate.txt` of `current_release` gives 2026_03, 02-Sep-2026.
  The MMseqs2 three-round command searched the unmodified `uniprot_sprot.fasta`. Its hit table,
  `botryllus_mmseqs2_sp.m8`, has 3.57 million lines, 292,113 distinct Swiss-Prot targets, and
  22,043 of the 45,339 Botryllus proteins with at least one hit.
- **DIAMOND:** 575,748 is the database size; 533,805 is the number of distinct Swiss-Prot proteins
  among the hits at E ≤ 1e-3, top 25 per query.
- **Correction:** report v2 said 570,421 AlphaFold models. The database both Foldseek runs used
  has 590,183 (its lookup file, and "Target database size: 590183" in the search log). It is larger
  than 575,748 because the AlphaFold database keeps models of retired and merged accessions.
- **jackhmmer in Biomni's session** ran the other way round: Swiss-Prot immune seed proteins as
  queries against the Botryllus proteome. It is not the jackhmmer run against Swiss-Prot above.

Sources Biomni gave: the `relnotes.txt` of each release under
https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/ and
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt
