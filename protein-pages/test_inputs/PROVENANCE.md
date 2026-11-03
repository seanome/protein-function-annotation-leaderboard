# Where the 51-protein test set comes from

The test set is 51 proteins of the new *Botryllus schlosseri* proteome (FUN ids, such as
`FUN000003_FUN000003`): 50 drawn at random and the Botryllus histocompatibility factor (BHF),
added by name. They were cut by `make_test_inputs.py` (run by `make_test_inputs.sbatch`,
Sherlock job 49793015 for the 50, then job 49841596 with BHF added, both 2026-10-27) out of
tool outputs that already existed. The only tool run for it is hmmscan. Job 49841596 kept
job 49793015's hmmscan rows and ran hmmscan only on the proteins that run had not scanned.

## Which proteins

The 50 come from the 7,216 FUN proteins that 2024-kmerseek-analysis notebook 281 searched
with phmmer, jackhmmer and MMseqs2 three rounds, keeping each hit's position on the
Botryllus protein. At E-value 0.001, counting the four programs that have output for them
(MMseqs2, phmmer, jackhmmer, kmerseek), 1,365 of the 7,216 had no hit, 280 had a hit from
one program, 220 from two and 5,351 from three. The script drew 16 with no hit, 16 with one
program (spread over the programs) and 18 with two or more, with random seed 20261027.
`selection.tsv` on Sherlock lists them.

BHF is `FUN008084_FUN008084`. Its 252 residues are identical to UniProt S5RI84,
"Histocompatibility factor" of *B. schlosseri* (Voskoboynik et al. 2013, PubMed 23888037),
checked on 2026-10-27 against `https://rest.uniprot.org/uniprotkb/S5RI84.fasta`. It is added
after the random draw, so the other 50 are the same as before.

hmmscan then adds hits to some of these, so the counts on the built site differ from the
counts at selection.

## Files

The final input folder is `$SCRATCH/protein-pages-final/inputs`.


| File | Source on Sherlock (`$SCRATCH` = `/scratch/users/olgabot`) | Change made |
|---|---|---|
| `proteins.fasta` | `botryllus-mhc/qpos/chunks/chunk_*.fasta` | the 51 only; trailing `*` removed |
| `mmseqs2_3rounds.m8`, `phmmer.domtblout`, `phmmer.out`, `jackhmmer.domtblout`, `jackhmmer.out` | rerun by `rerun_with_alignments.sbatch` (job 49844850) on the 51, with the same images, database and settings as notebook 281's `botryllus-mhc/281_tool_search_query_positions.sbatch` (`-E 10`; jackhmmer `-N 3`; MMseqs2 `-s 7 --num-iterations 3 -e 10`), keeping the alignments (HMMER `-o --notextw`; MMseqs2 `qaln,taln`) | none. The rerun gave the same rows as notebook 281's tables: phmmer 29,448, jackhmmer 66,196 and MMseqs2 5,513 rows, all in both. |
| `kmerseek.tsv.gz` | `botryllus-kmerseek-swissprot/results/matches/` (the whole proteome against all of Swiss-Prot 2026_03, 38 alphabet and k-mer size settings), cut by `extract_kmerseek.sbatch` (job 50499774, 2026-10-31) | regions of the 51 with region E-value at most 10, from the 35 settings whose E-value fit worked (34,052 rows; every one of the 51 has some); a subset of the columns. Replaces the first file, cut from the dark-set run's `dark-regions/`, which only held proteins with no other hit, at one setting. |
| `mmseqs2.m8` | `botryllus_mmseqs2_sp.m8` from Olga, 2026-10-31: Biomni's MMseqs2 run of the whole proteome against all of Swiss-Prot 2026_03 (Biomni says three rounds; the settings were not saved with the table), kept at E-value 0.001 or below | rows of the 51, matched by sequence MD5 to `full_inputs/proteome_lengths.tsv.gz`; columns put in build.py's order (query, target, fident, alnlen, qstart, qend, tstart, tend, evalue, bits) |
| `prostt5_foldseek.m8.gz` | `botryllus_prostt5_foldseek.m8` from Biomni, downloaded by Olga 2026-11-01: ProstT5 3Di letters, Foldseek 10, against 590,183 AlphaFold models of all Swiss-Prot (Foldseek's prebuilt database; sea-squirt entries not removed), E-value cutoff 10 | rows of the 51, matched by sequence MD5, by `../full_inputs/make_foldseek_tables.py`: 16,579 rows on 46 proteins (5 have none; Foldseek kept at most 1,000 hits per protein); columns put in build.py's order, with `alntmscore` empty (ProstT5 mode gives no TM-score) and the accession as `theader`, keeping an isoform's number (`AF-C3VPR6-2-F1-model_v6` gives `C3VPR6-2`) |
| `foldseek.m8.gz` | `botryllus_foldseek_struct.m8` from Biomni, given by Olga 2026-11-03: Foldseek on 3-D structures Biomni predicted for the 17,990 proteins of 200 aa or less, against the same 590,183 AlphaFold models of Swiss-Prot; 2,582,805 rows on 13,946 proteins. Its columns are query, target, fident, alnlen, qstart, qend, tstart, tend, evalue, bits and two Biomni has not named, which are left out | rows of the 51, matched by sequence MD5, by `../full_inputs/make_foldseek_tables.py ... struct`: 3,069 rows on 21 of the 25 drawn proteins of 200 aa or less (6 with a hit at E-value 0.001 or below) (the other 26 drawn proteins are longer, so were not run) |
| `structures/FUN*.pdb` except FUN008084 | `botryllus_structures_le200aa.zip` from Biomni, given by Olga 2026-11-03: one predicted structure per protein of 200 aa or less (17,990). The files do not name the predicting program; they start `PARENT N/A` and give confidence on a 0 to 1 scale | the 25 drawn proteins of 200 aa or less, matched by sequence MD5 and checked residue by residue, by `../structures/pick_drawn_structures.py`, which adds a first REMARK line saying where each came from. All 17,990 structures' residues match the proteome FASTA (checked 2026-11-03) |
| `hhsearch/*.hhr.gz` | `protein-pages-final/inputs/hhsearch/`, written by `hhsearch.sbatch` (array job 50673585, 2026-11-01, after `hhsuite_install.sbatch`, job 50673583): HH-suite 3.3.0 (BioContainers image `hhsuite:3.3.0--h503566f_15`); per protein, hhblits against UniRef30 2023_02 (`-n 3 -e 1e-3`) builds the alignment, then hhsearch against PDB70 (2022-03-13 build) with default hit lists | gzipped. 8,825 hits on the 51; 35 of the 51 have a hit with Prob of at least 50. BHF's alignment has 4 sequences and its best hit is Prob 18.8. Tasks took 8-29 min on 8 CPUs; the longest (task 30) peaked at 31.9 GB of its 32 GB, so a rerun should ask for 48 GB |
| `structures/FUN008084.pdb` | BHF's structure predicted by ESMFold2 (Biohub, MIT license; `esm` 3.4.2, weights `biohub/ESMFold2` and `biohub/ESMC-6B`), seed 0, on one H100 on Sherlock, 2026-11-02 (job 50776858, `structures/esmfold2_fold.sbatch`; 20 s to fold). Its first line, a REMARK added here, names the predictor | none. BHF's sequence is public (UniProt S5RI84). Neighbouring CA atoms a median 3.84 angstrom apart; mean pLDDT 76.7 of 100; pTM 0.25, so the predictor is unsure how the parts sit against each other. It replaces a file from Olga's Downloads (2026-11-01, predictor not recorded) whose CA atoms were a median 7.5 angstrom apart |
| `biomni_best.json` | Biomni's whole-proteome tables in `../full_inputs/` (DIAMOND best Swiss-Prot hit, jackhmmer from 24 immune seed proteins, FANTASIA, UniFIRE), the results that give no position on the protein | `make_biomni_best.py ../full_inputs data/proteins.fasta`, 2026-11-01: one entry per drawn protein, matched by sequence MD5; only accessions, names, scores and counts, no sequence. 15 of the 51 have a DIAMOND hit, 8 a jackhmmer hit from an immune seed, 37 a FANTASIA answer (3 below distance 0.15), 22 UniFIRE predictions |
| `interproscan.tsv` | InterProScan 5.59, all member databases, from nf-core/proteinannotator (olgabot/botryllus-mhc PR 13; `$SCRATCH/botryllus-annotation-full/results/proteinannotator/interproscan/*.xml`), turned into InterProScan's TSV layout by `interproscan_xml_to_tsv.py` (job 50500198) | rows of the 51, matched by sequence MD5; 767 rows on 38 proteins |
| `hmmscan.domtblout` | run here: hmmscan 3.4 (`-E 10 --domE 10`) against `$SCRATCH/data/pfam/Pfam-A.hmm` (30,134 models, Pfam-A 38.2), in three jobs that each kept the earlier rows and scanned only new proteins: 49793015 (the first 50 and their targets at E <= 0.001), 49842799 (BHF and targets up to E < 1) and 49845538 (targets up to E <= 10) | on the 51 proteins and the 25,729 Swiss-Prot proteins they match at E-value 10 or below, so build.py can name the Pfam families of a search hit's target |
| `hmmscan.out` | hmmscan's main output on the 51 proteins (job 49845108, same options), for the alignment of each Pfam model to the protein | none |
| `swissprot.fasta` | `botryllus-mhc-swissprot/db/swissprot/uniprot_sprot.fasta.gz` (release 2026_03) | only the 26,602 proteins any file names, plus the 4,225 kmerseek and MMseqs2 targets added on 2026-10-31 (from `botryllus-kmerseek-swissprot/input/uniprot_sprot_2026_03.fasta.gz`; their Pfam families are not in `hmmscan.domtblout`), with UniProt's whole headers (name, species, gene); the 3,486 ProstT5 + Foldseek targets it lacked were added on 2026-11-01 from UniProt REST (current release) by `fetch_missing_targets.py`; 4 not in release 2026_03 were taken from the dark-set reference (`invertebrate-dark-set/reference/minus_Ascidiacea/reference.fasta`), with no entry name; the 2,831 isoforms these targets name (for example `C3VPR6-2`) were added the same way on 2026-11-01, and 15 isoforms that current UniProt no longer lists (12 of P05661, and P46549-2, P90980-3, Q9BWV1-2; 17 hits on 5 proteins) could not be fetched |
| `protein_names.tsv` | written by hand: BHF's name, from UniProt S5RI84 | |
| `run_settings.tsv` | written by hand from the commands in the sbatch scripts above | |
| `kmerseek_alphabets.rs` | `seanome/kmerseek` `src/rust/alphabets.rs` at commit fc3ce72c (2026-08-25, the last change to that file on the checkout read) | none |

MMseqs2 three rounds, phmmer and jackhmmer searched Swiss-Prot 2025_01 without its 120 sea-squirt (Ascidiacea)
entries, 572,850 proteins, the reference of the dark-set run (`invertebrate-dark-set/reference/minus_Ascidiacea`).
kmerseek, Biomni's MMseqs2 and DIAMOND searched all of Swiss-Prot 2026_03. The table below says how each was checked.

## The copy in this repository, `test_inputs/data/`

`cut_test_data.py` made `data/` from the final input folder, so anyone can build the test
site without Sherlock (`python build.py test_inputs/data --out site`). It keeps what build.py
reads and gzips the files still over 20 MB:

| File | Cut |
|---|---|
| `hmmscan.domtblout.gz` | every row of the 51 proteins; for the Swiss-Prot proteins, only domains with i-Evalue at most 0.001, the only ones build.py uses for a target's families (152,140 of 1,182,305 rows) |
| `jackhmmer.out.gz` | the last of the three rounds of each query, the only one build.py reads |
| `phmmer.out.gz` | none, gzipped |
| every other file | none |

One file in `data/` was made here, not cut from Sherlock:

| File | Source | How |
|---|---|---|
| `target_pfam_interpro.tsv` | InterPro REST API, `https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/<accession>/` for a protein and `https://www.ebi.ac.uk/interpro/api/protein/uniprot/<accession>/?isoforms=<isoform>` for an isoform (P12345-2), one request per protein, at most 3 at a time; HTTP 204 means no Pfam match. Read on 2026-11-01 (InterPro's release of that day). | `fetch_target_pfam.py test_inputs/data`, run on a laptop. It lists the Swiss-Prot proteins that a drawn hit (E-value at most 10) on the 51 pages points to and that have no domain in `hmmscan.domtblout.gz` or `interproscan.tsv`. Those two cover few targets: the hmmscan rows of Swiss-Prot proteins were cut to i-Evalue 0.001, and the 4,225 kmerseek and MMseqs2 targets and 3,486 ProstT5 + Foldseek targets added later were never scanned. 10,660 proteins were missing (2,846 of them isoforms that ProstT5 + Foldseek matched). InterPro gives at least one Pfam match for 10,113 and none for 547; 30,601 rows, one per Pfam match, with its start and end on the protein. After it, no drawn target on the 51 pages lacks domain data. |

The site built from `data/` is the same, page for page, as the site built from the full
folder (checked on 2026-10-27 by comparing the pages with the build time and asset version
removed).

## Not run for these proteins

No output exists on the FUN proteome for Foldseek on predicted 3-D structures, and the proteome has no predicted structures (`plddt.tsv`) or expression table.
Those tracks read "not run". ProstT5 + Foldseek was run (`prostt5_foldseek.m8.gz` above).

## Which Swiss-Prot each search used

Checked on Sherlock, 2026-11-02, from the `reldate.txt` beside each database and by counting
proteins. The pages show this in a table under each drawing (`SWISSPROT_USED` in `build.py`).

| Programs | Database | Proteins | How it was checked |
|---|---|---|---|
| kmerseek; the sequences in `swissprot.fasta` | Swiss-Prot release 2026_03 (2 September 2026) | 575,748 | `botryllus-mhc-swissprot/db/swissprot/reldate.txt`; headers of `uniprot_sprot_2026_03.fasta.gz` counted |
| MMseqs2 three rounds, phmmer, jackhmmer | `invertebrate-dark-set/reference/minus_Ascidiacea/reference.fasta`: Swiss-Prot release 2025_01 without its 120 sea-squirt (Ascidiacea) proteins, built 2026-09-08 | 572,850 | Biomni matched the count to 2025_01's release notes (572,970 entries); checked 2026-11-03 on Sherlock against `uniprot_sprot.fasta.gz` of `previous_releases/release-2025_01`: every one of the 572,850 accessions is in 2025_01, and 2025_01 has exactly 120 more |
| MMseqs2 (Biomni) | Swiss-Prot release 2026_03, all of it | 575,748 | Biomni's note of 2026-11-03 (`full_inputs/biomni_databases_2026-11-03.md`): downloaded from `current_release` on 2026-10-31, which was 2026_03; its hit table has the same 3,566,205 rows and 22,043 proteins as `botryllus_mmseqs2_sp.m8` |
| DIAMOND (Biomni) | Swiss-Prot release 2026_03, all of it | 575,748 | Biomni's note of 2026-11-03: the same file as its MMseqs2 run. 533,805 is the number of Swiss-Prot proteins among its hits (E-value 0.001 or below, at most 25 per protein), not the database size, as `biomni_report_v3.md` had it |
| ProstT5 + Foldseek (Biomni) | AlphaFold models of Swiss-Prot, Foldseek's prebuilt AlphaFold/Swiss-Prot database, downloaded 2026-10-31 | 590,183 | Biomni's note of 2026-11-03: "Target database size: 590183" in the search log. More than 575,748 because AlphaFold keeps models of retired and merged accessions. `biomni_report_v3.md`'s 570,421 was wrong |
