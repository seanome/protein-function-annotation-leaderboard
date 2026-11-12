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
| `kmerseek.tsv` | `invertebrate-dark-set/results-0.4/dark-regions/botryllus_dark_regions.csv.gz` | rows of the 51 with `moltype` hp_lehninger2 and `ksize` 17 only (the file also holds 5 other settings); commas changed to tabs; column names kept. This file holds regions only for proteins the dark-set run found no other hit for, so most of the 50 have no kmerseek row. It has no GO column. |
| `hmmscan.domtblout` | run here: hmmscan 3.4 (`-E 10 --domE 10`) against `$SCRATCH/data/pfam/Pfam-A.hmm` (30,134 models, Pfam-A 38.2), in three jobs that each kept the earlier rows and scanned only new proteins: 49793015 (the first 50 and their targets at E <= 0.001), 49842799 (BHF and targets up to E < 1) and 49845538 (targets up to E <= 10) | on the 51 proteins and the 25,729 Swiss-Prot proteins they match at E-value 10 or below, so build.py can name the Pfam families of a search hit's target |
| `hmmscan.out` | hmmscan's main output on the 51 proteins (job 49845108, same options), for the alignment of each Pfam model to the protein | none |
| `swissprot.fasta` | `botryllus-mhc-swissprot/db/swissprot/uniprot_sprot.fasta.gz` (release 2026_03) | only the 26,602 proteins any file names, with UniProt's whole headers (name, species, gene); 4 not in release 2026_03 were taken from the dark-set reference (`invertebrate-dark-set/reference/minus_Ascidiacea/reference.fasta`), with no entry name |
| `protein_names.tsv` | written by hand: BHF's name, from UniProt S5RI84 | |
| `run_settings.tsv` | written by hand from the commands in the sbatch scripts above | |
| `kmerseek_alphabets.rs` | `seanome/kmerseek` `src/rust/alphabets.rs` at commit fc3ce72c (2026-08-25, the last change to that file on the checkout read) | none |

All searches above were against Swiss-Prot without its 120 sea-squirt (Ascidiacea) entries,
572,850 proteins, the reference of the dark-set run (`invertebrate-dark-set/reference/minus_Ascidiacea`).

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

The site built from `data/` is the same, page for page, as the site built from the full
folder (checked on 2026-10-27 by comparing the pages with the build time and asset version
removed).

## Not run for these proteins

No output exists on the FUN proteome for MMseqs2 one round, Foldseek, ProstT5 + Foldseek,
hhsearch or InterProScan, and the proteome has no predicted structures (`plddt.tsv`) or
expression table. Those tracks read "not run". An older Botryllus gene set (ids like
`g21357.t1`) was searched with MMseqs2 one round and ProstT5 + Foldseek, but it is a
different set of gene models and is not used.
