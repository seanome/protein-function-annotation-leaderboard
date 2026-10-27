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

| File | Source on Sherlock (`$SCRATCH` = `/scratch/users/olgabot`) | Change made |
|---|---|---|
| `proteins.fasta` | `botryllus-mhc/qpos/chunks/chunk_*.fasta` | the 51 only; trailing `*` removed |
| `mmseqs2_3rounds.m8` | `botryllus-mhc/qpos/out/mmseqs2.chunk_*.tsv.gz` (`-s 7 --num-iterations 3 -e 10`) | columns put in build.py's order. That run did not write `fident`, `alnlen` or `theader`, so those columns are empty and the tooltip says "not in the file". |
| `phmmer.domtblout` | `botryllus-mhc/qpos/out/phmmer.chunk_*.domtblout.gz` (`-E 10`) | rows of the 51 only |
| `jackhmmer.domtblout` | `botryllus-mhc/qpos/out/jackhmmer.chunk_*.domtblout.gz` (`-N 3 -E 10`) | rows of the 51 only, all rounds (build.py keeps the last) |
| `kmerseek.tsv` | `invertebrate-dark-set/results-0.4/dark-regions/botryllus_dark_regions.csv.gz` | rows of the 51 with `moltype` hp_lehninger2 and `ksize` 17 only (the file also holds 5 other settings); commas changed to tabs; column names kept. This file holds regions only for proteins the dark-set run found no other hit for, so most of the 50 have no kmerseek row. It has no GO column. |
| `hmmscan.domtblout` | run here: hmmscan 3.4 (`-E 10 --domE 10`) against `$SCRATCH/data/pfam/Pfam-A.hmm` (30,134 models, Pfam-A 38.2) | on the 51 proteins and on the 16,262 Swiss-Prot proteins they match at E-value 0.001 or below, so build.py can name the Pfam family of a search hit's target |
| `swissprot.fasta` | `botryllus-mhc-swissprot/db/swissprot/uniprot_sprot.fasta.gz` (release 2026_03) | only the 26,596 proteins any file names; 4 not in release 2026_03 were taken from the dark-set reference (`invertebrate-dark-set/reference/minus_Ascidiacea/reference.fasta`), with no entry name |
| `run_settings.tsv` | written by hand from the commands in the sbatch scripts above | |
| `kmerseek_alphabets.rs` | `seanome/kmerseek` `src/rust/alphabets.rs` at commit fc3ce72c (2026-08-25, the last change to that file on the checkout read) | none |

All searches above were against Swiss-Prot without its 120 sea-squirt (Ascidiacea) entries,
572,850 proteins, the reference of the dark-set run (`invertebrate-dark-set/reference/minus_Ascidiacea`).

## Not run for these proteins

No output exists on the FUN proteome for MMseqs2 one round, Foldseek, ProstT5 + Foldseek,
hhsearch or InterProScan, and the proteome has no predicted structures (`plddt.tsv`) or
expression table. Those tracks read "not run". An older Botryllus gene set (ids like
`g21357.t1`) was searched with MMseqs2 one round and ProstT5 + Foldseek, but it is a
different set of gene models and is not used.
