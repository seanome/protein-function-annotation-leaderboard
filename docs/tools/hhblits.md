# HH-suite (hhsearch)

Tool id `hhblits`. Numbers from seanome/2024-kmerseek-analysis commit `3aa8d649`.

## What it compares

Two profiles against each other. Here each profile was built from a single protein, not from a family alignment, because no UniRef30 database was given.

## What it needs

- sequence only: yes
- AlphaFold model: no
- GPU: no

## Reference searched

Single-protein profiles of the relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as an HH-suite database.

## How to read its score

Bit score from hhsearch: higher ranks first.

## Strong and weak strata

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is `single_seq`: Fmax 0.313 over 2,435 true instances. Below, each stratum compares this setting with the best setting of any other tool in that same stratum, so the other tool's setting can change from line to line. Only strata with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- no same-family domain in the relative: Fmax 0.032 (n = 321); best other: Foldseek `3di_aa`, 0.027.
- identity 40–60% to the closest same-family domain: Fmax 0.248 (n = 563); best other: phmmer `default`, 0.252.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.066 (n = 271); best other: kmerseek `wwmj5_k11_lcFalse`, 0.330.
- regions 31–60 aa long: Fmax 0.085 (n = 420); best other: kmerseek `gbmr4_k21_lcFalse`, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting `single_seq`: best feature type DOMAIN (domain), Fmax 0.294 (n = 1,107); worst MOTIF (short motif), Fmax 0.003 (n = 190).

## Known failure modes {#failure-modes}

- The run names this tool hhblits (its tool id), but it ran hhsearch on profiles of single proteins. HH-suite's usual use, with family alignments from UniRef30, was not run, so these numbers say nothing about it.

## Cost per 1,000 queries

Setting `single_seq`, searching zebrafish: 21 CPU-hours, measured on Sherlock. Peak memory 23.4 GB. That is about $0.935 at $0.044625 per CPU-hour.
This search completed 2 times; the runs took 21 to 21.6 CPU-hours per 1,000 queries, and the newest is shown.
Across the 15 relatives searched with this setting: 2.18 to 30.9 CPU-hours per 1,000 queries; peak memory up to 23.4 GB.

CPU-hours count the search only, not building the relative's database.

## Settings in the table

- `single_seq`: profiles built from single proteins, no family database

## Version, container, command

- Version: HH-suite 3.3.0
- Container: `quay.io/biocontainers/hhsuite:3.3.0--h503566f_15 (sha256:4bf9bb5229de18f522a94f4443c19fdcbb0f0cb0e6ea92f5390aa170bcb0a24f)`

```bash
ffindex_apply \
    human-38fdd87d_hhdb/db_a3m.ffdata human-38fdd87d_hhdb/db_a3m.ffindex \
    -d results.ffdata -i results.ffindex \
    -- hhsearch -i stdin -d zebrafish_hhdb/db -blasttab /dev/stdout -cpu 16 -v 0
```

One hhsearch per query protein, run through ffindex_apply.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.248 (n = 563), 60–100% 0.283 (n = 1,289). If no, 20–30% 0.321 (n = 93), 30–40% 0.190 (n = 165). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `single_seq`; strata with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.230 (n = 963), 121–250 aa 0.298 (n = 504), 251 aa or more 0.417 (n = 277). If no, 16–30 aa 0.066 (n = 271), 31–60 aa 0.085 (n = 420). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `single_seq`.)
3. **Can you spend about 21 CPU-hours per 1,000 proteins for each relative you search?** That is the measured cost against zebrafish; the range over relatives is under "Cost per 1,000 queries" above.
