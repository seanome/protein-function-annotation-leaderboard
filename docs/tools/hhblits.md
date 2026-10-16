# HH-suite (hhsearch)

Two profiles against each other.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.248 (n = 563), 60–100% 0.283 (n = 1,289). If no, 20–30% 0.321 (n = 93), 30–40% 0.190 (n = 165). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting single-protein profiles; subsets with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.230 (n = 963), 121–250 aa 0.298 (n = 504), 251 aa or more 0.417 (n = 277). If no, 16–30 aa 0.066 (n = 271), 31–60 aa 0.085 (n = 420). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting single-protein profiles.)
3. **Can you spend about 24.8 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" below.

## Best Fmax and cost

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is single-protein profiles (`single_seq`): Fmax 0.313 over 2,435 true instances. Searching costs about $1.11 per 1,000 queries (24.8 core-hours), for each relative searched.

## Known failure modes

- Aligner calls run past the domain, because a call spans the whole aligned stretch. Over correct calls, they start a median 5 to 21 residues early. Over all calls, 45% to 60% start more than 100 residues early. The number: median 5 to 21 residues early; 45% to 60% more than 100 residues early. Status: open. What to do: Read a call's ends as the alignment's ends, not the domain's. To score where a domain sits, place the target's domain through the alignment. Source: [seanome/2024-kmerseek-analysis PR 146](https://github.com/seanome/2024-kmerseek-analysis/pull/146); [PR 3, comment 6193138553](https://github.com/seanome/protein-function-annotation-leaderboard/pull/3#issuecomment-6193138553).

## Where it is strong and weak

Each line compares this setting with the best setting of any other tool on the same subset of the proteins, so the other tool's setting can change from line to line. Only subsets with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- no same-family domain in the relative: Fmax 0.032 (n = 321); best other: Foldseek, shape letters + amino acids, 0.027.
- identity 40–60% to the closest same-family domain: Fmax 0.248 (n = 563); best other: phmmer, default options, 0.252.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.066 (n = 271); best other: kmerseek, wwmj5, k = 11, 0.330.
- regions 31–60 aa long: Fmax 0.085 (n = 420); best other: kmerseek, gbmr4, k = 21, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting single-protein profiles: best feature type domain, Fmax 0.294 (n = 1,107); worst short motif, Fmax 0.003 (n = 190).


## What it compares

Two profiles against each other. Here each profile was built from a single protein, not from a family alignment, because no UniRef30 database was given.

## What it needs

- sequence only: yes
- 3D structure prediction: no
- GPU: no

## Reference searched

Single-protein profiles of the relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as an HH-suite database.

## How to read its score

Bit score from hhsearch: higher ranks first.

## Caveats of how this run was set up {#caveats}

- The run names this tool hhblits (its tool id), but it ran hhsearch on profiles of single proteins. HH-suite's usual use, with family alignments from UniRef30, was not run, so these numbers say nothing about it.

## Cost per 1,000 queries

Setting `single_seq`, searching one relative (mean over the 15 relatives it searched): 24.8 core-hours reserved per 1,000 queries (16 cores for the task's wall time), measured on Sherlock. Peak memory 23.4 GB. That is about $1.11 at $0.044625 per CPU-hour.
Building its databases took 30.9 core-hours over 7 databases. A database is built once and reused, so this is not in the per-query cost.
This covers the 7 database builds the task records hold; the search used 15 relatives.

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Licence

- Licence: GPL-3.0 ([source](https://github.com/soedinglab/hh-suite))
- Commercial use allowed: Yes. If you distribute a changed version, you must release its source under GPL-3.0 too.

## Settings in the table

- single-protein profiles (`single_seq`): profiles built from single proteins, no family database

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

Tool id `hhblits`. Numbers from seanome/2024-kmerseek-analysis commit `4a4b1171`.

