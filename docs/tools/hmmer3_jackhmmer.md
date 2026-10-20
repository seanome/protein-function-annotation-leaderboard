# jackhmmer

Tool id `hmmer3_jackhmmer`. Numbers from seanome/2024-kmerseek-analysis commit `7eca428b`.

## What it compares

A profile of a protein family: it searches with one protein, builds a profile from the hits, and searches again with that profile, three rounds in all.

## What it needs

- sequence only: yes
- 3D structure prediction: no
- GPU: no

## Reference searched

The relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as a FASTA file.

## How to read its score

Domain bit score: higher ranks first. The E-value cutoff was 10.

## Strong and weak strata

Under the Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, its best setting is `n3`: Fmax 0.290 over 2,435 true instances. Below, each stratum compares this setting with the best setting of any other tool in that same stratum, so the other tool's setting can change from line to line. Only strata with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- regions 121–250 aa long: Fmax 0.333 (n = 504); best other: phmmer `default`, 0.329.
- no same-family domain in the relative: Fmax 0.023 (n = 321); best other: HH-suite (hhsearch) `single_seq`, 0.032.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.029 (n = 271); best other: kmerseek `wwmj5_k11_lcFalse`, 0.330.
- regions 31–60 aa long: Fmax 0.066 (n = 420); best other: kmerseek `gbmr4_k21_lcFalse`, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting `n3`: best feature type ACT_SITE (active site residue), Fmax 0.287 (n = 176); worst MOTIF (short motif), Fmax 0.002 (n = 190).

## Known failure modes {#failure-modes}

- Each round can pull in proteins from a different family, which then pull in more. How often that happened here: not measured.

## Cost per 1,000 queries

Setting `n3`, searching one relative (mean over the 15 relatives it searched): 2.45 core-hours reserved per 1,000 queries (8 cores for the task's wall time), measured on Sherlock. Peak memory 5.2 GB. That is about $0.109 at $0.044625 per CPU-hour.
Building its databases for setting `n3`: not measured (no task record holds a database build for this tool and setting).

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Licence

- Licence: BSD 3-clause (HMMER) ([source](https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE))
- Commercial use allowed: Yes.

## Settings in the table

- `n3`: three rounds of search (-N 3)

## Version, container, command

- Version: HMMER 3.4
- Container: `quay.io/biocontainers/hmmer:3.4--hb6cb901_4 (sha256:7a2b317b8d2fd3650b4924a8482cddeb940d4a0746c6a1501ff03ac1b7439e0c)`

```bash
jackhmmer \
    -N 3 \
    --domtblout /dev/stdout \
    --tblout /dev/null \
    -o /dev/stderr \
    --noali \
    -E 10.0 \
    --cpu 8 \
    UP000005640_9606.fasta UP000000437_7955.fasta
```

Query: 998 human proteins. Target: the zebrafish proteome.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.217 (n = 563), 60–100% 0.264 (n = 1,289). If no, 20–30% 0.327 (n = 93), 30–40% 0.175 (n = 165). (Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, setting `n3`; strata with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.210 (n = 963), 121–250 aa 0.333 (n = 504), 251 aa or more 0.427 (n = 277). If no, 16–30 aa 0.029 (n = 271), 31–60 aa 0.066 (n = 420). (Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, setting `n3`.)
3. **Can you spend about 2.45 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" above.
