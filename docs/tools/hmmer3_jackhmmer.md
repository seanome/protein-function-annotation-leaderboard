# jackhmmer

Tool id `hmmer3_jackhmmer`. Numbers from seanome/2024-kmerseek-analysis commit `6af1255b`.

## What it compares

A profile of a protein family: it searches with one protein, builds a profile from the hits, and searches again with that profile, three rounds in all.

## What it needs

- sequence only: yes
- AlphaFold model: no
- GPU: no

## Reference searched

The relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as a FASTA file.

## How to read its score

Domain bit score: higher ranks first. The E-value cutoff was 10.

## Strong and weak strata

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is `n3`: Fmax 0.290 over 2,435 true instances. Below, each stratum is compared with the best other tool and setting in the same stratum. Only strata with at least 30 true instances are compared.

Where it compares best with the other tools (its Fmax, then the best other):

- regions 121–250 aa long: Fmax 0.333 (n = 504); best other: phmmer `default`, 0.329.
- no same-family domain in the relative: Fmax 0.023 (n = 321); best other: HH-suite (hhsearch) `single_seq`, 0.032.

Where it compares worst:

- regions 16–30 aa long: Fmax 0.029 (n = 271); best other: kmerseek `wwmj5_k11_lcFalse`, 0.330.
- regions 31–60 aa long: Fmax 0.066 (n = 420); best other: kmerseek `gbmr4_k21_lcFalse`, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting `n3`: best feature type ACT_SITE (active site residue), Fmax 0.287 (n = 176); worst MOTIF (short motif), Fmax 0.002 (n = 190).

## Known failure modes {#failure-modes}

- Each round can pull in proteins from a different family, which then pull in more. How often that happened here: not measured.

## Cost per 1,000 queries

Setting `n3`, searching zebrafish: 1.83 CPU-hours and 5.2 GB peak memory, measured on Sherlock (Stanford Sherlock, SLURM partition hns, CPU nodes, no GPU); about $0.0818 at $0.044625 per CPU-hour.
Across the 15 relatives searched with this setting: 0.0402 to 2.43 CPU-hours per 1,000 queries; peak memory up to 5.2 GB.

CPU-hours count the search only, not building the relative's database.

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

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.217, 60–100% 0.264. If no, 20–30% 0.327, 30–40% 0.175. (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `n3`; strata with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.210, 121–250 aa 0.333, 251 aa or more 0.427. If no, 16–30 aa 0.029, 31–60 aa 0.066. (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `n3`.)
3. **Can you spend about 1.83 CPU-hours per 1,000 proteins for each relative you search?** That is the measured search cost against zebrafish; the range over relatives is under "Cost per 1,000 queries" above.
