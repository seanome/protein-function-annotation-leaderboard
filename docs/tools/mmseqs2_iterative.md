# MMseqs2, three rounds

Tool id `mmseqs2_iterative`. Numbers from seanome/2024-kmerseek-analysis commit `548ec5b8`.

## What it compares

A profile of a protein family, built by MMseqs2 from its own hits over three rounds of search, as jackhmmer does.

## What it needs

- sequence only: yes
- AlphaFold model: no
- GPU: no

## Reference searched

The relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as an MMseqs2 database.

## How to read its score

Bit score: higher ranks first. The E-value cutoff was 10.

## Strong and weak strata

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is `s7`: Fmax 0.324 over 2,435 true instances. Below, each stratum compares this setting with the best setting of any other tool in that same stratum, so the other tool's setting can change from line to line. Only strata with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- identity 30–40% to the closest same-family domain: Fmax 0.209 (n = 165); best other: phmmer `default`, 0.199.
- AlphaFold confidence (pLDDT) 90–100: Fmax 0.484 (n = 186); best other: MMseqs2 `s7`, 0.480.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.028 (n = 271); best other: kmerseek `wwmj5_k11_lcFalse`, 0.330.
- regions 31–60 aa long: Fmax 0.085 (n = 420); best other: kmerseek `gbmr4_k21_lcFalse`, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting `s7`: best feature type DOMAIN (domain), Fmax 0.321 (n = 1,107); worst MOTIF (short motif), Fmax 0.004 (n = 190).

## Known failure modes {#failure-modes}

- At most 1,000 target proteins are kept per query (--max-seqs 1000). How often that cut hits: not measured.
- Each round can pull in proteins from a different family. How often: not measured.

## Cost per 1,000 queries

Setting `s7`, searching zebrafish: 0.0743 CPU-hours, measured on Sherlock. Peak memory 1.1 GB. That is about $0.00332 at $0.044625 per CPU-hour.
This search completed 2 times; the runs took 0.0743 to 0.119 CPU-hours per 1,000 queries, and the newest is shown.
Across the 15 relatives searched with this setting: 0.0236 to 0.0769 CPU-hours per 1,000 queries; peak memory up to 1.1 GB.

CPU-hours count the search only, not building the relative's database.

## Settings in the table

- `s7`: sensitivity 7 (-s 7), three rounds

## Version, container, command

- Version: MMseqs2 18.8cc5c
- Container: `quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0 (sha256:3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef)`

```bash
mmseqs search \
    human-38fdd87d_mmdb/db zebrafish_mmdb/db result mmseqs_tmp \
    --threads 16 \
    -s 7 \
    --num-iterations 3 \
    --max-seqs 1000 \
    -e 10.0
```

Followed by mmseqs convertalis to write the table of hits.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.241 (n = 563), 60–100% 0.300 (n = 1,289). If no, 20–30% 0.337 (n = 93), 30–40% 0.209 (n = 165). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `s7`; strata with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.245 (n = 963), 121–250 aa 0.319 (n = 504), 251 aa or more 0.462 (n = 277). If no, 16–30 aa 0.028 (n = 271), 31–60 aa 0.085 (n = 420). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `s7`.)
3. **Can you spend about 0.0743 CPU-hours per 1,000 proteins for each relative you search?** That is the measured cost against zebrafish; the range over relatives is under "Cost per 1,000 queries" above.
