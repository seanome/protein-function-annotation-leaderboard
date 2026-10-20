# MMseqs2

Amino-acid letters: it first keeps target proteins that share similar short words with the query, then aligns only those, which makes it much faster than aligning every pair.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.249 (n = 563), 60–100% 0.297 (n = 1,289). If no, 20–30% 0.325 (n = 93), 30–40% 0.188 (n = 165). (Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, setting sensitivity 7; subsets with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.243 (n = 963), 121–250 aa 0.305 (n = 504), 251 aa or more 0.462 (n = 277). If no, 16–30 aa 0.027 (n = 271), 31–60 aa 0.089 (n = 420). (Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, setting sensitivity 7.)
3. **Can you spend about 0.0252 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" below.

## Best Fmax and cost

Under the Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, its best setting is sensitivity 7 (`s7`): Fmax 0.327 over 2,435 true instances. Searching costs about $0.00113 per 1,000 queries (0.0252 core-hours), for each relative searched.

## Where it is strong and weak

Each line compares this setting with the best setting of any other tool on the same subset of the proteins, so the other tool's setting can change from line to line. Only subsets with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- identity 40–60% to the closest same-family domain: Fmax 0.249 (n = 563); best other: phmmer, default options, 0.252.
- AlphaFold confidence (pLDDT) 90–100: Fmax 0.480 (n = 186); best other: MMseqs2, three rounds, sensitivity 7, 3 rounds, 0.484.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.027 (n = 271); best other: kmerseek, wwmj5, k = 11, 0.330.
- regions 31–60 aa long: Fmax 0.089 (n = 420); best other: kmerseek, gbmr4, k = 21, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting sensitivity 7: best feature type active site, Fmax 0.312 (n = 176); worst short motif, Fmax 0.000 (n = 190).

## Can its weak hits be trusted? {#weak-hits}

<span class="few">Not run on the dark set.</span> The weak-hit test so far covers phmmer, jackhmmer, MMseqs2 (3 rounds) and kmerseek ([how it was measured](../methods.html#weak-hits)).


## What it compares

Amino-acid letters: it first keeps target proteins that share similar short words with the query, then aligns only those, which makes it much faster than aligning every pair.

## What it needs

- sequence only: yes
- 3D structure prediction: no
- GPU: no

## Reference searched

The relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as an MMseqs2 database.

## How to read its score

Bit score: higher ranks first. The E-value cutoff was 10.

## Known failure modes {#failure-modes}

- At most 1,000 target proteins are kept per query (--max-seqs 1000), so a query from a very large family can lose hits. How often: not measured.

## Cost per 1,000 queries

Setting `s7`, searching one relative (mean over the 15 relatives it searched): 0.0252 core-hours reserved per 1,000 queries (16 cores for the task's wall time), measured on Sherlock. Peak memory 1.1 GB. That is about $0.00113 at $0.044625 per CPU-hour.
Building its databases took 0.00889 core-hours over 7 databases. A database is built once and reused, so this is not in the per-query cost.
This covers the 7 database builds the task records hold; the search used 15 relatives.

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Licence

- Licence: MIT ([source](https://github.com/soedinglab/MMseqs2))
- Commercial use allowed: Yes.

## Settings in the table

- sensitivity 7 (`s7`): sensitivity 7 (-s 7), one round

## Version, container, command

- Version: MMseqs2 18.8cc5c
- Container: `quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0 (sha256:3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef)`

```bash
mmseqs search \
    human-38fdd87d_mmdb/db zebrafish_mmdb/db result mmseqs_tmp \
    --threads 16 \
    -s 7 \
    --max-seqs 1000 \
    -e 10.0
```

Followed by mmseqs convertalis to write the table of hits.

Tool id `mmseqs2_seqseq`. Numbers from seanome/2024-kmerseek-analysis commit `7eca428b`.

