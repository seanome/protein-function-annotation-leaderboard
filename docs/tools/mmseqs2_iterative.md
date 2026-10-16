# MMseqs2, three rounds

A profile of a protein family, built by MMseqs2 from its own hits over three rounds of search, as jackhmmer does.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.241 (n = 563), 60–100% 0.300 (n = 1,289). If no, 20–30% 0.337 (n = 93), 30–40% 0.209 (n = 165). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting sensitivity 7, 3 rounds; subsets with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.245 (n = 963), 121–250 aa 0.319 (n = 504), 251 aa or more 0.462 (n = 277). If no, 16–30 aa 0.028 (n = 271), 31–60 aa 0.085 (n = 420). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting sensitivity 7, 3 rounds.)
3. **Can you spend about 0.0787 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" below.

## Best Fmax and cost

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is sensitivity 7, 3 rounds (`s7`): Fmax 0.324 over 2,435 true instances. Searching costs about $0.00351 per 1,000 queries (0.0787 core-hours), for each relative searched.

## Where it is strong and weak

Each line compares this setting with the best setting of any other tool on the same subset of the proteins, so the other tool's setting can change from line to line. Only subsets with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- identity 30–40% to the closest same-family domain: Fmax 0.209 (n = 165); best other: phmmer, default options, 0.199.
- AlphaFold confidence (pLDDT) 90–100: Fmax 0.484 (n = 186); best other: MMseqs2, sensitivity 7, 0.480.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.028 (n = 271); best other: kmerseek, wwmj5, k = 11, 0.330.
- regions 31–60 aa long: Fmax 0.085 (n = 420); best other: kmerseek, gbmr4, k = 21, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting sensitivity 7, 3 rounds: best feature type domain, Fmax 0.321 (n = 1,107); worst short motif, Fmax 0.004 (n = 190).

## Can its weak hits be trusted? {#weak-hits}

When this tool's best hit on a protein is weak, is it right? It was measured on the dark set of four species, each searched against Swiss-Prot without its own clade. Dark proteins have no phmmer, jackhmmer or MMseqs2 hit at E ≤ 1; placed proteins are the rest. A best hit is right when the target shares a Pfam family with the query and the hit covers a Pfam domain of that family on the query. n counts the best hits whose target has any Pfam family. Test half of the proteins only. Grey: fewer than 30 judged hits ([how it was measured](../methods.html#weak-hits)).

**Dark proteins**, best hit right:

| species | E ≤ 0.001 | 0.001 < E ≤ 1 | 1 < E ≤ 10 |
|---|---:|---:|---:|
| human | <span class="few">no judged hits (n = 0)</span> | <span class="few">no judged hits (n = 0)</span> | 5.9% (5 of 85; 95% CI 2.5–13.0%) |
| C. elegans | <span class="few">no judged hits (n = 0)</span> | <span class="few">no judged hits (n = 0)</span> | 4.2% (4 of 96; 95% CI 1.6–10.2%) |
| S. cerevisiae | <span class="few">no judged hits (n = 0)</span> | <span class="few">no judged hits (n = 0)</span> | 3.2% (1 of 31; 95% CI 0.6–16.2%) |
| B. burgdorferi | <span class="few">no judged hits (n = 0)</span> | <span class="few">no judged hits (n = 0)</span> | <span class="few">0.0% (0 of 13; 95% CI 0.0–22.8%)</span> |

**Placed proteins**, best hit right:

| species | E ≤ 0.001 | 0.001 < E ≤ 1 | 1 < E ≤ 10 |
|---|---:|---:|---:|
| human | 94.4% (7,941 of 8,414; 95% CI 93.9–94.9%) | 53.3% (57 of 107; 95% CI 43.9–62.4%) | 12.9% (18 of 140; 95% CI 8.3–19.4%) |
| C. elegans | 92.4% (4,471 of 4,840; 95% CI 91.6–93.1%) | 30.1% (34 of 113; 95% CI 22.4–39.1%) | 15.4% (36 of 234; 95% CI 11.3–20.6%) |
| S. cerevisiae | 92.3% (2,018 of 2,187; 95% CI 91.1–93.3%) | 46.5% (20 of 43; 95% CI 32.5–61.1%) | 21.3% (10 of 47; 95% CI 12.0–34.9%) |
| B. burgdorferi | 93.5% (301 of 322; 95% CI 90.2–95.7%) | <span class="few">10.0% (1 of 10; 95% CI 1.8–40.4%)</span> | <span class="few">6.2% (1 of 16; 95% CI 1.1–28.3%)</span> |

Setting in the dark-set run: MMseqs2, sensitivity 7, 3 rounds, E ≤ 10.

[PR 133](https://github.com/seanome/2024-kmerseek-analysis/pull/133)'s looser test counts a shared Pfam family or clan anywhere on the two proteins. By that test, on dark proteins at E 1–10:

- human 5.9% (5 of 85)
- C. elegans 13.5% (13 of 96)
- S. cerevisiae 3.2% (1 of 31)
- B. burgdorferi 0.0% (0 of 13)


## What it compares

A profile of a protein family, built by MMseqs2 from its own hits over three rounds of search, as jackhmmer does.

## What it needs

- sequence only: yes
- 3D structure prediction: no
- GPU: no

## Reference searched

The relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as an MMseqs2 database.

## How to read its score

Bit score: higher ranks first. The E-value cutoff was 10.

## Known failure modes {#failure-modes}

- At most 1,000 target proteins are kept per query (--max-seqs 1000). How often that cut hits: not measured.
- Each round can pull in proteins from a different family. How often: not measured.

## Cost per 1,000 queries

Setting `s7`, searching one relative (mean over the 15 relatives it searched): 0.0787 core-hours reserved per 1,000 queries (16 cores for the task's wall time), measured on Sherlock. Peak memory 1.1 GB. That is about $0.00351 at $0.044625 per CPU-hour.
Building its databases took 0.00889 core-hours over 7 databases. A database is built once and reused, so this is not in the per-query cost.
This covers the 7 database builds the task records hold; the search used 15 relatives.

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Licence

- Licence: MIT ([source](https://github.com/soedinglab/MMseqs2))
- Commercial use allowed: Yes.

## Settings in the table

- sensitivity 7, 3 rounds (`s7`): sensitivity 7 (-s 7), three rounds

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

Tool id `mmseqs2_iterative`. Numbers from seanome/2024-kmerseek-analysis commit `4a4b1171`.

