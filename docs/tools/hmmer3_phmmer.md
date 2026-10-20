# phmmer

Amino-acid letters, one query protein against each protein of the relative, scored with a substitution matrix the way BLAST does but with HMMER's model of insertions and deletions.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.252 (n = 563), 60–100% 0.308 (n = 1,289). If no, 20–30% 0.359 (n = 93), 30–40% 0.199 (n = 165). (Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, setting default options; subsets with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.250 (n = 963), 121–250 aa 0.329 (n = 504), 251 aa or more 0.477 (n = 277). If no, 16–30 aa 0.042 (n = 271), 31–60 aa 0.095 (n = 420). (Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, setting default options.)
3. **Can you spend about 0.417 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" below.

## Best Fmax and cost

Under the Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, its best setting is default options (`default`): Fmax 0.333 over 2,435 true instances. Searching costs about $0.0186 per 1,000 queries (0.417 core-hours), for each relative searched.

## Where it is strong and weak

Each line compares this setting with the best setting of any other tool on the same subset of the proteins, so the other tool's setting can change from line to line. Only subsets with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- identity 20–30% to the closest same-family domain: Fmax 0.359 (n = 93); best other: MMseqs2, three rounds, sensitivity 7, 3 rounds, 0.337.
- AlphaFold confidence (pLDDT) 50–70: Fmax 0.325 (n = 748); best other: MMseqs2, sensitivity 7, 0.310.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.042 (n = 271); best other: kmerseek, wwmj5, k = 11, 0.330.
- regions 31–60 aa long: Fmax 0.095 (n = 420); best other: kmerseek, gbmr4, k = 21, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting default options: best feature type domain, Fmax 0.304 (n = 1,107); worst short motif, Fmax 0.000 (n = 190).

## Can its weak hits be trusted? {#weak-hits}

When this tool's best hit on a protein is weak, is it right? It was measured on the dark set of four species, each searched against Swiss-Prot without its own clade. Dark proteins have no phmmer, jackhmmer or MMseqs2 hit at E ≤ 1; placed proteins are the rest. A best hit is right when the target shares a Pfam family with the query and the hit covers a Pfam domain of that family on the query. n counts the best hits whose target has any Pfam family. Test half of the proteins only. Grey: fewer than 30 judged hits ([how it was measured](../methods.html#weak-hits)).

**Dark proteins**, best hit right:

| species | E ≤ 0.001 | 0.001 < E ≤ 1 | 1 < E ≤ 10 |
|---|---:|---:|---:|
| human | <span class="few">no judged hits (n = 0)</span> | <span class="few">no judged hits (n = 0)</span> | 2.2% (9 of 401; 95% CI 1.2–4.2%) |
| C. elegans | <span class="few">no judged hits (n = 0)</span> | <span class="few">no judged hits (n = 0)</span> | 1.8% (13 of 729; 95% CI 1.0–3.0%) |
| S. cerevisiae | <span class="few">no judged hits (n = 0)</span> | <span class="few">no judged hits (n = 0)</span> | 3.0% (4 of 132; 95% CI 1.2–7.5%) |
| B. burgdorferi | <span class="few">no judged hits (n = 0)</span> | <span class="few">no judged hits (n = 0)</span> | 0.0% (0 of 65; 95% CI 0.0–5.6%) |

**Placed proteins**, best hit right:

| species | E ≤ 0.001 | 0.001 < E ≤ 1 | 1 < E ≤ 10 |
|---|---:|---:|---:|
| human | 96.4% (8,110 of 8,417; 95% CI 95.9–96.7%) | 19.6% (118 of 601; 95% CI 16.7–23.0%) | <span class="few">4.8% (1 of 21; 95% CI 0.8–22.7%)</span> |
| C. elegans | 94.5% (4,616 of 4,884; 95% CI 93.8–95.1%) | 14.9% (151 of 1,011; 95% CI 12.9–17.3%) | <span class="few">3.4% (1 of 29; 95% CI 0.6–17.2%)</span> |
| S. cerevisiae | 95.8% (2,097 of 2,189; 95% CI 94.9–96.6%) | 17.6% (36 of 204; 95% CI 13.0–23.5%) | <span class="few">0.0% (0 of 10; 95% CI 0.0–27.8%)</span> |
| B. burgdorferi | 98.1% (313 of 319; 95% CI 96.0–99.1%) | 4.6% (5 of 109; 95% CI 2.0–10.3%) | <span class="few">0.0% (0 of 3; 95% CI 0.0–56.1%)</span> |

Setting in the dark-set run: phmmer, default options, E ≤ 10.

[PR 133](https://github.com/seanome/2024-kmerseek-analysis/pull/133)'s looser test counts a shared Pfam family or clan anywhere on the two proteins. By that test, on dark proteins at E 1–10:

- human 4.0% (16 of 401)
- C. elegans 3.4% (25 of 729)
- S. cerevisiae 3.8% (5 of 132)
- B. burgdorferi 0.0% (0 of 65)


## What it compares

Amino-acid letters, one query protein against each protein of the relative, scored with a substitution matrix the way BLAST does but with HMMER's model of insertions and deletions.

## What it needs

- sequence only: yes
- 3D structure prediction: no
- GPU: no

## Reference searched

The relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as a FASTA file.

## How to read its score

Domain bit score: higher ranks first. The E-value cutoff was 10.

## Known failure modes {#failure-modes}

- Failure modes specific to phmmer: not measured in this run.

## Cost per 1,000 queries

Setting `default`, searching one relative (mean over the 15 relatives it searched): 0.417 core-hours reserved per 1,000 queries (8 cores for the task's wall time), measured on Sherlock. Peak memory 4.9 GB. That is about $0.0186 at $0.044625 per CPU-hour.
Building its databases for setting `default`: not measured (no task record holds a database build for this tool and setting).

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Licence

- Licence: BSD 3-clause (HMMER) ([source](https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE))
- Commercial use allowed: Yes.

## Settings in the table

- default options (`default`): phmmer's default options, E-value cutoff 10

## Version, container, command

- Version: HMMER 3.4
- Container: `quay.io/biocontainers/hmmer:3.4--hb6cb901_4 (sha256:7a2b317b8d2fd3650b4924a8482cddeb940d4a0746c6a1501ff03ac1b7439e0c)`

```bash
phmmer \
    --domtblout /dev/stdout \
    --tblout /dev/null \
    -o /dev/stderr \
    --noali \
    -E 10.0 \
    --cpu 8 \
    UP000005640_9606.fasta UP000000437_7955.fasta
```

Query: 998 human proteins. Target: the zebrafish proteome.

Tool id `hmmer3_phmmer`. Numbers from seanome/2024-kmerseek-analysis commit `7eca428b`.

