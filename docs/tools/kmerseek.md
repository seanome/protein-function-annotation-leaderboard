# kmerseek

Short words of k letters (k-mers) shared between two proteins.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.005 (n = 563), 60–100% 0.297 (n = 1,289). If no, 20–30% 0.000 (n = 93), 30–40% 0.000 (n = 165). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting gbmr4, k = 21; subsets with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.208 (n = 963), 121–250 aa 0.098 (n = 504), 251 aa or more 0.050 (n = 277). If no, 16–30 aa 0.281 (n = 271), 31–60 aa 0.232 (n = 420). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting gbmr4, k = 21.)
3. **Can you spend about 0.0552 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" below.

## Best Fmax and cost

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is gbmr4, k = 21 (`gbmr4_k21_lcFalse`): Fmax 0.192 over 2,435 true instances. Searching costs about $0.00246 per 1,000 queries (0.0552 core-hours), for each relative searched.

## Where it is strong and weak

Each line compares this setting with the best setting of any other tool on the same subset of the proteins, so the other tool's setting can change from line to line. Only subsets with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- regions 16–30 aa long: Fmax 0.281 (n = 271); best other: HH-suite (hhsearch), single-protein profiles, 0.066.
- regions 31–60 aa long: Fmax 0.232 (n = 420); best other: phmmer, default options, 0.095.

Where it is furthest behind the best other tool:

- regions 251 aa or more long: Fmax 0.050 (n = 277); best other: phmmer, default options, 0.477.
- identity 20–30% to the closest same-family domain: Fmax 0.000 (n = 93); best other: phmmer, default options, 0.359.

Swiss-Prot features, raw calls, zebrafish, setting wwmj5, k = 11: best feature type membrane-spanning, Fmax 0.362 (n = 1,124); worst coiled coil, Fmax 0.000 (n = 145).

## Can its weak hits be trusted? {#weak-hits}

When this tool's best hit on a protein is weak, is it right? It was measured on the dark set of four species, each searched against Swiss-Prot without its own clade. Dark proteins have no phmmer, jackhmmer or MMseqs2 hit at E ≤ 1; placed proteins are the rest. A best hit is right when the target shares a Pfam family with the query and the hit covers a Pfam domain of that family on the query. n counts the best hits whose target has any Pfam family. Test half of the proteins only. Grey: fewer than 30 judged hits ([how it was measured](../methods.html#weak-hits)).

**Dark proteins**, best hit right:

| species | E ≤ 0.001 | 0.001 < E ≤ 1 | 1 < E ≤ 10 |
|---|---:|---:|---:|
| human | <span class="few">0.0% (0 of 6; 95% CI 0.0–39.0%)</span> | 0.0% (0 of 75; 95% CI 0.0–4.9%) | 0.0% (0 of 292; 95% CI 0.0–1.3%) |
| C. elegans | <span class="few">no judged hits (n = 0)</span> | 0.0% (0 of 33; 95% CI 0.0–10.4%) | 0.0% (0 of 38; 95% CI 0.0–9.2%) |
| S. cerevisiae | <span class="few">no judged hits (n = 0)</span> | <span class="few">0.0% (0 of 15; 95% CI 0.0–20.4%)</span> | 0.0% (0 of 106; 95% CI 0.0–3.5%) |
| B. burgdorferi | <span class="few">no judged hits (n = 0)</span> | <span class="few">0.0% (0 of 3; 95% CI 0.0–56.1%)</span> | 0.0% (0 of 51; 95% CI 0.0–7.0%) |

**Placed proteins**, best hit right:

| species | E ≤ 0.001 | 0.001 < E ≤ 1 | 1 < E ≤ 10 |
|---|---:|---:|---:|
| human | 94.3% (6,010 of 6,374; 95% CI 93.7–94.8%) | 57.5% (541 of 941; 95% CI 54.3–60.6%) | 9.5% (106 of 1,120; 95% CI 7.9–11.3%) |
| C. elegans | <span class="few">85.7% (6 of 7; 95% CI 48.7–97.4%)</span> | 15.9% (31 of 195; 95% CI 11.4–21.7%) | 2.0% (7 of 349; 95% CI 1.0–4.1%) |
| S. cerevisiae | 95.3% (832 of 873; 95% CI 93.7–96.5%) | 49.6% (236 of 476; 95% CI 45.1–54.1%) | 9.5% (78 of 823; 95% CI 7.7–11.7%) |
| B. burgdorferi | 95.3% (203 of 213; 95% CI 91.6–97.4%) | 50.7% (34 of 67; 95% CI 39.1–62.3%) | 9.9% (13 of 131; 95% CI 5.9–16.2%) |

The kmerseek setting for each species was picked on the other half of the proteins, as the one with the highest share right at E 1–10 on dark proteins:

- human: alphabet wass14, k = 6, scaled 1
- C. elegans: alphabet hp_lehninger_hpc3, k = 25, scaled 10
- S. cerevisiae: alphabet hp_pbotc_1st_ed2, k = 23, scaled 1
- B. burgdorferi: alphabet uniprot18, k = 5, scaled 1

The pick is close to arbitrary: on the tune half no setting was ahead by more than two right hits, and in two species every setting had none, so the setting with the most judged hits was taken ([rule and numbers](https://github.com/seanome/2024-kmerseek-analysis/blob/olgabot/dark-set-weak-hit-accuracy/nextflow-runs/invertebrate-dark-set/weak-hits/README.md)).

[PR 133](https://github.com/seanome/2024-kmerseek-analysis/pull/133)'s looser test counts a shared Pfam family or clan anywhere on the two proteins. By that test, on dark proteins at E 1–10:

- human 0.7% (2 of 292)
- C. elegans 0.0% (0 of 38)
- S. cerevisiae 0.9% (1 of 106)
- B. burgdorferi 0.0% (0 of 51)


## What it compares

Short words of k letters (k-mers) shared between two proteins. First, each amino acid is replaced by a letter from a smaller alphabet, for example H for hydrophobic and P for polar. Two proteins can then match where their residues differ but their pattern of residue types is the same.

## What it needs

- sequence only: yes
- 3D structure prediction: no
- GPU: no

## Reference searched

An index of the relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), built once per alphabet and k-mer length.

## How to read its score

region_enrichment: the number of k-mers the region shares with the target divided by the number expected by chance (fold enrichment). Higher ranks first. A region was kept when either test passed: its own score was at least 1.3, or the whole query protein's p-value was at most 0.05.

## Known failure modes {#failure-modes}

- A feature shorter than k residues contains no complete k-mer, so it cannot be found. Figure 2 above shows how Fmax changes with region length.
- The task script records a kmerseek search that exits with an error as a search with no hits, so a crash would look like a miss. How often that happened in this run: not measured.
- Setting hp_thomas_dill2 k33 was scored, but no trace row holds its search task, so its cost is not measured, and the scoring task for it was not found either.

## Cost per 1,000 queries

Setting `gbmr4_k21_lcFalse`, searching one relative (mean over the 15 relatives it searched): 0.0552 core-hours reserved per 1,000 queries (4 cores for the task's wall time), measured on Sherlock. Peak memory 5.2 GB. That is about $0.00246 at $0.044625 per CPU-hour.
Building its databases took 0.147 core-hours over 6 databases. A database is built once and reused, so this is not in the per-query cost.
This covers the 6 database builds the task records hold; the search used 15 relatives.

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Licence

- Licence: MIT ([source](https://github.com/seanome/kmerseek))
- Commercial use allowed: Yes.

## Settings in the table

- hp_thomas_dill2, k = 24 (`hp_thomas_dill2_k24_lcFalse`): alphabet hp_thomas_dill2 (2 letters, hydrophobic ACFILMVWY and polar DEGHKNPQRST (Thomas and Dill 1996)), k = 24, low-complexity k-mers kept
- hp_thomas_dill2, k = 27 (`hp_thomas_dill2_k27_lcFalse`): alphabet hp_thomas_dill2 (2 letters, hydrophobic ACFILMVWY and polar DEGHKNPQRST (Thomas and Dill 1996)), k = 27, low-complexity k-mers kept
- hp_thomas_dill2, k = 30 (`hp_thomas_dill2_k30_lcFalse`): alphabet hp_thomas_dill2 (2 letters, hydrophobic ACFILMVWY and polar DEGHKNPQRST (Thomas and Dill 1996)), k = 30, low-complexity k-mers kept
- hp_thomas_dill2, k = 33 (`hp_thomas_dill2_k33_lcFalse`): alphabet hp_thomas_dill2 (2 letters, hydrophobic ACFILMVWY and polar DEGHKNPQRST (Thomas and Dill 1996)), k = 33, low-complexity k-mers kept
- gbmr4, k = 15 (`gbmr4_k15_lcFalse`): alphabet gbmr4 (4 letters: ADKERNTSQ, YFLIVMCWH, G, P (Solis and Rackovsky 2000)), k = 15, low-complexity k-mers kept
- gbmr4, k = 17 (`gbmr4_k17_lcFalse`): alphabet gbmr4 (4 letters: ADKERNTSQ, YFLIVMCWH, G, P (Solis and Rackovsky 2000)), k = 17, low-complexity k-mers kept
- gbmr4, k = 19 (`gbmr4_k19_lcFalse`): alphabet gbmr4 (4 letters: ADKERNTSQ, YFLIVMCWH, G, P (Solis and Rackovsky 2000)), k = 19, low-complexity k-mers kept
- gbmr4, k = 21 (`gbmr4_k21_lcFalse`): alphabet gbmr4 (4 letters: ADKERNTSQ, YFLIVMCWH, G, P (Solis and Rackovsky 2000)), k = 21, low-complexity k-mers kept
- wwmj5, k = 11 (`wwmj5_k11_lcFalse`): alphabet wwmj5 (5 letters: CMFILVWY, ATH, GP, DE, SNQRK (Wang and Wang 1999)), k = 11, low-complexity k-mers kept
- wwmj5, k = 12 (`wwmj5_k12_lcFalse`): alphabet wwmj5 (5 letters: CMFILVWY, ATH, GP, DE, SNQRK (Wang and Wang 1999)), k = 12, low-complexity k-mers kept
- wwmj5, k = 14 (`wwmj5_k14_lcFalse`): alphabet wwmj5 (5 letters: CMFILVWY, ATH, GP, DE, SNQRK (Wang and Wang 1999)), k = 14, low-complexity k-mers kept
- wwmj5, k = 15 (`wwmj5_k15_lcFalse`): alphabet wwmj5 (5 letters: CMFILVWY, ATH, GP, DE, SNQRK (Wang and Wang 1999)), k = 15, low-complexity k-mers kept
- protein20, k = 6 (`protein20_k6_lcFalse`): alphabet protein20 (the 20 amino acids, not reduced), k = 6, low-complexity k-mers kept
- protein20, k = 7 (`protein20_k7_lcFalse`): alphabet protein20 (the 20 amino acids, not reduced), k = 7, low-complexity k-mers kept
- protein20, k = 8 (`protein20_k8_lcFalse`): alphabet protein20 (the 20 amino acids, not reduced), k = 8, low-complexity k-mers kept

## Version, container, command

- Version: kmerseek 0.3.1 (printed in the search log)
- Container: `docker.io/olgabot/kmerseek:2026-08-24-reduced-alphabets`

```bash
kmerseek search \
    --alphabet hp_thomas_dill2 \
    --ksize    24 \
    --query    UP000005640_9606.fasta \
    --target   zebrafish.hp_thomas_dill2.k24.lcfalse.kmerseek.rocksdb \
    --threshold         0.0 \
    --min-shared-kmers  2 \
    --max-query-pvalue  0.05 \
    --min-region-score  1.3
```

Shown for setting hp_thomas_dill2 k24. Other settings change --alphabet and --ksize only.

Tool id `kmerseek`. Numbers from seanome/2024-kmerseek-analysis commit `4a4b1171`.

