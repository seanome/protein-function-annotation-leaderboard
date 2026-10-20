# ProstT5 + Foldseek

Tool id `prostt5`. Numbers from seanome/2024-kmerseek-analysis commit `7eca428b`.

## What it compares

Shape letters (3Di) predicted from the amino-acid sequence by the ProstT5 protein language model, then searched with Foldseek, so no AlphaFold model is needed.

## What it needs

- sequence only: yes
- AlphaFold model: no
- GPU: yes

## Reference searched

3Di letters predicted by ProstT5 for the relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as a Foldseek database.

## How to read its score

Bit score: higher ranks first. The E-value cutoff was 10.

## Strong and weak strata

Under the Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, its best setting is `3di_from_seq`: Fmax 0.223 over 2,435 true instances. Below, each stratum compares this setting with the best setting of any other tool in that same stratum, so the other tool's setting can change from line to line. Only strata with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- no same-family domain in the relative: Fmax 0.020 (n = 321); best other: HH-suite (hhsearch) `single_seq`, 0.032.
- identity 30–40% to the closest same-family domain: Fmax 0.144 (n = 165); best other: MMseqs2, three rounds `s7`, 0.209.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.000 (n = 271); best other: kmerseek `wwmj5_k11_lcFalse`, 0.330.
- regions 31–60 aa long: Fmax 0.032 (n = 420); best other: kmerseek `gbmr4_k21_lcFalse`, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting `3di_from_seq`: best feature type ACT_SITE (active site residue), Fmax 0.278 (n = 176); worst TRANSMEM (membrane-spanning segment), Fmax 0.000 (n = 1,124).

## Known failure modes {#failure-modes}

- Proteins longer than 6,000 residues are skipped before the 3Di letters are predicted, and are listed in a skipped-proteins file.
- Predicting the 3Di letters ran on one GPU (foldseek createdb --gpu 1). For the 998 queries that step is counted in the GPU column and the cost; for each relative's proteome it is a one-off database build and is not counted.

## Cost per 1,000 queries

Setting `3di_from_seq`, searching one relative (mean over the 15 relatives it searched): 0.109 core-hours reserved per 1,000 queries (16 cores for the task's wall time), measured on Sherlock. Peak memory 1.5 GB. That is about $0.00486 at $0.044625 per CPU-hour.
Building its databases took 5.31 core-hours over 3 databases, plus 0.664 GPU-hours (1 GPU via --gres=gpu:1 (model not recorded)). A database is built once and reused, so this is not in the per-query cost.

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Licence

- Licence: GPL-3.0 (Foldseek); MIT (ProstT5 model weights) ([source](https://huggingface.co/Rostlab/ProstT5))
- Commercial use allowed: Yes. Foldseek's GPL-3.0 applies to changed versions you distribute.

## Settings in the table

- `3di_from_seq`: 3Di letters predicted from sequence by ProstT5

## Version, container, command

- Version: Foldseek 10.941cd33 with ProstT5 weights
- Container: `quay.io/biocontainers/foldseek:10.941cd33--h5021889_1 (sha256:1156a052f31b2afb85257c02e83a962f559c9752273fe1064ab735f90ac29d1a)`

```bash
foldseek createdb filtered.fasta zebrafish_prostt5/db \
    --prostt5-model weights --gpu 1 --threads 8

foldseek search \
    human-38fdd87d_prostt5/db zebrafish_prostt5/db result tmp \
    --threads 16 \
    -e 10.0 \
    --max-seqs 1000
```

The first command (one per proteome, on a GPU) predicts the 3Di letters. The second is the search the cost columns measure.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.139 (n = 563), 60–100% 0.209 (n = 1,289). If no, 20–30% 0.247 (n = 93), 30–40% 0.144 (n = 165). (Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, setting `3di_from_seq`; strata with at least 30 true instances.)
2. **Do you have a GPU to predict the 3Di letters for your proteins?** The run used one GPU for that step; its cost is not counted here.
3. **Can you spend about 0.109 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" above.
