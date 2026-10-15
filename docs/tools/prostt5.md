# ProstT5 + Foldseek

Tool id `prostt5`. Numbers from seanome/2024-kmerseek-analysis commit `6af1255b`.

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

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is `3di_from_seq`: Fmax 0.223 over 2,435 true instances. Below, each stratum is compared with the best other tool and setting in the same stratum. Only strata with at least 30 true instances are compared.

Where it compares best with the other tools (its Fmax, then the best other):

- no same-family domain in the relative: Fmax 0.020 (n = 321); best other: HH-suite (hhsearch) `single_seq`, 0.032.
- identity 30–40% to the closest same-family domain: Fmax 0.144 (n = 165); best other: MMseqs2, three rounds `s7`, 0.209.

Where it compares worst:

- regions 16–30 aa long: Fmax 0.000 (n = 271); best other: kmerseek `wwmj5_k11_lcFalse`, 0.330.
- regions 31–60 aa long: Fmax 0.032 (n = 420); best other: kmerseek `gbmr4_k21_lcFalse`, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting `3di_from_seq`: best feature type ACT_SITE (active site residue), Fmax 0.278 (n = 176); worst TRANSMEM (membrane-spanning segment), Fmax 0.000 (n = 1,124).

## Known failure modes {#failure-modes}

- Proteins longer than 6,000 residues are skipped before the 3Di letters are predicted, and are listed in a skipped-proteins file.
- Predicting the 3Di letters ran on one GPU (foldseek createdb --gpu 1). That step is not counted in the cost columns, which count the search only.

## Cost per 1,000 queries

Setting `3di_from_seq`, searching zebrafish: 0.0861 CPU-hours and 1.5 GB peak memory, measured on Sherlock (Stanford Sherlock, SLURM partition hns, CPU nodes, no GPU); about $0.00384 at $0.044625 per CPU-hour.
Across the 15 relatives searched with this setting: 0.0541 to 0.123 CPU-hours per 1,000 queries; peak memory up to 1.5 GB.

CPU-hours count the search only, not building the relative's database.

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

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.139, 60–100% 0.209. If no, 20–30% 0.247, 30–40% 0.144. (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `3di_from_seq`; strata with at least 30 true instances.)
2. **Do you have a GPU to predict the 3Di letters for your proteins?** The run used one GPU for that step; its cost is not counted here.
3. **Can you spend about 0.0861 CPU-hours per 1,000 proteins for each relative you search?** That is the measured search cost against zebrafish; the range over relatives is under "Cost per 1,000 queries" above.
