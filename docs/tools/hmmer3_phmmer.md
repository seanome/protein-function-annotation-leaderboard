# phmmer

Tool id `hmmer3_phmmer`. Numbers from seanome/2024-kmerseek-analysis commit `6af1255b`.

## What it compares

Amino-acid letters, one query protein against each protein of the relative, scored with a substitution matrix the way BLAST does but with HMMER's model of insertions and deletions.

## What it needs

- sequence only: yes
- AlphaFold model: no
- GPU: no

## Reference searched

The relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as a FASTA file.

## How to read its score

Domain bit score: higher ranks first. The E-value cutoff was 10.

## Strong and weak strata

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is `default`: Fmax 0.333 over 2,435 true instances. Below, each stratum is compared with the best other tool and setting in the same stratum. Only strata with at least 30 true instances are compared.

Where it compares best with the other tools (its Fmax, then the best other):

- identity 20–30% to the closest same-family domain: Fmax 0.359 (n = 93); best other: MMseqs2, three rounds `s7`, 0.337.
- AlphaFold confidence (pLDDT) 50–70: Fmax 0.325 (n = 748); best other: MMseqs2 `s7`, 0.310.

Where it compares worst:

- regions 16–30 aa long: Fmax 0.042 (n = 271); best other: kmerseek `wwmj5_k11_lcFalse`, 0.330.
- regions 31–60 aa long: Fmax 0.095 (n = 420); best other: kmerseek `gbmr4_k21_lcFalse`, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting `default`: best feature type DOMAIN (domain), Fmax 0.304 (n = 1,107); worst MOTIF (short motif), Fmax 0.000 (n = 190).

## Known failure modes {#failure-modes}

- Failure modes specific to phmmer: not measured in this run.

## Cost per 1,000 queries

Setting `default`, searching zebrafish: 0.381 CPU-hours and 4.9 GB peak memory, measured on Sherlock (Stanford Sherlock, SLURM partition hns, CPU nodes, no GPU); about $0.017 at $0.044625 per CPU-hour.
Across the 15 relatives searched with this setting: 0.0189 to 0.381 CPU-hours per 1,000 queries; peak memory up to 4.9 GB.

CPU-hours count the search only, not building the relative's database.

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

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.252, 60–100% 0.308. If no, 20–30% 0.359, 30–40% 0.199. (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `default`; strata with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.250, 121–250 aa 0.329, 251 aa or more 0.477. If no, 16–30 aa 0.042, 31–60 aa 0.095. (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `default`.)
3. **Can you spend about 0.381 CPU-hours per 1,000 proteins for each relative you search?** That is the measured search cost against zebrafish; the range over relatives is under "Cost per 1,000 queries" above.
