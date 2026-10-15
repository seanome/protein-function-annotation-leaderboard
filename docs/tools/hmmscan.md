# hmmscan with Pfam

Tool id `hmmscan`. Numbers from seanome/2024-kmerseek-analysis commit `6af1255b`.

## What it compares

Each query protein against the curated profile models of every Pfam family, without searching any relative.

## What it needs

- sequence only: yes
- AlphaFold model: no
- GPU: no

## Reference searched

Pfam-A 38.2 family models (30,134 models).

## How to read its score

Domain bit score: higher ranks first. E-value cutoffs 10 (sequence and domain).

## Strong and weak strata

Under the Pfam truth set, raw calls, no relative (hmmscan searches the Pfam models), its best setting is `direct`: Fmax 0.886 over 2,435 true instances. Below, each stratum is compared with the best other tool and setting in the same stratum. Only strata with at least 30 true instances are compared.

No other tool was scored in the same strata, so no comparison: not run.

Swiss-Prot feature types: not run.

## Known failure modes {#failure-modes}

- Pfam domains are defined by these same models, so scoring hmmscan against the Pfam truth set grades the models against their own output. Read its Pfam row as the best a search could do, not as a competitor.
- It was scored against the Pfam truth set only. Swiss-Prot features and Pfam-N: not run.

## Cost per 1,000 queries

Setting `direct`, scanning Pfam's models: 0.198 CPU-hours and 0.5 GB peak memory, measured on Sherlock (Stanford Sherlock, SLURM partition hns, CPU nodes, no GPU); about $0.00883 at $0.044625 per CPU-hour.
Across the 1 relatives searched with this setting: 0.198 to 0.198 CPU-hours per 1,000 queries; peak memory up to 0.5 GB.

CPU-hours count the search only, not building the relative's database.

## Version, container, command

- Version: HMMER 3.4
- Container: `quay.io/biocontainers/hmmer:3.4--hb6cb901_4 (sha256:7a2b317b8d2fd3650b4924a8482cddeb940d4a0746c6a1501ff03ac1b7439e0c)`

```bash
hmmscan \
    --domtblout /dev/stdout \
    -o /dev/null \
    --noali \
    -E 10.0 --domE 10.0 \
    --cpu 4 \
    Pfam-A.hmm UP000005640_9606.fasta
```

Query: the same 998 human proteins.

## Will it work for you?

1. **Is every family you expect to find already in Pfam?** If no, this tool cannot find it: it only reports families that have a Pfam model.
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.785, 121–250 aa 0.802, 251 aa or more 0.871. If no, 16–30 aa 0.580, 31–60 aa 0.573. (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `direct`.)
3. **Can you spend about 0.198 CPU-hours per 1,000 proteins for each relative you search?** That is the measured search cost against none: Pfam models (hmmscan only); the range over relatives is under "Cost per 1,000 queries" above.
