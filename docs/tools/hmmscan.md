# hmmscan with Pfam

Each query protein against the curated profile models of every Pfam family, without searching any relative.

## Will it work for you?

1. **Is every family you expect to find already in Pfam?** If no, this tool cannot find it: it only reports families that have a Pfam model.
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.785 (n = 963), 121–250 aa 0.802 (n = 504), 251 aa or more 0.871 (n = 277). If no, 16–30 aa 0.580 (n = 271), 31–60 aa 0.573 (n = 420). (Pfam truth set, raw calls, no relative searched, setting Pfam models.)
3. **Can you spend about 0.389 core-hours per 1,000 proteins?** That is the measured cost of scanning them against Pfam's models.

## Best Fmax and cost

Under the Pfam truth set, raw calls, no relative (hmmscan searches the Pfam models), its best setting is Pfam models (`direct`): Fmax 0.886 over 2,435 true instances. Searching costs about $0.0173 per 1,000 queries (0.389 core-hours), for scanning Pfam's models.

## Known failure modes

None recorded yet.

## Where it is strong and weak

Each line compares this setting with the best setting of any other tool on the same subset of the proteins, so the other tool's setting can change from line to line. Only subsets with at least 30 true instances are compared.

No other tool was scored on the same subsets, so no comparison: not run.

Swiss-Prot feature types: not run.


## What it compares

Each query protein against the curated profile models of every Pfam family, without searching any relative.

## What it needs

- sequence only: yes
- 3D structure prediction: no
- GPU: no

## Reference searched

Pfam-A 38.2 family models (30,134 models).

## How to read its score

Domain bit score: higher ranks first. E-value cutoffs 10 (sequence and domain).

## Caveats of how this run was set up {#caveats}

- Pfam domains are defined by these same models, so scoring hmmscan against the Pfam truth set grades the models against their own output. Read its Pfam row as the best a search could do, not as a competitor.
- It was scored against the Pfam truth set only. Swiss-Prot features and Pfam-N: not run.

## Cost per 1,000 queries

Setting `direct`, scanning Pfam's models: 0.389 core-hours reserved per 1,000 queries (4 cores for the task's wall time), measured on Sherlock. Peak memory 0.5 GB. That is about $0.0173 at $0.044625 per CPU-hour.
Building its databases for setting `direct`: not measured (no task record holds a database build for this tool and setting).

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Licence

- Licence: BSD 3-clause (HMMER); CC0 1.0 (Pfam-A models, per the Pfam 38.2 release notes) ([source](https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE))
- Commercial use allowed: Yes.

## Settings in the table

- Pfam models (`direct`): scans Pfam's models; no relative searched

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

Tool id `hmmscan`. Numbers from seanome/2024-kmerseek-analysis commit `7eca428b`.

