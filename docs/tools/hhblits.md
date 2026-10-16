# HH-suite (hhsearch)

Tool id `hhblits`. Numbers from seanome/2024-kmerseek-analysis commit `4a4b1171`.

## What it compares

Two profiles against each other. Here each profile was built from a single protein, not from a family alignment, because no UniRef30 database was given.

## What it needs

- sequence only: yes
- 3D structure prediction: no
- GPU: no

## Reference searched

Single-protein profiles of the relative's proteome (Quest for Orthologs reference proteomes, release 2020_04), as an HH-suite database.

## How to read its score

Bit score from hhsearch: higher ranks first.

## Strong and weak strata

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is `single_seq`: Fmax 0.313 over 2,435 true instances. Below, each stratum compares this setting with the best setting of any other tool in that same stratum, so the other tool's setting can change from line to line. Only strata with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- no same-family domain in the relative: Fmax 0.032 (n = 321); best other: Foldseek `3di_aa`, 0.027.
- identity 40–60% to the closest same-family domain: Fmax 0.248 (n = 563); best other: phmmer `default`, 0.252.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.066 (n = 271); best other: kmerseek `wwmj5_k11_lcFalse`, 0.330.
- regions 31–60 aa long: Fmax 0.085 (n = 420); best other: kmerseek `gbmr4_k21_lcFalse`, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting `single_seq`: best feature type DOMAIN (domain), Fmax 0.294 (n = 1,107); worst MOTIF (short motif), Fmax 0.003 (n = 190).

## Known failure modes {#failure-modes}

- The run names this tool hhblits (its tool id), but it ran hhsearch on profiles of single proteins. HH-suite's usual use, with family alignments from UniRef30, was not run, so these numbers say nothing about it.

## Cost per 1,000 queries

Setting `single_seq`, searching one relative (mean over the 15 relatives it searched): 24.8 core-hours reserved per 1,000 queries (16 cores for the task's wall time), measured on Sherlock. Peak memory 23.4 GB. That is about $1.11 at $0.044625 per CPU-hour.
Building its databases took 30.9 core-hours over 7 databases. A database is built once and reused, so this is not in the per-query cost.

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Hardware

Hardware of the search (setting `single_seq`, every task): Sherlock partition hns, from the trace. Host name: not recorded in the trace, so the node is the one SLURM's job record (sacct) gives. 15 tasks on 13 nodes. CPU: AMD 7543 (4 tasks); AMD 8224P (4 tasks); AMD 9754 (3 tasks); AMD 7742 (2 tasks); not recorded, as scontrol no longer finds the node (2 tasks). GPU: none allocated.

Hardware of the database build: Sherlock partition hns, from the trace. Host name: not recorded in the trace, so the node is the one SLURM's job record (sacct) gives. 7 tasks on 7 nodes. CPU: AMD 8224P (4 tasks); AMD 7543 (1 task); AMD 9754 (1 task); not recorded, as scontrol no longer finds the node (1 task). GPU: none allocated.

Sources, read 2026-10-16: the partition is the `queue` column of the run's Nextflow traces. The traces have no `hostname` column, so the node is the one SLURM's job record (`sacct`) gives for the trace's job id. CPU and GPU are that node's features in SLURM (`scontrol show node`: CPU_MNF, CPU_SKU, GPU_SKU), named as on Sherlock's [node features page](https://www.sherlock.stanford.edu/docs/advanced-topics/node-features/). Partition hns is owned by a research group, and Sherlock's [tech specs page](https://www.sherlock.stanford.edu/docs/tech/) lists no CPU models for such partitions. Script: `scripts/read_hardware.py`.

## Licence

- Licence: GPL-3.0 ([source](https://github.com/soedinglab/hh-suite))
- Commercial use allowed: Yes. If you distribute a changed version, you must release its source under GPL-3.0 too.

## Settings in the table

- `single_seq`: profiles built from single proteins, no family database

## Version, container, command

- Version: HH-suite 3.3.0
- Container: `quay.io/biocontainers/hhsuite:3.3.0--h503566f_15 (sha256:4bf9bb5229de18f522a94f4443c19fdcbb0f0cb0e6ea92f5390aa170bcb0a24f)`

```bash
ffindex_apply \
    human-38fdd87d_hhdb/db_a3m.ffdata human-38fdd87d_hhdb/db_a3m.ffindex \
    -d results.ffdata -i results.ffindex \
    -- hhsearch -i stdin -d zebrafish_hhdb/db -blasttab /dev/stdout -cpu 16 -v 0
```

One hhsearch per query protein, run through ffindex_apply.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.248 (n = 563), 60–100% 0.283 (n = 1,289). If no, 20–30% 0.321 (n = 93), 30–40% 0.190 (n = 165). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `single_seq`; strata with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.230 (n = 963), 121–250 aa 0.298 (n = 504), 251 aa or more 0.417 (n = 277). If no, 16–30 aa 0.066 (n = 271), 31–60 aa 0.085 (n = 420). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `single_seq`.)
3. **Can you spend about 24.8 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" above.
