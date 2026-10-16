# MMseqs2

Tool id `mmseqs2_seqseq`. Numbers from seanome/2024-kmerseek-analysis commit `4a4b1171`.

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

## Strong and weak strata

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is `s7`: Fmax 0.327 over 2,435 true instances. Below, each stratum compares this setting with the best setting of any other tool in that same stratum, so the other tool's setting can change from line to line. Only strata with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- identity 40–60% to the closest same-family domain: Fmax 0.249 (n = 563); best other: phmmer `default`, 0.252.
- AlphaFold confidence (pLDDT) 90–100: Fmax 0.480 (n = 186); best other: MMseqs2, three rounds `s7`, 0.484.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.027 (n = 271); best other: kmerseek `wwmj5_k11_lcFalse`, 0.330.
- regions 31–60 aa long: Fmax 0.089 (n = 420); best other: kmerseek `gbmr4_k21_lcFalse`, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting `s7`: best feature type ACT_SITE (active site residue), Fmax 0.312 (n = 176); worst MOTIF (short motif), Fmax 0.000 (n = 190).

## Known failure modes {#failure-modes}

- At most 1,000 target proteins are kept per query (--max-seqs 1000), so a query from a very large family can lose hits. How often: not measured.

## Cost per 1,000 queries

Setting `s7`, searching one relative (mean over the 15 relatives it searched): 0.0252 core-hours reserved per 1,000 queries (16 cores for the task's wall time), measured on Sherlock. Peak memory 1.1 GB. That is about $0.00113 at $0.044625 per CPU-hour.
Building its databases took 0.00889 core-hours over 7 databases. A database is built once and reused, so this is not in the per-query cost.

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Hardware

Hardware of the search (setting `s7`, every task): Sherlock partition hns, from the trace. Host name: not recorded in the trace, so the node is the one SLURM's job record (sacct) gives. 15 tasks on 9 nodes. CPU: AMD 7543 (5 tasks); AMD 9754 (5 tasks); AMD 8224P (3 tasks); not recorded, as scontrol no longer finds the node (2 tasks). GPU: none allocated.

Hardware of the database build: Sherlock partition hns, from the trace. Host name: not recorded in the trace, so the node is the one SLURM's job record (sacct) gives. 7 tasks on 5 nodes. CPU: AMD 8224P (3 tasks); AMD 7543 (1 task); AMD 7742 (1 task); AMD 9754 (1 task); not recorded, as scontrol no longer finds the node (1 task). GPU: none allocated.

Sources, read 2026-10-16: the partition is the `queue` column of the run's Nextflow traces. The traces have no `hostname` column, so the node is the one SLURM's job record (`sacct`) gives for the trace's job id. CPU and GPU are that node's features in SLURM (`scontrol show node`: CPU_MNF, CPU_SKU, GPU_SKU), named as on Sherlock's [node features page](https://www.sherlock.stanford.edu/docs/advanced-topics/node-features/). Partition hns is owned by a research group, and Sherlock's [tech specs page](https://www.sherlock.stanford.edu/docs/tech/) lists no CPU models for such partitions. Script: `scripts/read_hardware.py`.

## Licence

- Licence: MIT ([source](https://github.com/soedinglab/MMseqs2))
- Commercial use allowed: Yes.

## Settings in the table

- `s7`: sensitivity 7 (-s 7), one round

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

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.249 (n = 563), 60–100% 0.297 (n = 1,289). If no, 20–30% 0.325 (n = 93), 30–40% 0.188 (n = 165). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `s7`; strata with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.243 (n = 963), 121–250 aa 0.305 (n = 504), 251 aa or more 0.462 (n = 277). If no, 16–30 aa 0.027 (n = 271), 31–60 aa 0.089 (n = 420). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting `s7`.)
3. **Can you spend about 0.0252 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" above.
