# Foldseek

Predicted 3-D shapes.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.164 (n = 563), 60–100% 0.241 (n = 1,289). If no, 20–30% 0.268 (n = 93), 30–40% 0.167 (n = 165). (Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, setting shape letters + amino acids; subsets with at least 30 true instances.)
2. **Do you have a 3D structure prediction for every protein you want to annotate?** If not, Foldseek cannot search the proteins without one.
3. **Can you spend about 0.129 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" below.

## Best Fmax and cost

Under the Pfam truth set, raw calls, zebrafish (435 Mya) as the relative searched, its best setting is shape letters + amino acids (`3di_aa`): Fmax 0.250 over 2,435 true instances. Searching costs about $0.00574 per 1,000 queries (0.129 core-hours), for each relative searched.

## Known failure modes

- In the ELM motif run, Foldseek was scored on orthologs that have no AlphaFold model. It cannot search those, so its share of motifs found is low for a reason that is not the search. The number: The source does not count the orthologs with no AlphaFold model.. Status: open. What to do: Score a structure tool only on proteins that have a structure, and print how many proteins were left out. Source: [seanome/2024-kmerseek-analysis results/leaderboard/inventory.md at 6df20c1f](https://github.com/seanome/2024-kmerseek-analysis/blob/6df20c1f290fc123781df8d58097bf7b836535ac/results/leaderboard/inventory.md).
- Aligner calls run past the domain, because a call spans the whole aligned stretch. Over correct calls, they start a median 5 to 21 residues early. Over all calls, 45% to 60% start more than 100 residues early. The number: median 5 to 21 residues early; 45% to 60% more than 100 residues early. Status: open. What to do: Read a call's ends as the alignment's ends, not the domain's. To score where a domain sits, place the target's domain through the alignment. Source: [seanome/2024-kmerseek-analysis PR 146](https://github.com/seanome/2024-kmerseek-analysis/pull/146); [PR 3, comment 6193138553](https://github.com/seanome/protein-function-annotation-leaderboard/pull/3#issuecomment-6193138553).

## Where it is strong and weak

Each line compares this setting with the best setting of any other tool on the same subset of the proteins, so the other tool's setting can change from line to line. Only subsets with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- no same-family domain in the relative: Fmax 0.027 (n = 321); best other: HH-suite (hhsearch), single-protein profiles, 0.032.
- identity 30–40% to the closest same-family domain: Fmax 0.167 (n = 165); best other: MMseqs2, three rounds, sensitivity 7, 3 rounds, 0.209.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.004 (n = 271); best other: kmerseek, wwmj5, k = 11, 0.330.
- regions 31–60 aa long: Fmax 0.034 (n = 420); best other: kmerseek, gbmr4, k = 21, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting shape letters + amino acids: best feature type active site, Fmax 0.285 (n = 176); worst membrane-spanning, Fmax 0.000 (n = 1,124).


## What it compares

Predicted 3-D shapes. Each residue of an AlphaFold model gets one of 20 shape letters (3Di) from its local fold. The shape letters are aligned together with the amino acids.

## What it needs

- sequence only: no
- 3D structure prediction: yes
- GPU: no

## Reference searched

AlphaFold models of the relative's proteome (Quest for Orthologs reference proteomes, release 2020_04 accessions), as a Foldseek database.

## How to read its score

Bit score: higher ranks first. The E-value cutoff was 10.

## Caveats of how this run was set up {#caveats}

- A protein with no AlphaFold model cannot be searched. Five relatives have no Foldseek rows at all (chimpanzee, cow, dog, opossum, western clawed frog): not run.
- AlphaFold models of long proteins come in overlapping pieces (F1, F2, ...). The task script moves piece n's coordinates 200 x (n - 1) residues along the protein; if a piece started elsewhere, its calls would land in the wrong place. Checked against the real piece starts: not measured.

## Cost per 1,000 queries

Setting `3di_aa`, searching one relative (mean over the 10 relatives it searched): 0.129 core-hours reserved per 1,000 queries (16 cores for the task's wall time), measured on Sherlock. Peak memory 1.5 GB. That is about $0.00574 at $0.044625 per CPU-hour.
Building its databases took 0.267 core-hours over 2 databases. A database is built once and reused, so this is not in the per-query cost.
This covers the 2 database builds the task records hold; the search used 10 relatives.

The per-query cost counts the search only. The run recorded cost per tool and setting, summed over all relatives, so the cost for one relative is the mean.

## Licence

- Licence: GPL-3.0 ([source](https://github.com/steineggerlab/foldseek))
- Commercial use allowed: Yes. If you distribute a changed version, you must release its source under GPL-3.0 too.

## Settings in the table

- shape letters + amino acids (`3di_aa`): 3Di shape letters and amino acids together

## Version, container, command

- Version: Foldseek 10.941cd33
- Container: `quay.io/biocontainers/foldseek:10.941cd33--h5021889_1 (sha256:1156a052f31b2afb85257c02e83a962f559c9752273fe1064ab735f90ac29d1a)`

```bash
foldseek search \
    human-38fdd87d_fsdb/db zebrafish_fsdb/db result foldseek_tmp \
    --threads 16 \
    -e 10.0 \
    --max-seqs 1000
```

Followed by foldseek convertalis to write the table of hits.

Tool id `foldseek`. Numbers from seanome/2024-kmerseek-analysis commit `7eca428b`.

