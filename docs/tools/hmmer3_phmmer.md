# phmmer

Amino-acid letters, one query protein against each protein of the relative, scored with a substitution matrix the way BLAST does but with HMMER's model of insertions and deletions.

## Will it work for you?

1. **Is the closest same-family domain in the relative you will search at least 40% identical, for most of your proteins?** If yes, Fmax here is 40–60% 0.252 (n = 563), 60–100% 0.308 (n = 1,289). If no, 20–30% 0.359 (n = 93), 30–40% 0.199 (n = 165). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting default options; subsets with at least 30 true instances.)
2. **Are the features you care about longer than 60 residues?** If yes, Fmax here is 61–120 aa 0.250 (n = 963), 121–250 aa 0.329 (n = 504), 251 aa or more 0.477 (n = 277). If no, 16–30 aa 0.042 (n = 271), 31–60 aa 0.095 (n = 420). (Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, setting default options.)
3. **Can you spend about 0.417 core-hours per 1,000 proteins for each relative you search?** That is the measured mean per relative; see "Cost per 1,000 queries" below.

## Best Fmax and cost

Under the Pfam truth set, raw calls, zebrafish (430 Mya) as the relative searched, its best setting is default options (`default`): Fmax 0.333 over 2,435 true instances. Searching costs about $0.0186 per 1,000 queries (0.417 core-hours), for each relative searched.

## Known failure modes

- HMMER rounds E-values to two significant figures. Counting E ≤ cutoff on those rounded values, instead of E < cutoff, made the dark-set counts 2 to 222 lower. The dark set is the proteins that phmmer, jackhmmer and MMseqs2 all fail to place. The number: 2 to 222. Status: open. What to do: Say whether each cutoff is strict (E < 1) or not (E ≤ 1), and use the same rule in every count that is compared. Source: [seanome/2024-kmerseek-analysis PR 57, comment 6151921382](https://github.com/seanome/2024-kmerseek-analysis/pull/57#issuecomment-6151921382); [seanome/2024-kmerseek-analysis nextflow-runs/invertebrate-dark-set/main.nf at ec91da53](https://github.com/seanome/2024-kmerseek-analysis/blob/ec91da53b26d88474679989205045de45e38a162/nextflow-runs/invertebrate-dark-set/main.nf).
- A call is the whole alignment, not the target's domain placed on the query, so a long alignment over repeats never lands on one repeat. On a 23-aa LRR_1 repeat, every MMseqs2 and phmmer call was 47 aa or longer, so none landed. The number: 47 aa or longer, on a 23-aa repeat. Status: open. What to do: Place the target's domain on the query through the alignment before scoring. Notebook 280 switched to this rule; the leaderboard's own calls are still whole alignments. Source: [seanome/2024-kmerseek-analysis PR 121, comment 6154499159](https://github.com/seanome/2024-kmerseek-analysis/pull/121#issuecomment-6154499159); [seanome/2024-kmerseek-analysis PR 121, comment 6154499159](https://github.com/seanome/2024-kmerseek-analysis/pull/121#issuecomment-6154499159); [seanome/2024-kmerseek-analysis PR 121, comment 6169102356](https://github.com/seanome/2024-kmerseek-analysis/pull/121#issuecomment-6169102356).
- Aligner calls run past the domain, because a call spans the whole aligned stretch. Over correct calls, they start a median 5 to 21 residues early. Over all calls, 45% to 60% start more than 100 residues early. The number: median 5 to 21 residues early; 45% to 60% more than 100 residues early. Status: open. What to do: Read a call's ends as the alignment's ends, not the domain's. To score where a domain sits, place the target's domain through the alignment. Source: [seanome/2024-kmerseek-analysis PR 146](https://github.com/seanome/2024-kmerseek-analysis/pull/146); [PR 3, comment 6193138553](https://github.com/seanome/protein-function-annotation-leaderboard/pull/3#issuecomment-6193138553).

## Where it is strong and weak

Each line compares this setting with the best setting of any other tool on the same subset of the proteins, so the other tool's setting can change from line to line. Only subsets with at least 30 true instances are compared.

Where it is closest to, or furthest ahead of, the best other tool:

- identity 20–30% to the closest same-family domain: Fmax 0.359 (n = 93); best other: MMseqs2, three rounds, sensitivity 7, 3 rounds, 0.337.
- AlphaFold confidence (pLDDT) 50–70: Fmax 0.325 (n = 748); best other: MMseqs2, sensitivity 7, 0.310.

Where it is furthest behind the best other tool:

- regions 16–30 aa long: Fmax 0.042 (n = 271); best other: kmerseek, wwmj5, k = 11, 0.330.
- regions 31–60 aa long: Fmax 0.095 (n = 420); best other: kmerseek, gbmr4, k = 21, 0.232.

Swiss-Prot features, raw calls, zebrafish, setting default options: best feature type domain, Fmax 0.304 (n = 1,107); worst short motif, Fmax 0.000 (n = 190).


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

## Caveats of how this run was set up {#caveats}

- None recorded for this run.

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

Tool id `hmmer3_phmmer`. Numbers from seanome/2024-kmerseek-analysis commit `4a4b1171`.

