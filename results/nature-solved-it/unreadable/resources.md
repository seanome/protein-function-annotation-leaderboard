# Compute for the solvers-unreadable run, predicted from job 47961379

## What was measured

The only measured history is Sherlock job 47961379 (2026-10-15, `scripts/why/sherlock_labels.sbatch`).
Slurm's record of it, read with `sacct` on 2026-10-19, is in `sacct_47961379.txt`.

| Task | What ran | Proteins | CPUs | Wall clock | CPU time (TotalCPU) | Peak memory (MaxRSS) | Memory asked |
|---|---|---|---|---|---|---|---|
| 0 | MMseqs2 easy-search against Swiss-Prot | 50,928 (human 20,652 + oyster 30,272 + 4) | 16 | 8 min 26 s | 2 h 01 min 26 s | 4.24 GiB | 32 GB |
| 1 | BUSCO metazoa_odb12, human | 20,652 | 16 | 1 min 54 s | 9 min 49 s | 2.49 GiB | 32 GB |
| 2 | BUSCO metazoa_odb12, oyster | 30,272 | 16 | 1 min 58 s | 11 min 05 s | 2.77 GiB | 32 GB |

From these, `scripts/nature_solved_it/unreadable/predict_resources.py` predicts:

- MMseqs2: 0.0397 CPU-hours and 0.00276 wall-clock hours per 1,000 proteins at 16 CPUs,
  scaled in proportion to the number of proteins.
- BUSCO: the straight line through the two measured runs. CPU-hours = 0.1181 + 2.2 × 10⁻⁶ ×
  proteins. Wall-clock hours = 0.0293 + 1.16 × 10⁻⁷ × proteins. Most of a BUSCO run is a
  fixed cost.
- Memory: 4.24 GiB was the largest peak. Every label task keeps the 32 GB that worked.

## Not measured

- BUSCO on `fungi_odb12`, `viridiplantae_odb12` and `eukaryota_odb12` has never run here. It is
  predicted with the `metazoa_odb12` line.
- MMseqs2 memory above 50,928 queries. The largest eukaryotic reference proteome in UniProt
  2026_03 has 246,042 proteins (*Araneus ventricosus*).
- MMseqs2's fixed cost per run. Job 47961379 was one run for both proteomes, so the time
  spent preparing Swiss-Prot cannot be split from the time per protein. One run per proteome
  repeats that cost, so small proteomes may take longer than predicted.
- Download time on Sherlock. The fetch step uses the slowest paged UniProt download of the
  leaderboard's `fetch.py`, run on the laptop on 2026-10-15: 57,874 rows in 5 min 46 s, or
  167 rows a second.

## Model organisms, predicted (`resources_model_organisms.tsv`)

| Species | Proteins | UniProt entries to download | BUSCO lineages | CPU-hours | Wall clock (h) | Download (h) |
|---|---|---|---|---|---|---|
| *Homo sapiens* | 20,652 | 147,520 | metazoa + eukaryota | 1.15 | 0.12 | 0.25 |
| *Mus musculus* | 21,860 | 54,855 | metazoa + eukaryota | 1.20 | 0.12 | 0.09 |
| *Danio rerio* | 26,567 | 69,241 | metazoa + eukaryota | 1.41 | 0.14 | 0.12 |
| *Drosophila melanogaster* | 13,814 | 21,947 | metazoa + eukaryota | 0.85 | 0.10 | 0.04 |
| *Caenorhabditis elegans* | 19,792 | 26,629 | metazoa + eukaryota | 1.11 | 0.12 | 0.04 |
| *Saccharomyces cerevisiae* | 6,066 | 6,067 | fungi + eukaryota | 0.50 | 0.08 | 0.01 |
| *Arabidopsis thaliana* | 27,496 | 39,272 | viridiplantae + eukaryota | 1.45 | 0.14 | 0.07 |

Proteins are the counts in the 2026_03 reference proteomes README. For human the README says
147,503, but the FASTA, replaced on 2026-09-15, holds 20,652 (counted by the leaderboard), so
20,652 is used.

## The whole run, before the atlas list exists

The atlas size is not known yet. Prompt 1 aims at 60 diseases. As a ceiling, take 80 atlas
species with a reference proteome, 80 controls, the 7 model organisms and the oyster: 168
proteomes. Eukaryotic reference proteomes in UniProt 2026_03 have a median of 13,888 proteins
and a 90th percentile of 31,785.

| If every proteome had | CPU-hours per proteome | CPU-hours, 168 proteomes | Wall clock, 20 tasks at once | Wall clock, 4 tasks at once |
|---|---|---|---|---|
| 13,888 proteins (median) | 0.85 | 143 | 0.8 h | 4.2 h |
| 31,785 proteins (90th percentile) | 1.64 | 275 | 1.3 h | 6.5 h |

The largest single proteome (246,042 proteins) is predicted at 11.1 CPU-hours and 0.8 h of
wall clock. The run fits in two compute days with a wide margin, even at 4 tasks at once. The
queue wait is not counted.

## Asks

`submit_all.sh` sets the label array's `--time` to 3 times the largest predicted task, at
least 1 hour, rounded up to the hour, from `predict_resources.py --time-ask label`. The fetch
array gets the same rule but at least 2 hours, because its speed on Sherlock is unmeasured. A
fetch task that runs out of time can be resubmitted; files already downloaded and checked are
not fetched again. Setup asks 4 CPUs, 16 GB and 3 hours, because since 2026-10-20 it builds the run's own BUSCO
image. The one earlier build on Sherlock took about 10 minutes (its lock file is dated
07:22 and the finished image 07:32 on 2026-10-15); nothing else about setup is measured.
Collect asks 30 minutes, unmeasured.
