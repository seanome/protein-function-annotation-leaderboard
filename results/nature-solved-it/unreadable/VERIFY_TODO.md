# Second-agent check, still to do

Written 2026-10-20 for an agent that has not seen this analysis. Report a table of what held
and what failed; do not fix anything silently.

1. **Shares from the raw tables.** On Sherlock (read-only; compute in a job, not on a login
   node), for each of the 42 folders `$SCRATCH/nature-solved-it/unreadable/out/<upid>/`,
   recompute the unlabelled share from `mmseqs_vs_sprot.tsv.gz`, the proteome FASTA and the
   UniProt table in `$SCRATCH/nature-solved-it/unreadable/inputs/`, with the leaderboard's
   definition: neither an InterPro match (any non-empty InterPro or member-database column)
   nor an MMseqs2 hit with `fident > 0.30` (a hit to its own entry counts), over the
   accessions in the FASTA. Compare with `n_unlabelled` and `unlabelled_pct` in
   `run_2026-10-20/out/<upid>/summary.json` and `verdict.json`.
2. **Medians.** Recompute the group medians, all proteomes and BUSCO of at least 90%
   (`busco_completeness_pct`, kingdom dataset), from `verdict.json`, leaving out the
   "leaderboard reference" row. Expected: model 2.57 / 2.57, atlas 2.02 / 1.66, controls
   6.63 / 5.505; n with BUSCO of at least 90%: 5, 8, 4.
3. **The rule.** Apply the three sentences in `PLAN.md` to those medians and say which holds.
4. **The control draw.** Rerun `scripts/nature_solved_it/unreadable/draw_controls.py` with
   seed 20261019 on `atlas_proteomes.csv`, `model_organisms.csv`,
   `atlas_named_species.txt`, the committed `reference_proteomes_README.txt.gz` and a copy of
   the committed `timetree_cache.json` (out-dir a scratch folder), and check that the 17
   controls in `controls.csv` come out again.
5. **BUSCO.** For a few proteomes, read the C: value in
   `run_2026-10-20/out/<upid>/busco_<lineage>/short_summary*.txt` and compare with
   `verdict.json`.
6. **CPU-hours.** Recompute the total (38.18 CPU-hours) from `run_2026-10-20/sacct.txt`.
