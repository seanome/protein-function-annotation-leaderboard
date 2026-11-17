# Protein annotation leaderboard

![claims](https://seanome.github.io/protein-function-annotation-leaderboard/badges/claims-kmerseek.svg)

Which tool finds the known domains and features of a protein, and what does it cost to run?
A static site, served by GitHub Pages from `docs/`:
https://seanome.github.io/protein-function-annotation-leaderboard/

## Where the numbers come from

`data/` holds four files copied from one commit of
[seanome/2024-kmerseek-analysis](https://github.com/seanome/2024-kmerseek-analysis), folder
`results/leaderboard/`, where `make leaderboard-tables` writes them:

| file | holds |
|---|---|
| `results.parquet` | Fmax, AUPRC, median IoU and counts per truth set, dedup, tool, setting, relative and stratum |
| `compute.parquet` | CPU-hours and peak memory of each search, from the run's Nextflow task records |
| `coverage.md` | tool by stratum: measured or not run |
| `provenance.json` | sha256 of every input and output |

`data/SOURCE.json` records the commit, and every page prints it. The site never fills a gap:
where the table has no row, the page shows "—" and "not run".

## Build and check

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

```bash
make tables SOURCE_COMMIT=<commit of 2024-kmerseek-analysis>
```

```bash
make site
```

```bash
make check
```

`make check` renders every filter combination in headless Chrome (the installed Google
Chrome, through Playwright) and compares every number and every "—" with the parquet files,
then runs the prose checker on every page.

## Layout

- `scripts/copy_tables.py`: copies `results/leaderboard/` from a commit, checks the sha256 sums
- `scripts/build_site.py`: writes `docs/` from `data/`
- `scripts/tools.py`: what each tool compares, needs, and was run with (read from the run itself)
- `docs/assets/`: stylesheet, page script, Fantasque Sans Mono 1.8.0 (SIL Open Font License)
- `tests/check_site.py`, `tests/check_prose_all.py`: the checks above

## Claims

`claims/<tool>.yaml` lists a tool's claims, each with the result that would retire it
(schema: `claims/schema.json`). `scripts/claims.py` sets each status from its test: table
tests are recomputed from `data/results.parquet`, record tests must quote the experiment
inventory (`data/experiments.json`, from seanome/2024-kmerseek-analysis PR 97) exactly. The
`claims` GitHub Action runs `python scripts/claims.py --check` on every push; the build writes
`docs/badges/claims-<tool>.svg`.

## Contributions

[CONTRIBUTING.md](CONTRIBUTING.md) has both ways in. In short:

- `data/query_sets/<name>/`: the fixed query sets a tool author runs on, with their answers
  (`scripts/make_query_set.py` copied them from the region benchmark's midi-plus run).
- `submissions/<name>/`: one folder per tool submission. `.github/workflows/score-submission.yml`
  scores a pull request that adds one, with `scripts/score_submission.py`, and posts the rows.
- `scorer/`: byte-for-byte copies of the leaderboard's scorer; `tests/check_scorer.py` checks
  their sha256 against `scorer/SOURCE.json`.
- `/contribute/truth/` (`scripts/contribute_pages.py`, `docs/assets/truth-form.js`): the
  wet-lab form. It opens a filled-in GitHub issue; `scripts/truth_contributed.py` turns accepted
  issues into `data/truth_contributed.tsv` ([docs/curation.md](docs/curation.md)).
- `docs/assets/counts.js`: the contributed-rows count for the organism page, and the Copy
  command counter for the organism and Why pages.

### Usage counts

The Copy command counter needs somewhere to store a number, and a GitHub Pages site cannot
store anything: a page can only write to GitHub with a token, and a token in a public page is
public. So the counter is off. `docs/data/counter.json` has `"endpoint": null`, the button sends
nothing, and the Why page says copies are not counted. Turning it on needs a small service
outside GitHub that keeps one number per (organism, tool set) and nothing about the visitor.
Setting `endpoint` to that service's address turns the counter on; `counts.js` then sends one
request per click, with no cookies, no referrer and no identifier, and reads the totals back
from `<endpoint>/counts`.
