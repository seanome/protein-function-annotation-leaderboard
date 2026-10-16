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

## The organism pages (home page)

`index.html` asks one thing: "I work on". The species' page gives four facts, three tools in
the order to run them, and a command to copy. Pre-built pages (`organism/<species>.html`, each
with `.html.json`) cover every species on the Earth BioGenome Project target list that has a
UniProt reference proteome, plus human, mouse, yeast, E. coli and Arabidopsis. Any other
species is built in the browser by `organism/live.html#taxid=<NCBI taxon id>`, from the same
code (`docs/assets/organism.js`, which the build also runs in Node).

| fact | source |
|---|---|
| nearest species with at least 100 reviewed entries, and the million years to it | UniProt 2026_03 counts over NCBI Taxonomy; TimeTree 5 (pairwise for pre-built pages, the shared ancestor's age in the browser, because TimeTree's API cannot be called from a web page) |
| reviewed entries in the species' order and class | UniProt 2026_03 |
| share of the proteome with an AlphaFold model; mean pLDDT (estimate, from model counts in 5-point bands) | UniProt cross-references; AlphaFold DB search API |
| on the Earth BioGenome Project target list | GoaT, attribute `long_list` |

The lookups are refreshed with `make organism-lookups` (downloads go to `$HOME/data/organism-page`)
and land in `data/organism/`, with `sources.json` recording each address, release and sha256.
The leaderboard moved to `leaderboard.html`; an old link with filters after `#` is sent there.

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
