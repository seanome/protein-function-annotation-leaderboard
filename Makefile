# Source of the tables: a commit of seanome/2024-kmerseek-analysis that holds
# results/leaderboard/, written there by `make leaderboard-tables`.
SOURCE_REPO ?= seanome/2024-kmerseek-analysis
SOURCE_COMMIT ?= 4a4b1171
PILEUP_COMMIT ?= 39d81ec9
# No-search baseline rows and their cost: the Sherlock run (job 48547674) on
# olgabot/no-search-baselines, PR 158 of the source repo.
NOSEARCH_COMMIT ?= 5b644dc8
NOSEARCH_PATH ?= results/leaderboard
# Best-hit transfer rows (DIAMOND, MMseqs2 against Swiss-Prot minus Mammalia): the Sherlock
# run of head job 48547660, PR 160 of the source repo. Rows, cost rows and cards were
# committed at 7ede75df; later commits add the median IoU to call_rules.parquet and the two
# figures, and leave those three files byte for byte as they were. f0802ebe redraws the
# figures with TimeTree 5 divergence times (PR 148 at 7eca428b). sha256 of each file is
# pinned in data/best_hit_pins.json.
BESTHIT_COMMIT ?= f0802ebe
PY ?= .venv/bin/python

.PHONY: tables site check serve

tables:
	$(PY) scripts/copy_tables.py --repo $(SOURCE_REPO) --commit $(SOURCE_COMMIT) --pileup-commit $(PILEUP_COMMIT) --nosearch-commit $(NOSEARCH_COMMIT) --nosearch-path $(NOSEARCH_PATH) --besthit-commit $(BESTHIT_COMMIT) --out data

site:
	$(PY) scripts/build_site.py --data data --out docs

check:
	$(PY) scripts/claims.py --data data --check
	$(PY) tests/check_claims.py
	$(PY) tests/check_build.py --data data
	$(PY) tests/check_site.py --site docs
	$(PY) tests/check_pileup.py --site docs --data data
	$(PY) tests/check_explainers.py --site docs
	$(PY) tests/check_figures.py --site docs --data data
	$(PY) tests/check_prose_all.py --site docs

serve:
	$(PY) -m http.server -d docs 8765
