# Source of the tables: a commit of seanome/2024-kmerseek-analysis that holds
# results/leaderboard/, written there by `make leaderboard-tables`.
SOURCE_REPO ?= seanome/2024-kmerseek-analysis
SOURCE_COMMIT ?= 7eca428b
PILEUP_COMMIT ?= 39d81ec9
# No-search baseline rows and their cost: the Sherlock run (job 48547674) on
# olgabot/no-search-baselines, PR 158 of the source repo.
NOSEARCH_COMMIT ?= 5b644dc8
NOSEARCH_PATH ?= results/leaderboard
PY ?= .venv/bin/python

.PHONY: tables site check serve

tables:
	$(PY) scripts/copy_tables.py --repo $(SOURCE_REPO) --commit $(SOURCE_COMMIT) --pileup-commit $(PILEUP_COMMIT) --nosearch-commit $(NOSEARCH_COMMIT) --nosearch-path $(NOSEARCH_PATH) --out data

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
	$(PY) tests/check_scope40.py --site docs --data data
	$(PY) tests/check_prose_all.py --site docs

serve:
	$(PY) -m http.server -d docs 8765
