# Source of the tables: a commit of seanome/2024-kmerseek-analysis that holds
# results/leaderboard/, written there by `make leaderboard-tables`.
SOURCE_REPO ?= seanome/2024-kmerseek-analysis
SOURCE_COMMIT ?= 7eca428b
PILEUP_COMMIT ?= 39d81ec9
# results/leaderboard/query_strata/ (the query-side subsets), from the analysis repo's branch
# olgabot/leaderboard-query-strata. It branched before 7eca428b, so its results.parquet still has
# the typed divergence times; only query_strata/ is read from it.
QUERY_STRATA_COMMIT ?= 26a778cb
PY ?= .venv/bin/python

.PHONY: tables site check serve

tables:
	$(PY) scripts/copy_tables.py --repo $(SOURCE_REPO) --commit $(SOURCE_COMMIT) --pileup-commit $(PILEUP_COMMIT) --query-strata-commit $(QUERY_STRATA_COMMIT) --out data

site:
	$(PY) scripts/build_site.py --data data --out docs

check:
	$(PY) scripts/claims.py --data data --check
	$(PY) tests/check_claims.py
	$(PY) tests/check_build.py --data data
	$(PY) tests/check_site.py --site docs
	$(PY) tests/check_pileup.py --site docs --data data
	$(PY) tests/check_explainers.py --site docs
	$(PY) tests/check_prose_all.py --site docs

serve:
	$(PY) -m http.server -d docs 8765
