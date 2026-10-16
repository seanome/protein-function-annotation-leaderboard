# Source of the tables: a commit of seanome/2024-kmerseek-analysis that holds
# results/leaderboard/, written there by `make leaderboard-tables`.
SOURCE_REPO ?= seanome/2024-kmerseek-analysis
SOURCE_COMMIT ?= 4a4b1171
PILEUP_COMMIT ?= 39d81ec9
# The alphabet page's tables: results/leaderboard/alphabets/, written by
# `make leaderboard-alphabets` (branch olgabot/leaderboard-alphabets until it merges).
ALPHABETS_COMMIT ?= 98e9bccb
PY ?= .venv/bin/python

.PHONY: tables alphabet-tables site check serve

tables:
	$(PY) scripts/copy_tables.py --repo $(SOURCE_REPO) --commit $(SOURCE_COMMIT) --pileup-commit $(PILEUP_COMMIT) --alphabets-commit $(ALPHABETS_COMMIT) --out data

alphabet-tables:
	$(PY) scripts/copy_tables.py --repo $(SOURCE_REPO) --commit $(SOURCE_COMMIT) --alphabets-commit $(ALPHABETS_COMMIT) --alphabets-only --out data

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
