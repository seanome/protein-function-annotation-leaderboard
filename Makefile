# Source of the tables: a commit of seanome/2024-kmerseek-analysis that holds
# results/leaderboard/, written there by `make leaderboard-tables`.
SOURCE_REPO ?= seanome/2024-kmerseek-analysis
SOURCE_COMMIT ?= 548ec5b8
PY ?= .venv/bin/python

.PHONY: tables site check serve

tables:
	$(PY) scripts/copy_tables.py --repo $(SOURCE_REPO) --commit $(SOURCE_COMMIT) --out data

site:
	$(PY) scripts/build_site.py --data data --out docs

check:
	$(PY) tests/check_site.py --site docs
	$(PY) tests/check_prose_all.py --site docs

serve:
	$(PY) -m http.server -d docs 8765
