# Source of the tables: a commit of seanome/2024-kmerseek-analysis that holds
# results/leaderboard/, written there by `make leaderboard-tables`.
SOURCE_REPO ?= seanome/2024-kmerseek-analysis
SOURCE_COMMIT ?= 4a4b1171
PILEUP_COMMIT ?= 39d81ec9
PY ?= .venv/bin/python

.PHONY: tables hardware site check serve

# Needs `ssh sherlock` (open the socket first with `ssh -fN sherlock`).
hardware:
	$(PY) scripts/read_hardware.py --data data

tables:
	$(PY) scripts/copy_tables.py --repo $(SOURCE_REPO) --commit $(SOURCE_COMMIT) --pileup-commit $(PILEUP_COMMIT) --out data

site:
	$(PY) scripts/build_site.py --data data --out docs

check:
	$(PY) scripts/claims.py --data data --check
	$(PY) tests/check_site.py --site docs
	$(PY) tests/check_prose_all.py --site docs

serve:
	$(PY) -m http.server -d docs 8765
