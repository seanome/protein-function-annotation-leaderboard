# Source of the tables: a commit of seanome/2024-kmerseek-analysis that holds
# results/leaderboard/, written there by `make leaderboard-tables`.
SOURCE_REPO ?= seanome/2024-kmerseek-analysis
SOURCE_COMMIT ?= 4a4b1171
PILEUP_COMMIT ?= 39d81ec9
PY ?= .venv/bin/python

.PHONY: tables site check serve rows release quotes-online

tables:
	$(PY) scripts/copy_tables.py --repo $(SOURCE_REPO) --commit $(SOURCE_COMMIT) --pileup-commit $(PILEUP_COMMIT) --out data

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
	$(PY) tests/check_release.py

# The release table in the shape of schema/results_row.json, checked against it.
rows:
	$(PY) scripts/release_rows.py --data data --out build/release/results.parquet
	$(PY) scripts/validate_rows.py build/release/results.parquet

# Freeze one version: Zenodo record with a DOI, mirrored to a Hugging Face dataset.
#   make release VERSION=v1 DRY_RUN=1                       # build and print; send nothing
#   ZENODO_TOKEN=... HF_TOKEN=... make release VERSION=v1   # needs huggingface_hub==1.33.0
release:
	@test -n "$(VERSION)" || (echo "VERSION is needed, e.g. make release VERSION=v1 DRY_RUN=1"; exit 1)
	$(PY) scripts/release.py --version $(VERSION) $(if $(DRY_RUN),--dry-run)

# Fetch every author's page again and check each quote is still there (needs the network).
quotes-online:
	$(PY) scripts/fetch_quote.py --check-all

serve:
	$(PY) -m http.server -d docs 8765
