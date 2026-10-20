# Source of the tables: a commit of seanome/2024-kmerseek-analysis that holds
# results/leaderboard/, written there by `make leaderboard-tables`.
SOURCE_REPO ?= seanome/2024-kmerseek-analysis
SOURCE_COMMIT ?= 7eca428b
PILEUP_COMMIT ?= 39d81ec9
PY ?= .venv/bin/python
# Raw downloads of the organism lookups (NCBI taxdump, UniProt, GoaT, TimeTree, AlphaFold DB);
# about 120 MB, kept outside the repository.
ORGANISM_RAW ?= $(HOME)/data/organism-page

.PHONY: tables organism-lookups site check serve

tables:
	$(PY) scripts/copy_tables.py --repo $(SOURCE_REPO) --commit $(SOURCE_COMMIT) --pileup-commit $(PILEUP_COMMIT) --out data

# Lookups only, no tool runs. The alphafold step makes about 23 calls per pre-built species
# (about an hour for 1,952 species); re-running a step resumes where it stopped.
organism-lookups:
	$(PY) -I scripts/organism/fetch.py --raw $(ORGANISM_RAW) --step taxonomy
	$(PY) -I scripts/organism/fetch.py --raw $(ORGANISM_RAW) --step swissprot
	$(PY) -I scripts/organism/fetch.py --raw $(ORGANISM_RAW) --step proteomes
	$(PY) -I scripts/organism/fetch.py --raw $(ORGANISM_RAW) --step ebp
	$(PY) -I scripts/organism/fetch.py --raw $(ORGANISM_RAW) --step timetree
	$(PY) -I scripts/organism/lookups.py --raw $(ORGANISM_RAW) relatives
	$(PY) -I scripts/organism/fetch.py --raw $(ORGANISM_RAW) --step pairs
	$(PY) -I scripts/organism/fetch.py --raw $(ORGANISM_RAW) --step alphafold
	$(PY) -I scripts/organism/lookups.py --raw $(ORGANISM_RAW) tables

site:
	$(PY) scripts/build_site.py --data data --out docs

check:
	$(PY) scripts/claims.py --data data --check
	$(PY) tests/check_claims.py
	$(PY) tests/check_build.py --data data
	$(PY) tests/check_site.py --site docs
	$(PY) tests/check_pileup.py --site docs --data data
	$(PY) tests/check_explainers.py --site docs
	$(PY) tests/check_organism.py --site docs --data data --raw $(ORGANISM_RAW) --live
	node tests/check_organism_fallback.cjs docs
	$(PY) tests/check_prose_all.py --site docs

serve:
	$(PY) -m http.server -d docs 8765
