# Source of the tables: a commit of seanome/2024-kmerseek-analysis that holds
# results/leaderboard/, written there by `make leaderboard-tables`.
SOURCE_REPO ?= seanome/2024-kmerseek-analysis
SOURCE_COMMIT ?= 7eca428b
PILEUP_COMMIT ?= 39d81ec9
PY ?= .venv/bin/python

.PHONY: tables site check serve score rescore truth

# Score one submission folder the way the pull-request workflow does: make score SUB=submissions/<name>
score:
	$(PY) scripts/score_submission.py --submission $(SUB) --out scored/$(notdir $(SUB))

# At a release: rebuild truth_contributed from the curation queue, then score every
# submission again, so each one is graded on the newly accepted wet-lab rows too.
truth:
	$(PY) scripts/truth_contributed.py build

rescore: truth
	for d in submissions/*/; do $(PY) scripts/score_submission.py --submission $$d --out scored/$$(basename $$d) || exit 1; done

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
	$(PY) tests/check_scorer.py
	$(PY) tests/check_contribute.py
	$(PY) tests/check_prose_all.py --site docs

serve:
	$(PY) -m http.server -d docs 8765
