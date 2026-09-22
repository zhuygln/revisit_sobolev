# The education book. rtedu is deliberately not installed: it runs from
# $(CURDIR)/src (PI amendment 5). Every displayed number comes from code:
# the notebooks write data/generated/chNN.json, merged into data/results.json
# (transactional: both are removed before the notebooks run, amendment 4).
PY ?= ../.venv/bin/python
QUARTO ?= ../.venv/bin/quarto
export PYTHONPATH := $(CURDIR)/src
export QUARTO_PYTHON := $(abspath $(PY))
export SOURCE_DATE_EPOCH = 1788393600
export MPLBACKEND = Agg
NOTEBOOKS := $(sort $(wildcard notebooks/[0-9][0-9]_*.ipynb))

education: education-test education-notebooks education-videos education-render
education-fast: education-test education-notebooks education-render

education-test:
	$(PY) -m pytest tests -q

education-notebooks:
	rm -rf data/generated data/results.json notebooks/executed
	mkdir -p notebooks/executed
	for nb in $(NOTEBOOKS); do \
	  $(PY) -m jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 \
	    --output-dir notebooks/executed --output $$(basename $$nb) $$nb || exit 1; done
	$(PY) -c "from rtedu import results; print('merged', list(results.merge()))"
	$(PY) scripts/make_static_figures.py

education-videos:
	for s in scripts/make_video_*.py; do $(PY) $$s || exit 1; done

education-render:
	$(PY) scripts/check_no_literals.py
	$(QUARTO) render

education-check: education-test
	$(PY) scripts/check_no_literals.py
	$(PY) scripts/check_determinism.py

clean:
	rm -rf _book _freeze .quarto notebooks/executed data/generated

.PHONY: education education-fast education-test education-notebooks education-videos education-render education-check clean
