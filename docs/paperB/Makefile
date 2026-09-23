# Paper B build (skeleton). Every number in manuscript.tex is a macro from
# numbers.tex, which latex_tables.py regenerates from paperB/FROZEN.json;
# `check` verifies that the committed fragments are that regeneration.
PY ?= ../../.venv/bin/python
export SOURCE_DATE_EPOCH = 1788393600
export SOURCE_DATE_EPOCH_TEX_PRIMITIVES = 1

all: manuscript.pdf
	@$(MAKE) --no-print-directory check

manuscript.pdf: manuscript.tex numbers.tex tab_gates.tex tab_contrast.tex
	pdflatex -interaction=nonstopmode manuscript.tex
	pdflatex -interaction=nonstopmode manuscript.tex

freeze:
	$(PY) ../../paperB/freeze.py

tables:
	$(PY) latex_tables.py

figures:
	$(PY) figures.py

check:
	@$(PY) ../../paperB/freeze.py --check
	@$(PY) check_structure.py

clean:
	rm -f manuscript.{aux,bbl,blg,log,out,toc}

.PHONY: all check clean tables figures freeze
