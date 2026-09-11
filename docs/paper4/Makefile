# Paper IV build. Every number in manuscript.tex is a macro from numbers.tex,
# which latex_tables.py regenerates from paper4/FROZEN.json; `check` verifies
# that the committed fragments are that regeneration and that no literal
# result number has crept into the prose.
PY ?= ../../.venv/bin/python
export SOURCE_DATE_EPOCH = 1788393600
export SOURCE_DATE_EPOCH_TEX_PRIMITIVES = 1

all: manuscript.pdf
	@$(MAKE) --no-print-directory check

manuscript.pdf: manuscript.tex numbers.tex tab_lightcurve.tex tab_grid.tex tab_matrix.tex references.bib
	pdflatex -interaction=nonstopmode manuscript.tex
	bibtex manuscript
	pdflatex -interaction=nonstopmode manuscript.tex
	pdflatex -interaction=nonstopmode manuscript.tex

freeze:
	$(PY) ../../paper4/freeze.py

tables:
	$(PY) latex_tables.py

figures:
	$(PY) display_items.py

check:
	@$(PY) ../../paper4/freeze.py --check
	@$(PY) check_structure.py

clean:
	rm -f manuscript.{aux,bbl,blg,log,out,toc}

.PHONY: all check clean tables figures freeze
