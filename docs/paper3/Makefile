# Paper III build. Every number in manuscript.tex is a macro from numbers.tex,
# which latex_tables.py regenerates from paper3/FROZEN.json; `check` verifies
# that the committed numbers.tex is that regeneration and that no literal
# result number has crept into the prose.
PY ?= ../../.venv/bin/python

all: manuscript.pdf
	@$(MAKE) --no-print-directory check

manuscript.pdf: manuscript.tex numbers.tex tab_verdict.tex tab_grid.tex tab_scenarios.tex references.bib
	pdflatex -interaction=nonstopmode manuscript.tex
	bibtex manuscript
	pdflatex -interaction=nonstopmode manuscript.tex
	pdflatex -interaction=nonstopmode manuscript.tex

si.pdf: si.tex numbers.tex references.bib
	pdflatex -interaction=nonstopmode si.tex
	bibtex si
	pdflatex -interaction=nonstopmode si.tex
	pdflatex -interaction=nonstopmode si.tex

# Regenerate the number macros and table fragments from FROZEN.json.
tables:
	$(PY) latex_tables.py

# Regenerate the display items from the frozen derived JSONs.
figures:
	$(PY) display_items.py

check:
	@$(PY) check_structure.py

clean:
	rm -f manuscript.{aux,bbl,blg,log,out,toc} si.{aux,bbl,blg,log,out,toc}

.PHONY: all check clean tables figures
