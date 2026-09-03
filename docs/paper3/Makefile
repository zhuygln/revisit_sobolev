# Paper III build. Every number in manuscript.tex is a macro from numbers.tex,
# which latex_tables.py regenerates from paper3/FROZEN.json; `check` verifies
# that the committed numbers.tex is that regeneration and that no literal
# result number has crept into the prose.
PY ?= ../../.venv/bin/python
# Reproducible PDFs: pdflatex takes its dates (and the trailer /ID) from this.
export SOURCE_DATE_EPOCH = 1788393600
export SOURCE_DATE_EPOCH_TEX_PRIMITIVES = 1

all: manuscript.pdf si.pdf
	@$(MAKE) --no-print-directory check

manuscript.pdf: manuscript.tex numbers.tex tab_verdict.tex tab_grid.tex tab_scenarios.tex references.bib
	pdflatex -interaction=nonstopmode manuscript.tex
	bibtex manuscript
	pdflatex -interaction=nonstopmode manuscript.tex
	pdflatex -interaction=nonstopmode manuscript.tex

si.pdf: si.tex numbers.tex si_tab_robustness.tex si_tab_points.tex si_tab_chain.tex si_tab_tscale.tex si_tab_syserr.tex
	pdflatex -interaction=nonstopmode si.tex
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
