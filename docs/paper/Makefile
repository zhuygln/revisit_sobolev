all:
	pdflatex -interaction=nonstopmode manuscript.tex
	bibtex manuscript
	pdflatex -interaction=nonstopmode manuscript.tex
	pdflatex -interaction=nonstopmode manuscript.tex
	@$(MAKE) --no-print-directory check

# Content check. LaTeX compiles a stub section perfectly happily, and the page
# count can stay fixed while a whole appendix disappears -- both happened.
check:
	@python3 check_structure.py

clean:
	rm -f manuscript.{aux,bbl,blg,log,out,toc,pdf}

.PHONY: all check clean
