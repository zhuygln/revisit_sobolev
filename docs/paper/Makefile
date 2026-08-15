all:
	pdflatex -interaction=nonstopmode manuscript.tex
	bibtex manuscript
	pdflatex -interaction=nonstopmode manuscript.tex
	pdflatex -interaction=nonstopmode manuscript.tex

clean:
	rm -f manuscript.{aux,bbl,blg,log,out,toc,pdf}
