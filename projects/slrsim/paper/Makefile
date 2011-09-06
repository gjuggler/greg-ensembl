all: slrsim notes

slrsim: slrsim.tex slrsim.bib
	pdflatex slrsim; bibtex slrsim; pdflatex slrsim; pdflatex slrsim;

notes: notes.tex
	pdflatex notes.tex