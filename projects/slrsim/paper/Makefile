all: slrsim notes

slrsim: slrsim.tex slrsim.bib
	pdflatex slrsim; bibtex slrsim; pdflatex slrsim; pdflatex slrsim;
	pdflatex supp_slrsim; bibtex supp_slrsim; pdflatex supp_slrsim; pdflatex supp_slrsim;

notes: notes.tex
	pdflatex notes.tex