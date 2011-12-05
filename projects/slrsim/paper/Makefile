all: slrsim notes concat

slrsim: slrsim.tex slrsim.bib
	pdflatex slrsim; bibtex slrsim; pdflatex slrsim; pdflatex slrsim;
	pdflatex supp_slrsim; bibtex supp_slrsim; pdflatex supp_slrsim; pdflatex supp_slrsim;
	pdflatex slrsim_summary;

notes: notes.tex
	pdflatex notes.tex

concat: slrsim.pdf supp_slrsim.pdf
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=slrsim_full.pdf slrsim.pdf supp_slrsim.pdf