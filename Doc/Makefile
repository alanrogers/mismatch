#
all : main.ps 

main.ps : main.dvi
	dvips -o $@ main

main.dvi : main.aux
	latex main

main.aux : main.bbl misdefs.sty
	latex main

main.bbl : aux1
	bibtex main

aux1 : main.tex
	latex main
	touch aux1

