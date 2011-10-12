casper: casper.Rnw
	echo 'Sweave("casper.Rnw")' | R --vanilla --quiet
	R CMD pdflatex casper.tex
	bibtex casper
	R CMD pdflatex casper.tex
	R CMD pdflatex casper.tex
	rm casper.tex
	rm casper.log
	rm casper.aux
	rm casper.out
	rm casper.blg
	rm casper.bbl
	rm casper-fig*.pdf
	rm casper-fig*.eps
