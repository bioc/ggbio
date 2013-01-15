# Makefile to use knitr for package vignettes

# put all PDF targets here, separated by spaces
PDFS= ggbio.pdf

all: $(PDFS) 

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.toc
	rm -f figures/ggbio-*	ggbio-*
	rm -f ggbio.Rnw
	cp ggbio.Rnw.stub ggbio.Rnw
	
%.pdf: %.Rnw
	cp ggbio.Rnw.knitr ggbio.Rnw
	$(R_HOME)/bin/Rscript -e "library(knitr); knit2pdf('$*.Rnw')"

