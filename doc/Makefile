TEX = tex
TEXFLAGS = -recorder
PDFTEX = pdftex
GAP = ~/gap
BUILDMAN = $(GAP)/etc/buildman.pe
BIBTEX = bibtex
BIBFLAGS = -terse
OUTPUT = manual.ps manual.bbl manual.ind manual.blg manual.ilg code.tex \
	 notfound manual.pdf
FLS = manual.fls
MANUALINDEX = $(GAP)/doc/manualindex
DVIPS = dvips

manual: manual.dvi manual.pdf

manual.pdf: manual.tex matrixss.tex code.tex manual.ind
	$(PDFTEX) $(TEXFLAGS) $<

manual.dvi: manual.tex matrixss.tex code.tex manual.ind
	$(TEX) $(TEXFLAGS) $<

manual.ind: manual.tex manual.bbl manual.aux
	$(TEX) $(TEXFLAGS) $<
	$(TEX) $(TEXFLAGS) $<
	$(MANUALINDEX) $(basename $<)

manual.bbl: manual.aux report/schreiersims.bib 
	$(BIBTEX) $(BIBFLAGS) $(basename $<)

manual.aux: manual.tex
	$(TEX) $(TEXFLAGS) $<

code.tex: config.matrixss code.msk ../lib/*.gd ../lib/*.gi
	$(BUILDMAN) -f $<

manual.ps: manual.dvi
	$(DVIPS) -o $@ $<

clean:
	rm -f *~ `grep OUTPUT $(FLS) | cut -f 2 -d ' '` $(FLS) $(OUTPUT)
	$(MAKE) -C report clean
