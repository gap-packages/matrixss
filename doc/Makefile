TEX = tex
TEXFLAGS = -recorder
BUILDMAN = buildman.pe
BIBTEX = bibtex
BIBFLAGS = -terse
OUTPUT = manual.ps manual.bbl
FLS = manual.fls
GAPDOCPATH = /usr/share/gap/doc
MANUALINDEX = $(GAPDOCPATH)/manualindex
DVIPS = dvips

manual.dvi: manual.tex matrixss.tex code.tex manual.idx
	$(TEX) $(TEXFLAGS) $<

manual.idx: manual.tex manual.bbl manual.aux
	$(TEX) $(TEXFLAGS) $<
	$(TEX) $(TEXFLAGS) $<
	$(MANUALINDEX) $(basename $<)

manual.bbl: manual.aux report/schreiersims.bib 
	$(BIBTEX) $(BIBFLAGS) $(basename $<)

manual.aux: manual.tex
	$(TEX) $(TEXFLAGS) $<

code.tex: config.matrixss code.msk 
	$(BUILDMAN) -f $<

manual.ps: manual.dvi
	$(DVIPS) -o $@ $<

clean:
	rm -f *~ `grep OUTPUT $(FLS) | cut -f 2 -d ' '` $(FLS) $(OUTPUT)
	$(MAKE) -C report clean
