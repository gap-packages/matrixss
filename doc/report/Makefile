LATEX = latex
PDFTEX = pdflatex
TEXFLAGS = --recorder
DVIPS = dvips
MPOST = mpost
MFLAGS = --recorder
BIBFLAGS = -terse
INPUTS := $(wildcard *.tex)
DVI := $(INPUTS:.tex=.dvi)
PS := $(INPUTS:.tex=.ps)
PDF := $(INPUTS:.tex=.pdf)
MPX := $(INPUTS:.tex=.mpx)
FLS := $(INPUTS:.tex=.fls)
BBL := $(INPUTS:.tex=.bbl)
AUX := $(INPUTS:.tex=.aux)

all: $(DVI)

$(MPX): %.mpx: %.mp
	$(MPOST) $(MFLAGS) $<

$(BBL): %.bbl: %.aux
	$(BIBTEX) $(BIBFLAGS) $(basename $<)

$(AUX): %.aux: %.tex
	$(LATEX) $(TEXFLAGS) $<

$(DVI): %.dvi: %.tex %.aux %.bbl
	$(LATEX) $(TEXFLAGS) $<
	$(LATEX) $(TEXFLAGS) $<

$(PS): %.ps: %.dvi %.tex
	$(DVIPS) -o $@ $<

$(PDF): %.pdf: %.tex
	$(PDFTEX) $(TEXFLAGS) $<
	$(PDFTEX) $(TEXFLAGS) $<

clean:
	rm -f *~ `grep OUTPUT $(FLS) | cut -f 2 -d ' '` $(FLS) $(PS) $(PDF) $(BBL)
