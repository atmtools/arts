#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


# SUBDIRS = Figs

artsdocdir = $(prefix)/artsdoc

thisdocdir = $(artsdocdir)/uguide

EXTRA_DIST = \
	agu.bst \
	chapter_template.tex \
	formalism.tex \
	los1d.tex \
	main.tex \
	measerr.tex \
	reduction.tex \
	references.bib \
	rte.tex \
	sensor.tex \
	wfuns_atm.tex \
	wfuns_sensor.tex
	
thisdoc_DATA = \
	main.dvi \
	main.ps

main.dvi: agu.bst \
	  chapter_template.tex \
	  formalism.tex \
	  los1d.tex \
	  main.tex \
	  measerr.tex \
	  reduction.tex \
	  references.bib \
	  rte.tex \
	  sensor.tex \
	  wfuns_atm.tex \
	  wfuns_sensor.tex \
	  version.tex
	latex  artsref 
	bibtex artsref
	latex  artsref 
	latex  artsref 

main.ps: main.dvi
	dvips main

CLEANFILES = *.aux *.log
