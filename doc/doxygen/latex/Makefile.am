#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


artsdocdir = $(pkgdatadir)/doc

thisdocdir = $(artsdocdir)/doxygen/latex

#thisdoc_DATA         = refman.ps

#EXTRA_DIST           = $(wildcard *.tex *.eps doxygen.sty refman.dvi refman.ps)

CLEANFILES           = $(wildcard *.aux refman.idx refman.log refman.toc refman.ilg refman.ind)

MAINTAINERCLEANFILES = $(wildcard *.tex *.eps doxygen.sty refman.dvi refman.ps)

all:

# Cancel the implicit rule for .dvi and .ps files, since they
# introduces the wrong clean files: 
%.dvi: %.tex
%.ps: %.dvi


refman.dvi: $(wildcard *.tex)
	latex refman
	makeindex refman.idx
	latex refman

refman.ps: refman.dvi
	dvips -o refman.ps refman

refman.pdf: refman.ps
	ps2pdf refman.ps

