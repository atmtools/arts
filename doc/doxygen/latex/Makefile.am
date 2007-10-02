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

CLEANFILES           = $(wildcard *.aux) refman.idx refman.log refman.toc refman.ilg refman.ind)

MAINTAINERCLEANFILES = $(wildcard *.tex *.eps *.md5 *.ttf) doxygen.sty refman.dvi refman.ps

all:

if MAINTAINER_MODE

# Cancel the implicit rule for .dvi and .ps files, since they
# introduces the wrong clean files: 
%.dvi: %.tex
%.ps: %.dvi

ps: refman.ps

pdf: refman.pdf

ps_2on1: refman_2on1.ps

pdf_2on1: refman_2on1.pdf

refman.ps: refman.dvi
	dvips -o refman.ps refman.dvi

refman.pdf: refman.ps
	ps2pdf refman.ps refman.pdf

refman.dvi: refman.tex doxygen.sty
	echo "Running latex..."
	latex refman.tex
	echo "Running makeindex..."
	makeindex refman.idx
	echo "Rerunning latex...."
	latex refman.tex
	latex_count=5 ; \
	while egrep -s 'Rerun (LaTeX|to get cross-references right)' refman.log && [ $$latex_count -gt 0 ] ;\
	    do \
	      echo "Rerunning latex...." ;\
	      latex refman.tex ;\
	      latex_count=`expr $$latex_count - 1` ;\
	    done

refman_2on1.ps: refman.ps
	psnup -2 refman.ps >refman_2on1.ps

refman_2on1.pdf: refman_2on1.ps
	ps2pdf refman_2on1.ps refman_2on1.pdf

endif

