#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


SUBDIRS = Figs

artsdocdir = $(pkgdatadir)/doc

uguidedocdir = $(artsdocdir)/uguide

uguidedoc_DATA = arts_user.pdf arts_developer.pdf arts_theory.pdf

usertexfiles = $(wildcard *.tex)

developertexfiles = $(wildcard *.tex)

theorytexfiles = $(wildcard *.tex)

autotexfiles = 

EXTRA_DIST = \
	$(usertexfiles) \
	$(developertexfiles) \
	$(theorytexfiles) \
	arts_user.pdf \
	arts_developer.pdf \
	arts_theory.pdf


if MAINTAINER_MODE

# version.tex.in does not show up explicitly in EXTRA_DIST, since this is
# distributed automatically.

arts_user.pdf: $(top_srcdir)/configure.in $(texfiles) auto_version.tex
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_user && \
	makeindex arts_user.idx && \
	BSTINPUTS=.:$(srcdir):${BSTINPUTS}: BIBINPUTS=.:$(srcdir):${BIBINPUTS}: save_size=15000 bibtex arts_user && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_user && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_user && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_user | tee make.log || \
	(echo "Errors encountered while building user guide."; rm -f arts_user.pdf)

arts_developer.pdf: $(top_srcdir)/configure.in $(texfiles) auto_version.tex
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_developer && \
	makeindex arts_developer.idx && \
	BSTINPUTS=.:$(srcdir):${BSTINPUTS}: BIBINPUTS=.:$(srcdir):${BIBINPUTS}: save_size=15000 bibtex arts_developer && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_developer && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_developer && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_developer | tee make.log || \
	(echo "Errors encountered while building developer guide."; rm -f arts_developer.pdf)

arts_theory.pdf: $(top_srcdir)/configure.in $(texfiles) auto_version.tex
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_theory && \
	makeindex arts_theory.idx && \
	BSTINPUTS=.:$(srcdir):${BSTINPUTS}: BIBINPUTS=.:$(srcdir):${BIBINPUTS}: save_size=15000 bibtex arts_theory && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_theory && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_theory && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: save_size=15000 pdflatex arts_theory | tee make.log || \
	(echo "Errors encountered while building theory guide."; rm -f arts_theory.pdf)


CLEANFILES = $(wildcard *.aux) \
	     $(wildcard *.log) \
	     make.log \
	     $(wildcard *.bbl) \
	     $(wildcard *.blg) \
	     $(wildcard *.idx) \
	     $(wildcard *.ilg) \
	     $(wildcard *.ind) \
	     $(wildcard *.out) \
	     $(wildcard *.pdf) \
	     $(wildcard *.toc)

MAINTAINERCLEANFILES = \
	auto_version.tex

else

arts_user.pdf:
arts_developer.pdf:
arts_theory.pdf:
	@echo
	@echo "=========================================================="
	@echo
	@echo "The documentation is only re-generated in maintainer mode."
	@echo "To activate maintainer mode call configure with the option"
	@echo "--enable-maintainer-mode"
	@echo
	@echo "=========================================================="
	@echo

endif
