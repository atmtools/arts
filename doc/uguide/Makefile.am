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

uguidedoc_DATA = uguide.pdf

artsexecutable = $(top_builddir)/src/arts

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 	Don't forget to add new tex files to $(texfiles)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

texfiles = \
	absorption.tex \
	agendas.tex \
	agu.bst \
	app_units.tex \
	app_wsv.tex \
	atmosphere.tex \
	auto_version.tex \
	batch.tex \
	chapter_template.tex \
	clouds.tex \
	concept.tex \
	development.tex \
	fm_definitions.tex \
	formalism.tex \
	formats.tex \
	integration.tex \
	interpolation.tex \
	introduction.tex \
	lin_alg.tex \
	main.tex \
	matrix_vector.tex \
	measerr.tex \
	montecarlo.tex \
	old_absorption.tex \
	polarization.tex \
	ppath.tex \
	reduction.tex \
	references.bib \
	rte.tex \
	rte_theory.tex \
	scattering.tex \
	sensor.tex \
	surface.tex \
	symbol_defs.tex \
	uguide.tex \
	utilities.tex \
	wfuns_atm.tex \
	wfuns_sensor.tex

autotexfiles = 

EXTRA_DIST = \
	$(texfiles) \
	uguide.pdf


if MAINTAINER_MODE

# version.tex.in does not show up explicitly in EXTRA_DIST, since this is
# distributed automatically.

uguide.pdf: $(top_srcdir)/configure.in $(texfiles) auto_version.tex
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: pdflatex uguide && \
	makeindex uguide.idx && \
	BSTINPUTS=.:$(srcdir):${BSTINPUTS}: BIBINPUTS=.:$(srcdir):${BIBINPUTS}: bibtex uguide && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: pdflatex uguide && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: pdflatex uguide && \
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: pdflatex uguide | tee make.log || \
	(echo "Errors encountered while building user guide."; rm -f uguide.pdf)


CLEANFILES = $(wildcard *.aux) \
	     $(wildcard *.log) \
	     make.log \
	     uguide.bbl \
	     uguide.blg \
	     uguide.idx \
	     uguide.ilg \
	     uguide.ind \
	     uguide.out \
	     uguide.pdf \
	     uguide.toc

MAINTAINERCLEANFILES = \
	auto_version.tex

else

uguide.pdf:
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
