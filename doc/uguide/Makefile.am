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

uguidedoc_DATA = uguide.ps uguide_pdf.pdf


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 	Don't forget to add new tex files to $(texfiles)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

texfiles = \
	agu.bst \
	references.bib \
	chapter_template.tex \
	concept.tex \
	formalism.tex \
	development.tex \
	matrix_vector.tex \
	formats.tex \
	los1d.tex \
	main.tex \
	measerr.tex \
	reduction.tex \
	rte.tex \
	sensor.tex \
	uguide.tex \
	uguide_pdf.tex \
	utilities.tex \
	wfuns_atm.tex \
	wfuns_sensor.tex \
	app_wsv.tex

EXTRA_DIST = \
	$(texfiles) \
	uguide.dvi \
	uguide.ps \
	uguide_pdf.pdf

if ARTS_GENDOCS

# The dependence of uguide.dvi on configure.in is there to make sure
# that the documentation is rebuilt when the release number in
# configure.in is changed. Dependence on version.tex is intentionally
# omitted, since this would lead to the documentation being rebuilt
# after each run of the `configure' script. version.tex.in does not
# show up explicitly in EXTRA_DIST, since this is distributed automatically.

uguide.dvi: $(top_srcdir)/configure.in $(texfiles) auto_version.tex
	latex  uguide
	bibtex uguide
	latex  uguide
	latex  uguide
	rm -f  *.aux

uguide.ps: uguide.dvi
	dvips -f uguide > uguide.ps

uguide_pdf.pdf: $(top_srcdir)/configure.in $(texfiles) auto_version.tex
	pdflatex uguide_pdf
	bibtex uguide_pdf
	pdflatex  uguide_pdf
	pdflatex  uguide_pdf
	rm -f  *.aux


#html: main.ps
#	latex2html -dir html -split 0 -no_math -local_icons main
#
# This does not work and I don't know why. Somehow neither figures nor
# equations show up in the output, also the redifined sectioning
# command (levela, etc.) confuse the translator.


CLEANFILES = *.aux *.log

MAINTAINERCLEANFILES = \
	uguide_pdf.bbl \
	uguide_pdf.blg \
	uguide.bbl \
	uguide.blg \
	uguide.dvi \
	uguide.ps \
	uguide_pdf.pdf \
	main.toc \
	auto_version.tex

else

uguide.dvi: $(texfiles)
	@echo
	@echo "=========================================================="
	@echo
	@echo "The documentation is only re-generated in maintainer mode."
	@echo "To activate maintainer mode call configure with the option"
	@echo "--enable-maintainer-mode"
	@echo
	@echo "=========================================================="
	@echo

uguide.ps: uguide.dvi
	@echo
	@echo "=========================================================="
	@echo
	@echo "The documentation is only re-generated in maintainer mode."
	@echo "To activate maintainer mode call configure with the option"
	@echo "--enable-maintainer-mode"
	@echo
	@echo "=========================================================="
	@echo

uguide_pdf.pdf:
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
