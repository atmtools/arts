#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


SUBDIRS =

artsdocdir = $(pkgdatadir)/doc

referencedocdir = $(artsdocdir)/reference

referencedoc_DATA = reference.pdf

artsexecutable = $(top_builddir)/src/arts

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 	Don't forget to add new tex files to $(texfiles)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

texfiles = \
	app_groups.tex \
	app_methods.tex \
	main.tex \
	reference.tex

autotexfiles = \
	auto_groups.tex \
	auto_methods.tex

EXTRA_DIST = \
	$(texfiles) \
	$(autotexfiles) \
	reference.pdf \
	arts_groups_to_latex.sh \
	arts_methods_to_latex.sh


if MAINTAINER_MODE

# version.tex.in does not show up explicitly in EXTRA_DIST, since this is
# distributed automatically.
# 
# Here also, auto_groups.tex is added independently.

reference.pdf: $(top_srcdir)/configure.in $(texfiles) $(autotexfiles) auto_version.tex
	pdflatex reference
	makeindex reference.idx
#	bibtex reference
#	pdflatex  reference
	pdflatex  reference
	pdflatex  reference | tee make.log


auto_groups.tex: arts_groups_to_latex.sh $(artsexecutable)
	./arts_groups_to_latex.sh $(artsexecutable) > auto_groups.tex

auto_methods.tex: arts_methods_to_latex.sh $(artsexecutable)
	./arts_methods_to_latex.sh $(artsexecutable) > auto_methods.tex

CLEANFILES = *.aux \
	     *.log \
	     make.log \
	     reference.bbl \
	     reference.blg \
	     reference.idx \
	     reference.ilg \
	     reference.ind \
	     reference.out \
	     reference.pdf \
	     reference.toc

MAINTAINERCLEANFILES = \
	$(autotexfiles) \
	auto_version.tex

else

reference.pdf:
	@echo
	@echo "=========================================================="
	@echo
	@echo "The documentation is only re-generated in maintainer mode."
	@echo "To activate maintainer mode call configure with the option"
	@echo "--enable-maintainer-mode"
	@echo
	@echo "=========================================================="
	@echo

auto_groups.tex:
	@echo
	@echo "=========================================================="
	@echo
	@echo "The documentation is only re-generated in maintainer mode."
	@echo "To activate maintainer mode call configure with the option"
	@echo "--enable-maintainer-mode"
	@echo
	@echo "=========================================================="
	@echo

auto_methods.tex:
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
