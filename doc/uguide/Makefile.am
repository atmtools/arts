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

thisdocdir = $(artsdocdir)/uguide

thisdoc_DATA = main.ps


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 	Don't forget to add new tex files to $(texfiles)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

texfiles = \
agu.bst \
references.bib \
chapter_template.tex \
concept.tex \
development.tex \
formalism.tex \
los1d.tex \
main.tex \
measerr.tex \
reduction.tex \
rte.tex \
sensor.tex \
wfuns_atm.tex \
wfuns_sensor.tex 

EXTRA_DIST = \
$(texfiles) \
main.dvi \
main.ps

# The dependence of main.dvi on configure.in is there to make sure
# that the documentation is rebuilt when the release number in
# configure.in is changed. Dependence on version.tex is intentionally
# omitted, since this would lead to the documentation being rebuilt
# after each run of the `configure' script. version.tex.in does not
# show up explicitly in EXTRA_DIST, since this is distributed automatically.

main.dvi: $(top_srcdir)/configure.in $(texfiles)
	latex  main 
	bibtex main
	latex  main 
	latex  main 

main.ps: main.dvi
	dvips main

#html: main.ps
#	latex2html -dir html -split 0 -no_math -local_icons main
#
# This does not work and I don't know why. Somehow neither figures nor
# equations show up in the output, also the redifined sectioning
# command (levela, etc.) confuse the translator.


CLEANFILES = *.aux *.log
