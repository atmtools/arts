#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

SUBDIRS = data

artsdocdir = $(pkgdatadir)/doc

arts = $(top_srcdir)/src/arts

examplesdocdir = $(artsdocdir)/examples

arts_examples = \
	absorption_example.arts \
	absorption_lookup_example.arts \
	cont_simple_example.arts

examplesdoc_DATA = $(arts_examples)

EXTRA_DIST   = *.arts

CLEANFILES = *.xml *_example*.aa *_example*.ab *.rep *.out

#check:	$(arts) $(arts_examples)
#	$(arts) absorption_example.arts
#	$(arts) absorption_lookup_example.arts
#	$(arts) cont_simple_example.arts

