#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

artsdocdir = $(pkgdatadir)/doc

arts = $(top_srcdir)/src/arts

examplesdocdir = $(artsdocdir)/examples

arts_examples_in = \
	absorption_example.arts.in \
	amsu_example.arts.in \
	batch_example.arts.in \
	radiosonde_example.arts.in \
	rt_example.arts.in

arts_examples = \
	absorption_example.arts \
	amsu_example.arts \
	batch_example.arts \
	radiosonde_example.arts \
	rt_example.arts

examplesdoc_DATA = $(arts_examples)

EXTRA_DIST   = $(arts_examples_in)

CLEANFILES = *_example*.aa *_example*.ab *.rep *.out

check:	$(arts) $(arts_examples)
	$(arts) absorption_example.arts
	$(arts) amsu_example.arts
	$(arts) batch_example.arts
	$(arts) radiosonde_example.arts
	$(arts) rt_example.arts

