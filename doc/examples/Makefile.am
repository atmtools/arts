#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

artsdocdir = $(pkgdatadir)/doc

examplesdocdir = $(artsdocdir)/examples

arts_examples = \
	cont.arts \
	amsu.arts

artsdata_examples_in = \
	auto_example1.arts.in \
	auto_full.arts.in

artsdata_examples = \
	auto_example1.arts \
	auto_full.arts

if ARTS_DATA_PATH
examplesdoc_DATA = $(arts_examples) $(artsdata_examples)
else
examplesdoc_DATA = $(arts_examples)
endif

EXTRA_DIST   = $(arts_examples) $(artsdata_examples_in)
