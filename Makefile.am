#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

EXTRA_DIST = reconf configure make_parallel

DEFAULT_SUBDIRS = src aii ami m4

if ARTS_GENDOCS
SUBDIRS    = $(DEFAULT_SUBDIRS) doc
else
SUBDIRS    = $(DEFAULT_SUBDIRS)
endif

# This stuff does not work:
#
#documentation:
#	cd $(top_srcdir)/doc && $(MAKE)
#
#release:
#	$(MAKE) documentation
#	$(MAKE) dist
