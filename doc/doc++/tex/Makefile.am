#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

EXTRA_DIST  = docxx.sty 

noinst_DATA = artsref.tex artsref.dvi artsref.ps

toexamine = $(top_srcdir)/src/*.h $(top_srcdir)/src/*.cc

artsref.tex: $(toexamine)
	doc++ -v --tex --output artsref.tex  $(toexamine)

artsref.dvi: artsref.tex
	latex artsref
	latex artsref

artsref.ps: artsref.dvi
	dvips artsref

CLEANFILES = artsref.aux artsref.log

MAINTAINERCLEANFILES = artsref.tex artsref.dvi artsref.ps
