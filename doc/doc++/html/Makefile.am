#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

artsdocdir = $(pkgdatadir)/doc

thisdocdir = $(artsdocdir)/doc++/html

thisdoc_DATA = $(wildcard *.html *.gif *.class)

EXTRA_DIST  = $(wildcard *.html *.gif *.class)


# Unfortunately there has to be an explicit list here. 
# toexamine = $(top_srcdir)/src/*.h $(top_srcdir)/src/*.cc
# Does not work, because then the doc++ stuff is rebuilt each time the
# enduser installs the package, due to the automatically generated
# sources wsv.h, md.h, and md.cc
#
# Take toexamine from arts_SOURCES in file src/Makefile.am

toexamine = \
$(top_srcdir)/src/arts.h \
$(top_srcdir)/src/getopt.h \
$(top_srcdir)/src/sstream.h \
$(top_srcdir)/src/iotraits.h \
$(top_srcdir)/src/messages.h \
$(top_srcdir)/src/parameters.h \
$(top_srcdir)/src/vecmat.h \
$(top_srcdir)/src/exceptions.h \
$(top_srcdir)/src/file.h \
$(top_srcdir)/src/workspace.h \
$(top_srcdir)/src/wsv_group.h \
$(top_srcdir)/src/wsv_aux.h \
$(top_srcdir)/src/methods.h \
$(top_srcdir)/src/token.h \
$(top_srcdir)/src/parser.h \
$(top_srcdir)/src/make_vector.h \
$(top_srcdir)/src/make_array.h \
$(top_srcdir)/src/math_funcs.h \
$(top_srcdir)/src/atm_funcs.h \
$(top_srcdir)/src/absorption.h \
$(top_srcdir)/src/los.h \
$(top_srcdir)/src/hmatrix.h \
$(top_srcdir)/src/constants.cc \
$(top_srcdir)/src/getopt.c \
$(top_srcdir)/src/getopt1.c \
$(top_srcdir)/src/main.cc \
$(top_srcdir)/src/globals_1.cc \
$(top_srcdir)/src/globals_2.cc \
$(top_srcdir)/src/version.cc \
$(top_srcdir)/src/messages.cc \
$(top_srcdir)/src/parameters.cc \
$(top_srcdir)/src/file.cc \
$(top_srcdir)/src/workspace.cc \
$(top_srcdir)/src/workspace_aux.cc \
$(top_srcdir)/src/methods.cc \
$(top_srcdir)/src/methods_aux.cc \
$(top_srcdir)/src/token.cc \
$(top_srcdir)/src/parser.cc \
$(top_srcdir)/src/make_vector.cc \
$(top_srcdir)/src/math_funcs.cc \
$(top_srcdir)/src/atm_funcs.cc \
$(top_srcdir)/src/m_abs.cc \
$(top_srcdir)/src/m_io.cc \
$(top_srcdir)/src/m_los.cc \
$(top_srcdir)/src/m_wfs.cc \
$(top_srcdir)/src/absorption.cc \
$(top_srcdir)/src/hmatrix.cc

touse = \
$(toexamine) \
$(top_srcdir)/src/wsv.h \
$(top_srcdir)/src/md.h \
$(top_srcdir)/src/md.cc

index.html: $(toexamine)
	doc++ -v --filenames $(touse)

MAINTAINERCLEANFILES = *.html *.gif *.class
