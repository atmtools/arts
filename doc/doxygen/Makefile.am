#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

SUBDIRS        = . html # latex

# The latex docuement is huge, and seems not to be too useful.

EXTRA_DIST  = doxyfied

noinst_DATA = doxyfied

# Unfortunately there has to be an explicit list here. 
# toexamine = $(top_srcdir)/src/*.h $(top_srcdir)/src/*.cc
# Does not work, because then the doc++ stuff is rebuilt each time the
# enduser installs the package, due to the automatically generated
# sources wsv.h, md.h, and md.cc
#
# Take toexamine from arts_SOURCES in file src/Makefile.am, but without
# wsv_groups.h, wsv.h, md.h, and md.cc. 
# Additionally add dependence on:
# - configure.in
# - make_wsv_groups.cc 
# - make_wsv_h.cc 
# - make_md_h.cc    
# The dependece on configure.in is there to make sure that the
# documentation is recompiled if the release number in that file is
# changed.

toexamine = \
$(top_srcdir)/configure.in \
$(top_srcdir)/src/make_wsv_groups_h.cc \
$(top_srcdir)/src/make_wsv_h.cc \
$(top_srcdir)/src/make_md_h.cc \
$(top_srcdir)/src/arts.h \
$(top_srcdir)/src/getopt.h \
$(top_srcdir)/src/sstream.h \
$(top_srcdir)/src/iotraits.h \
$(top_srcdir)/src/messages.h \
$(top_srcdir)/src/parameters.h \
$(top_srcdir)/src/vecmat.h \
$(top_srcdir)/src/exceptions.h \
$(top_srcdir)/src/file.h \
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
$(top_srcdir)/src/vecmat.cc \
$(top_srcdir)/src/globals_2.cc \
$(top_srcdir)/src/version.cc \
$(top_srcdir)/src/messages.cc \
$(top_srcdir)/src/parameters.cc \
$(top_srcdir)/src/file.cc \
$(top_srcdir)/src/groups.cc \
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
$(top_srcdir)/src/m_hmatrix.cc \
$(top_srcdir)/src/absorption.cc \
$(top_srcdir)/src/species_data.cc \
$(top_srcdir)/src/partition_function_data.cc \
$(top_srcdir)/src/lineshapes.cc 


doxyfied: $(toexamine)
	doxygen Doxyfile > /dev/null
	touch doxyfied
