# Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

bin_PROGRAMS   = arts

# Set libraries:
#LDADD = -lhdf5 -lz 
LDADD = -lmfhdf -ldf -ljpeg -lz
# not needed at the moment: -ljpeg -lnsl -lm

#$(top_srcdir)/config.h \
# Does not work (source file in separate dir)

arts_SOURCES   = \
auto_wsv_groups.h \
auto_wsv.h \
auto_md.cc \
arts.h \
getopt.h \
sstream.h \
iotraits.h \
messages.h \
parameters.h \
exceptions.h \
file.h \
wsv_aux.h \
methods.h \
token.h \
parser.h \
make_vector.h \
make_array.h \
math_funcs.h \
atm_funcs.h \
absorption.h \
continua.h \
los.h \
array.h \
matpackI.h \
matpackII.h \
constants.cc \
getopt.c \
getopt1.c \
main.cc \
globals_2.cc \
version.cc \
messages.cc \
parameters.cc \
file.cc \
groups.cc \
workspace.cc \
workspace_aux.cc \
methods.cc \
methods_aux.cc \
token.cc \
parser.cc \
make_vector.cc \
math_funcs.cc \
atm_funcs.cc \
m_abs.cc \
m_batch.cc \
m_io.cc \
m_los.cc \
m_wfs.cc \
absorption.cc \
continua.cc \
auto_wsv_pointers.cc \
species_data.cc \
partition_function_data.cc \
lineshapes.cc \
matpackI.cc

arts_CLEAN=auto_md.cc

# Is this still valid? (PE 2001-09-19)
#
# When making changes to the list of sources above, please remember to
# also change the lists in the files:
#       doc/doxygen/Makefile.am


noinst_PROGRAMS = \
test_matpack \
make_auto_wsv_groups_h make_auto_wsv_h make_auto_wsv_pointers_cc make_auto_md_h make_auto_md_cc

test_matpack_SOURCES = \
arts.h matpackI.h array.h make_array.h mystring.h make_vector.h \
test_matpack.cc matpackI.cc make_vector.cc

make_auto_wsv_groups_h_SOURCES = \
arts.h file.h \
version.cc file.cc messages.cc groups.cc make_auto_wsv_groups_h.cc

make_auto_wsv_h_SOURCES = \
arts.h file.h absorption.h los.h hmatrix.h \
auto_wsv_groups.h wsv_aux.h \
version.cc groups.cc \
messages.cc file.cc workspace.cc workspace_aux.cc make_auto_wsv_h.cc

make_auto_wsv_pointers_cc_SOURCES = \
arts.h file.h wsv_aux.h \
version.cc groups.cc messages.cc file.cc workspace.cc make_auto_wsv_pointers_cc.cc

make_auto_md_h_SOURCES = \
arts.h auto_wsv.h token.h messages.h file.h \
absorption.h los.h hmatrix.h \
auto_wsv_groups.h wsv_aux.h \
methods.h make_array.h \
version.cc globals_2.cc \
token.cc messages.cc file.cc groups.cc workspace.cc workspace_aux.cc \
methods.cc methods_aux.cc make_auto_md_h.cc

make_auto_md_cc_SOURCES = \
arts.h auto_wsv.h token.h messages.h file.h \
absorption.h los.h hmatrix.h \
auto_wsv_groups.h wsv_aux.h \
methods.h make_array.h \
version.cc globals_2.cc \
token.cc messages.cc file.cc groups.cc workspace.cc workspace_aux.cc \
methods.cc methods_aux.cc make_auto_md_cc.cc

# make_auto_wsv_h_DEPENDENCIES = auto_wsv_groups.h

# make_auto_wsv_pointers_cc_DEPENDENCIES  = auto_wsv_groups.h auto_wsv.h

# make_auto_md_h_DEPENDENCIES  = auto_wsv.h

auto_wsv_groups.h: make_auto_wsv_groups_h
	$(top_builddir)/src/make_auto_wsv_groups_h

auto_wsv.h: make_auto_wsv_h
	$(top_builddir)/src/make_auto_wsv_h

auto_wsv_pointers.cc: make_auto_wsv_pointers_cc
	$(top_builddir)/src/make_auto_wsv_pointers_cc

auto_md.h: make_auto_md_h
	$(top_builddir)/src/make_auto_md_h

auto_md.cc: make_auto_md_cc
	$(top_builddir)/src/make_auto_md_cc

BUILT_SOURCES = \
auto_wsv_groups.h auto_wsv.h auto_wsv_pointers.cc auto_md.h auto_md.cc
