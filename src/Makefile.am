# Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

INCLUDES       = -I$(top_srcdir)/src/tnt
bin_PROGRAMS   = arts
SUBDIRS        = tnt

#$(top_srcdir)/config.h \
# Does not work (source file in separate dir)

arts_SOURCES   = \
wsv.h \
md.h \
md.cc \
arts.h \
getopt.h \
sstream.h \
iotraits.h \
messages.h \
parameters.h \
vecmat.h \
exceptions.h \
file.h \
workspace.h \
wsv_group.h \
wsv_aux.h \
methods.h \
token.h \
parser.h \
make_vector.h \
make_array.h \
math_funcs.h \
atm_funcs.h \
absorption.h \
los.h \
constants.cc \
getopt.c \
getopt1.c \
main.cc \
globals_1.cc \
globals_2.cc \
version.cc \
messages.cc \
parameters.cc \
file.cc \
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
m_io.cc \
m_los.cc \
m_wfs.cc \
absorption.cc

# When making changes to the list of sources above, please remember to
# also change the lists in the files:
# 	doc/doc++/html/Makefile.am and 
# 	doc/doc++/tex/Makefile.am


noinst_PROGRAMS = make_wsv_h make_md_h

make_wsv_h_SOURCES = \
arts.h vecmat.h messages.h file.h absorption.h los.h \
workspace.h wsv_group.h wsv_aux.h \
globals_1.cc \
messages.cc file.cc workspace.cc workspace_aux.cc make_wsv_h.cc version.cc

make_md_h_SOURCES = \
arts.h vecmat.h wsv.h token.h messages.h file.h \
absorption.h los.h \
workspace.h wsv_group.h wsv_aux.h \
methods.h make_array.h \
globals_1.cc globals_2.cc \
token.cc messages.cc file.cc workspace.cc workspace_aux.cc \
methods.cc methods_aux.cc make_md_h.cc version.cc

make_md_h_DEPENDENCIES = wsv.h

wsv.h: make_wsv_h
	./make_wsv_h

md.h: make_md_h
	./make_md_h

BUILT_SOURCES = wsv.h md.h
