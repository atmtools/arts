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

arts_SOURCES   = \
arts.h \
getopt.h \
messages.h parameters.h \
vecmat.h \
exceptions.h file.h \
workspace.h \
methods.h \
token.h \
parser.h \
md.h \
math_funcs.h \
getopt.c \
getopt1.c \
main.cc \
version.cc constants.cc \
messages.cc parameters.cc \
file.cc \
workspace.cc \
methods.cc \
token.cc \
parser.cc \
md.cc \
math_funcs.cc \
absorption.cc

noinst_PROGRAMS = make_wsv_h make_md_h

make_wsv_h_SOURCES = \
arts.h vecmat.h messages.h file.h workspace.h \
messages.cc file.cc workspace.cc make_wsv_h.cc

make_md_h_SOURCES = \
arts.h vecmat.h wsv.h token.h messages.h file.h workspace.h methods.h \
token.cc messages.cc file.cc workspace.cc methods.cc make_md_h.cc

wsv.h: make_wsv_h
	./make_wsv_h

md.h: make_md_h
	./make_md_h

BUILT_SOURCES = wsv.h md.h
