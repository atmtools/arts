# Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

bin_PROGRAMS = arts

# Set libraries:

# not needed at the moment: -ljpeg -lnsl -lm

arts_SOURCES = \
	absorption.cc \
	absorption.h \
	array.h \
	arts.h \
	atm_funcs.cc \
	atm_funcs.h \
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	auto_wsv_pointers.cc \
	complex.h \
	complex.cc \
	constants.cc \
	continua.cc \
	continua.h \
	exceptions.h \
	file.cc \
	file.h \
	getopt.c \
	getopt.h \
	getopt1.c \
	globals_2.cc \
	groups.cc \
	iotraits.h \
	lineshapes.cc \
	los.h \
	m_abs.cc \
	m_batch.cc \
	m_hdf.cc \
	m_io.cc \
	m_los.cc \
	m_wfs.cc \
	main.cc \
	make_array.h \
	make_vector.cc \
	make_vector.h \
	math_funcs.cc \
	math_funcs.h \
	matpackI.cc \
	matpackI.h \
	matpackII.h \
	messages.cc \
	messages.h \
	methods.cc \
	methods.h \
	methods_aux.cc \
	parameters.cc \
	parameters.h \
	parser.cc \
	parser.h \
	partition_function_data.cc \
	species_data.cc \
	sstream.h \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	wsv_aux.h \
	m_method_list.cc

methods.o: arts.h make_array.h auto_wsv.h methods.h auto_wsv_groups.h methods.cc
	$(CXX) $(DEFS) $(INCLUDES) -I.. $(AM_CPPFLAGS) $(CPPFLAGS) \
	       $(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -c methods.cc

species_data.o: arts.h make_array.h absorption.h species_data.cc
	$(CXX) $(DEFS) $(INCLUDES) -I.. $(AM_CPPFLAGS) $(CPPFLAGS) \
	       $(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -c species_data.cc

partition_function_data.o: arts.h make_array.h absorption.h partition_function_data.cc
	$(CXX) $(DEFS) $(INCLUDES) -I.. $(AM_CPPFLAGS) $(CPPFLAGS) \
	       $(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -c partition_function_data.cc


noinst_PROGRAMS = \
	make_auto_md_cc \
	make_auto_md_h \
	make_auto_wsv_groups_h \
	make_auto_wsv_h \
	make_auto_wsv_pointers_cc \
	test_matpack

test_matpack_SOURCES = \
	array.h \
	arts.h \
	make_array.h \
	make_vector.cc \
	make_vector.h \
	matpackI.cc \
	matpackI.h \
	matpackIII.h \
	mystring.h \
	test_matpack.cc

make_auto_wsv_groups_h_SOURCES = \
	arts.h \
	file.cc \
	file.h \
	groups.cc \
	make_auto_wsv_groups_h.cc \
	messages.cc \
	version.cc

make_auto_wsv_h_SOURCES = \
	absorption.h \
	arts.h \
	auto_wsv_groups.h \
	file.cc \
	file.h \
	groups.cc \
	los.h \
	make_auto_wsv_h.cc \
	matpackI.h \
	messages.cc \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	wsv_aux.h

make_auto_wsv_pointers_cc_SOURCES = \
	arts.h \
	auto_wsv_groups.h \
	file.cc \
	file.h \
	groups.cc \
	make_auto_wsv_pointers_cc.cc \
	matpackI.h \
	messages.cc \
	version.cc \
	workspace.cc \
	wsv_aux.h

make_auto_md_h_SOURCES = \
	absorption.h \
	arts.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	file.cc \
	file.h \
	globals_2.cc \
	groups.cc \
	los.h \
	make_array.h \
	make_auto_md_h.cc \
	matpackI.h \
	matpackIII.h \
	messages.cc \
	messages.h \
	methods.cc \
	methods.h \
	methods_aux.cc \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	wsv_aux.h

make_auto_md_cc_SOURCES = \
	absorption.h \
	arts.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	file.cc \
	file.h \
	globals_2.cc \
	groups.cc \
	los.h \
	make_array.h \
	make_auto_md_cc.cc \
	matpackI.h \
	matpackIII.h \
	messages.cc \
	messages.h \
	methods.cc \
	methods.h \
	methods_aux.cc \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	wsv_aux.h

auto_wsv.h: make_auto_wsv_h$(EXEEXT)
	$(top_builddir)/src/make_auto_wsv_h$(EXEEXT)

auto_wsv_groups.h: make_auto_wsv_groups_h$(EXEEXT)
	$(top_builddir)/src/make_auto_wsv_groups_h$(EXEEXT)

auto_wsv_pointers.cc: make_auto_wsv_pointers_cc$(EXEEXT)
	$(top_builddir)/src/make_auto_wsv_pointers_cc$(EXEEXT)

auto_md.h: make_auto_md_h$(EXEEXT)
	$(top_builddir)/src/make_auto_md_h$(EXEEXT)

auto_md.cc: make_auto_md_cc$(EXEEXT)
	$(top_builddir)/src/make_auto_md_cc$(EXEEXT)

BUILT_SOURCES = \
	auto_wsv_groups.h \
	auto_wsv.h \
	auto_wsv_pointers.cc \
	auto_md.h \
	auto_md.cc

CLEANFILES = \
	auto_wsv.h \
	auto_wsv_groups.h \
	auto_wsv_pointers.cc \
	auto_md.h \
	auto_md.cc

