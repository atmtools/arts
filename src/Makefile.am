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
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	auto_wsv_pointers.cc \
	check_input.cc \
	check_input.h \
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
	interpolation.cc \
	interpolation.h \
	iotraits.h \
	m_lineshapes.cc \
	m_agenda.cc \
        m_atmosphere.cc \
	m_basic_types.cc \
	m_cloudbox.cc \
	m_general.cc \
	m_hdf.cc \
	m_io.cc \
	m_ppath.cc \
	m_rte.cc \
	m_scatproperties.cc \
	m_xml.cc \
	m_clonesize.cc \
	m_copy.h \
	m_ignore.h \
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
	ppath.h \
	ppath.cc \
	poly_roots.h \
	poly_roots.cc \
	physics_funcs.cc \
	physics_funcs.h \
	rte.cc \
	rte.h \
	species_data.cc \
	sstream.h \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	wsv_aux.h \
	wsv_aux.cc \
	xml_io.h \
	xml_io.cc \
	agenda_class.h \
	agenda_class.cc \
	agenda_record.h \
	agenda_record.cc \
	agendas.cc \
	logic.h \
	logic.cc \
	lin_alg.h\
	lin_alg.cc\
	scatproperties.cc \
	scatproperties.h \
	supergeneric.h

arts_DEPENDENCIES= \
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	auto_wsv_pointers.cc

methods.o: arts.h make_array.h auto_wsv.h methods.h auto_wsv_groups.h methods.cc
	$(CXX) $(DEFS) $(INCLUDES) -I.. $(AM_CPPFLAGS) $(CPPFLAGS) \
	       $(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -c methods.cc

species_data.o: arts.h make_array.h absorption.h species_data.cc
	$(CXX) $(DEFS) $(INCLUDES) -I.. $(AM_CPPFLAGS) $(CPPFLAGS) \
	       $(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -c species_data.cc

arts: $(arts_DEPENDENCIES) $(arts_OBJECTS)
	@rm -f arts
	$(CXXLINK) $(arts_LDFLAGS) $(arts_OBJECTS) $(arts_LDADD) $(LIBS)

noinst_PROGRAMS = \
	make_auto_md_cc \
	make_auto_md_h \
	make_auto_wsv_groups_h \
	make_auto_wsv_h \
	make_auto_wsv_pointers_cc \
	test_matpack \
	test_poly \
	test_tensor \
	test_linalg \
	test_interpolation \
	test_scatproperties\
	test_xml

test_interpolation_SOURCES = \
	array.h \
	arts.h \
	make_array.h \
	make_vector.cc \
	make_vector.h \
	matpackI.cc \
	matpackI.h \
	matpackIII.h \
	matpackIV.h \
	matpackV.h \
	matpackVI.h \
	matpackVII.h \
	mystring.h \
	logic.h \
	logic.cc \
	check_input.h \
	check_input.cc \
	interpolation.h \
	interpolation.cc \
	test_interpolation.cc

test_linalg_SOURCES = \
	arts.h \
	matpackI.h \
	matpackI.cc \
	make_vector.h\
	make_vector.cc\
	lin_alg.cc \
	lin_alg.h \
	test_linalg.cc\
	logic.h\
	logic.cc

test_scatproperties_SOURCES = \
	arts.h \
	matpackI.h \
	matpackI.cc \
	make_vector.h\
	make_vector.cc\
	test_scatproperties.cc\
	scatproperties.h\
	scatproperties.cc\
	constants.cc\
	logic.h\
	math_funcs.cc\
	math_funcs.h\
	logic.cc

test_poly_SOURCES = \
	array.h \
	arts.h \
	make_array.h \
	matpackI.cc \
	matpackI.h \
	poly_roots.h \
	poly_roots.cc \
	test_poly.cc

test_tensor_SOURCES = \
	array.h \
	arts.h \
	make_array.h \
	make_vector.cc \
	make_vector.h \
	matpackI.cc \
	matpackI.h \
	matpackIII.h \
	matpackIV.h \
	matpackV.h \
	matpackVI.h \
	matpackVII.h \
	mystring.h \
	test_tensor.cc

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

test_xml_SOURCES = \
	xml_io.h \
	xml_io.cc \
	messages.cc \
	test_xml.cc

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
	wsv_aux.h \
	wsv_aux.cc

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
	wsv_aux.h \
	wsv_aux.cc


make_auto_wsv_h_DEPENDENCIES = auto_wsv_groups.h

make_auto_wsv_groups_h_DEPENDENCIES = 

make_auto_wsv_pointers_cc_DEPENDENCIES = auto_wsv_groups.h

make_auto_md_h_DEPENDENCIES = auto_wsv.h auto_wsv_groups.h

make_auto_md_cc_DEPENDENCIES = auto_wsv.h auto_wsv_groups.h

make_auto_wsv_h: $(make_auto_wsv_h_DEPENDENCIES) $(make_auto_wsv_h_OBJECTS)
	@rm -f make_auto_wsv_h
	$(CXXLINK) $(make_auto_wsv_h_LDFLAGS) \
	           $(make_auto_wsv_h_OBJECTS) \
	           $(make_auto_wsv_h_LDADD) $(LIBS)

make_auto_wsv_groups_h: $(make_auto_wsv_group_h_DEPENDENCIES) \
                        $(make_auto_wsv_groups_h_OBJECTS)
	@rm -f make_auto_wsv_groups_h
	$(CXXLINK) $(make_auto_wsv_groups_h_LDFLAGS) \
	           $(make_auto_wsv_groups_h_OBJECTS) \
	           $(make_auto_wsv_groups_h_LDADD) $(LIBS)

make_auto_wsv_pointers_cc: $(make_auto_wsv_pointers_cc_DEPENDENCIES) \
                           $(make_auto_wsv_pointers_cc_OBJECTS)
	@rm -f make_auto_wsv_pointers_cc
	$(CXXLINK) $(make_auto_wsv_pointers_cc_LDFLAGS) \
	           $(make_auto_wsv_pointers_cc_OBJECTS) \
	           $(make_auto_wsv_pointers_cc_LDADD) $(LIBS)

make_auto_md_h: $(make_auto_md_h_DEPENDENCIES) $(make_auto_md_h_OBJECTS)
	@rm -f make_auto_md_h
	$(CXXLINK) $(make_auto_md_h_LDFLAGS) \
	           $(make_auto_md_h_OBJECTS) \
	           $(make_auto_md_h_LDADD) $(LIBS)

make_auto_md_cc: $(make_auto_md_cc_DEPENDENCIES) $(make_auto_md_cc_OBJECTS)
	@rm -f make_auto_md_cc
	$(CXXLINK) $(make_auto_md_cc_LDFLAGS) \
	           $(make_auto_md_cc_OBJECTS) \
	           $(make_auto_md_cc_LDADD) $(LIBS)

auto_wsv.h: make_auto_wsv_h
	$(top_builddir)/src/make_auto_wsv_h

auto_wsv_groups.h: make_auto_wsv_groups_h
	$(top_builddir)/src/make_auto_wsv_groups_h

auto_wsv_pointers.cc: make_auto_wsv_pointers_cc
	$(top_builddir)/src/make_auto_wsv_pointers_cc

auto_md.h: make_auto_md_h
	$(top_builddir)/src/make_auto_md_h

auto_md.cc: make_auto_md_cc
	$(top_builddir)/src/make_auto_md_cc

CLEANFILES = \
	auto_wsv.h \
	auto_wsv_groups.h \
	auto_wsv_pointers.cc \
	auto_md.h \
	auto_md.cc

