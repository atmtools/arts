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
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	check_input.cc \
	check_input.h \
	cloudbox.h \
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
	gridded_fields.cc\
	gridded_fields.h\
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
	m_optproperties.cc \
	m_physics.cc \
	m_ppath.cc \
	m_refraction.cc \
	m_rte.cc \
	m_scatrte.cc \
	m_clonesize.cc \
	m_absorption.cc \
	m_copy.h \
	m_ignore.h \
	m_sensor.cc \
	m_xml.cc \
	m_xml.h \
	main.cc \
	make_array.h \
	make_vector.cc \
	make_vector.h \
	math_funcs.cc \
	math_funcs.h \
	matpackI.cc \
	matpackI.h \
	matpackII.h \
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
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
	refraction.h \
	refraction.cc \
	rte.cc \
	rte.h \
	special_interp.cc \
	special_interp.h \
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
	xml_io_private.h \
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
	supergeneric.h \
        zeemanproperties.cc \
	gas_abs_lookup.h \
	gas_abs_lookup.cc

methods.o: arts.h make_array.h auto_wsv.h methods.h auto_wsv_groups.h methods.cc
	$(CXX) $(DEFS) $(INCLUDES) -I.. $(AM_CPPFLAGS) $(CPPFLAGS) \
	       $(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -c methods.cc

species_data.o: arts.h make_array.h absorption.h species_data.cc
	$(CXX) $(DEFS) $(INCLUDES) -I.. $(AM_CPPFLAGS) $(CPPFLAGS) \
	       $(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -c species_data.cc

noinst_PROGRAMS = \
	make_auto_md_cc \
	make_auto_md_h \
	make_auto_wsv_groups_h \
	make_auto_wsv_h \
	make_auto_wsv_pointers_cc \
	test_binaryio \
	test_matpack \
	test_poly \
	test_tensor \
	test_linalg \
	test_interpolation \
	test_scatproperties\
	test_xml \
	test_complex \
	test_zeemanproperties

test_binaryio_SOURCES = \
	absorption.cc \
	absorption.h \
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	partition_function_data.cc \
	matpackI.cc \
	matpackI.h \
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
	xml_io.h \
	xml_io_private.h \
	xml_io.cc \
	messages.cc \
	species_data.cc \
	test_binaryio.cc

test_interpolation_SOURCES = \
	array.h \
	constants.cc \
	matpackI.cc \
	matpackI.h \
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
	interpolation.h \
	interpolation.cc \
	logic.cc \
	test_interpolation.cc \
	make_vector.cc \
	make_vector.h

test_linalg_SOURCES = \
	arts.h \
	matpackI.cc \
	matpackI.h \
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
	make_vector.h\
	make_vector.cc\
	lin_alg.cc \
	lin_alg.h \
	test_linalg.cc\
	logic.h\
	logic.cc

test_scatproperties_SOURCES = \
	arts.h \
	matpackI.cc \
	matpackI.h \
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
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
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
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
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
	mystring.h \
	test_matpack.cc

test_xml_SOURCES = \
	absorption.cc \
	absorption.h \
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	partition_function_data.cc \
	matpackI.cc \
	matpackI.h \
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
	xml_io.h \
	xml_io_private.h \
	xml_io.cc \
	messages.cc \
	species_data.cc \
	test_xml.cc

test_zeemanproperties_SOURCES = \
	absorption.cc \
	absorption.h \
	array.h \
	arts.h \
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	auto_wsv_pointers.cc \
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	check_input.cc \
	check_input.h \
	cloudbox.h \
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
	gridded_fields.cc\
	gridded_fields.h\
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
	m_optproperties.cc \
	m_physics.cc \
	m_ppath.cc \
	m_refraction.cc \
	m_rte.cc \
	m_scatrte.cc \
	m_sensor.cc \
	m_clonesize.cc \
	m_absorption.cc \
	m_copy.h \
	m_ignore.h \
	m_xml.cc \
	m_xml.h \
	make_array.h \
	make_vector.cc \
	make_vector.h \
	math_funcs.cc \
	math_funcs.h \
	matpackI.cc \
	matpackI.h \
	matpackII.h \
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
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
	refraction.h \
	refraction.cc \
	rte.cc \
	rte.h \
	special_interp.cc \
	special_interp.h \
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
	xml_io_private.h \
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
	supergeneric.h \
        zeemanproperties.cc \
	gas_abs_lookup.h \
	gas_abs_lookup.cc \
	test_zeemanproperties.cc

test_complex_SOURCES = \
	arts.h \
	complex.h \
	complex.cc \
	test_complex.cc

make_auto_wsv_groups_h_SOURCES = \
	arts.h \
	file.cc \
	file.h \
	groups.cc \
	make_auto_wsv_groups_h.cc \
	matpackI.cc \
	matpackI.h \
	messages.cc \
	version.cc

make_auto_wsv_h_SOURCES = \
	absorption.h \
	arts.h \
	auto_wsv_groups.h \
	file.cc \
	file.h \
	groups.cc \
	make_auto_wsv_h.cc \
	matpackI.cc \
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
	matpackI.cc \
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
	make_array.h \
	make_auto_md_h.cc \
	matpackI.cc \
	matpackI.h \
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
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
	make_array.h \
	make_auto_md_cc.cc \
	make_array.h \
	matpackI.cc \
	matpackI.h \
	matpackIII.cc \
	matpackIII.h \
	matpackIV.cc \
	matpackIV.h \
	matpackV.cc \
	matpackV.h \
	matpackVI.cc \
	matpackVI.h \
	matpackVII.cc \
	matpackVII.h \
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

