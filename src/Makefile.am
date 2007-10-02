# Copyright (C) 2000-2007 Stefan Buehler <sbuehler@ltu.se>
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

EXTRA_DIST = \
	check-code-cheaders.pl \
	check-code-tab-indent.sh \
	check-code-header-consistency.sh

bin_PROGRAMS = arts

arts_LDADD = libmatpack.a $(disort_ldflags)
arts_SOURCES = \
	main.cc \
	absorption.h \
	absorption.cc \
	abs_species_tags.h \
	abs_species_tags.cc \
	gas_abs_lookup.h \
	gas_abs_lookup.cc \
	lineshapes.cc \
	arts.cc \
	arts.h \
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	binio.h \
	binio.cc \
	check_input.cc \
	check_input.h \
	cloudbox.h \
	cloudbox.cc \
	constants.cc \
	continua.cc \
	continua.h \
	disort.cc \
	disort.h \
	fastem.cc \
	fastem.h \
	file.cc \
	file.h \
	geomag_calc.cc \
	geomag_calc.h \
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
	jacobian.cc \
	jacobian.h \
	legendre.cc \
	legendre.h \
	m_abs.cc \
	m_agenda.cc \
        m_atmosphere.cc \
	m_basic_types.cc \
	m_basic_types.h \
	m_cloudbox.cc \
	m_disort.cc \
	m_general.cc \
	m_general.h \
	m_optproperties.cc \
	m_physics.cc \
	m_ppath.cc \
	m_refraction.cc \
	m_rte.cc \
	m_scatrte.cc \
	m_abs_lookup.cc \
	m_copy.h \
	m_ignore.h \
	m_jacobian.cc \
	m_montecarlo.cc \
	m_sensor.cc \
	m_xml.cc \
	m_xml.h \
	m_batch.cc \
	mc_antenna.cc \
	mc_antenna.h \
	mc_interp.cc \
	mc_interp.h \
	messages.cc \
	messages.h \
	methods.cc \
	methods.h \
	methods_aux.cc \
	montecarlo.cc \
	montecarlo.h \
	optproperties.cc \
	optproperties.h \
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
	rng.cc \
	rng.h \
	rte.cc \
	rte.h \
	sensor.cc \
	sensor.h \
	scatrte.cc \
	scatrte.h \
	sorting.h \
	special_interp.cc \
	special_interp.h \
	species_data.cc \
	sstream.h \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	workspace_ng.cc \
	workspace_ng.h \
	wsv_aux.h \
	xml_io.h \
	xml_io_instantiation.h \
	xml_io_private.h \
	xml_io_basic_types.cc \
	xml_io_basic_types.h \
	xml_io_compound_types.cc \
	xml_io_compound_types.h \
	xml_io_array_types.cc \
	xml_io_array_types.h \
	xml_io.cc \
	agenda_class.h \
	agenda_class.cc \
	agenda_record.h \
	agenda_record.cc \
	agendas.cc \
	lin_alg.h \
	lin_alg.cc \
	supergeneric.h


libmatpack_a_SOURCES = \
	array.h \
	complex.h \
	complex.cc \
	exceptions.h \
	logic.h \
	logic.cc \
	make_array.h \
	make_vector.h \
	make_vector.cc \
	math_funcs.h \
	math_funcs.cc \
	matpack.h \
	matpackI.h \
	matpackI.cc \
	matpackII.h \
	matpackII.cc \
	matpackIII.h \
	matpackIII.cc \
	matpackIV.h \
	matpackIV.cc \
	matpackV.h \
	matpackV.cc \
	matpackVI.h \
	matpackVI.cc \
	matpackVII.h \
	matpackVII.cc \
	mystring.h

if ENABLE_DISORT

noinst_LIBRARIES = libmatpack.a libdisort.a

disort_ldflags = -lf2c libdisort.a -lf2c

libdisort_a_SOURCES = \
	disort_D1MACH.c \
	disort_DISORT.c \
	disort_ErrPack.c \
	disort_FAKEMAIN.c \
	disort_LINPAK.c \
	disort_R1MACH.c

test_disort_LDADD = $(disort_ldflags)
test_disort_SOURCES = \
	disort_DISOTEST.c

else

noinst_LIBRARIES = libmatpack.a

disort_ldflags = 

test_disort$(EXEEXT):

endif

continua.o: continua.cc arts.h matpackI.h array.h absorption.h messages.h continua.h
	if $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
		$(AM_CXXFLAGS) $(CXXFLAGS) $(CXXNOSTRICTALIASING) -MT $@ \
		-MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ $<; \
		then mv -f "$(DEPDIR)/$*.Tpo" "$(DEPDIR)/$*.Po"; \
		else rm -f "$(DEPDIR)/$*.Tpo"; exit 1; fi

methods.o: methods.cc arts.h make_array.h auto_wsv.h methods.h auto_wsv_groups.h
	if $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
		$(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -MT $@ \
		-MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ $<; \
		then mv -f "$(DEPDIR)/$*.Tpo" "$(DEPDIR)/$*.Po"; \
		else rm -f "$(DEPDIR)/$*.Tpo"; exit 1; fi

species_data.o: species_data.cc arts.h make_array.h absorption.h
	if $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
		$(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -MT $@ \
		-MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ $<; \
		then mv -f "$(DEPDIR)/$*.Tpo" "$(DEPDIR)/$*.Po"; \
		else rm -f "$(DEPDIR)/$*.Tpo"; exit 1; fi

partition_function_data.o: partition_function_data.cc arts.h make_array.h absorption.h
	if $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
		$(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -MT $@ \
		-MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ $<; \
		then mv -f "$(DEPDIR)/$*.Tpo" "$(DEPDIR)/$*.Po"; \
		else rm -f "$(DEPDIR)/$*.Tpo"; exit 1; fi

# OLE: Temporarely added to work around a bug in gcc 3.0
xml_io.o: xml_io.cc arts.h xml_io.h xml_io_instantiation.h xml_io_private.h xml_io_basic_types.h xml_io_compound_types.h xml_io_array_types.h bofstream.h bifstream.h
	if $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
		$(AM_CXXFLAGS) $(CXXFLAGS_FIXED) -MT $@ \
		-MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ $<; \
		then mv -f "$(DEPDIR)/$*.Tpo" "$(DEPDIR)/$*.Po"; \
		else rm -f "$(DEPDIR)/$*.Tpo"; exit 1; fi

noinst_PROGRAMS = \
	make_auto_md_cc \
	make_auto_md_h \
	make_auto_wsv_groups_h \
	make_auto_wsv_h \
	test_binaryio \
	test_disort \
	test_matpack \
	test_sparse \
	test_poly \
	test_readpp \
	test_tensor \
	test_legendre \
	test_linalg \
	test_interpolation \
	test_sensor \
	test_xml \
	test_complex \
	test_sorting \
	test_integration

test_binaryio_LDADD = libmatpack.a $(disort_ldflags)
test_binaryio_SOURCES = \
	absorption.h \
	absorption.cc \
	abs_species_tags.h \
	abs_species_tags.cc \
	gas_abs_lookup.h \
	gas_abs_lookup.cc \
	lineshapes.cc \
	arts.h \
	arts.cc \
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	binio.h \
	binio.cc \
	check_input.cc \
	check_input.h \
	cloudbox.h \
	cloudbox.cc \
	constants.cc \
	continua.cc \
	continua.h \
	disort.cc \
	disort.h \
	disort_DISORT.h \
	fastem.cc \
	fastem.h \
	file.cc \
	file.h \
	geomag_calc.cc \
	geomag_calc.h \
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
	jacobian.cc \
	jacobian.h \
	legendre.cc \
	legendre.h \
	m_abs.cc \
	m_agenda.cc \
        m_atmosphere.cc \
	m_basic_types.cc \
	m_basic_types.h \
	m_cloudbox.cc \
	m_disort.cc \
	m_general.cc \
	m_general.h \
	m_optproperties.cc \
	m_physics.cc \
	m_ppath.cc \
	m_refraction.cc \
	m_rte.cc \
	m_scatrte.cc \
	m_abs_lookup.cc \
	m_copy.h \
	m_ignore.h \
	m_jacobian.cc \
	m_montecarlo.cc \
	m_sensor.cc \
	m_xml.cc \
	m_xml.h \
	m_batch.cc \
	mc_antenna.cc \
	mc_antenna.h \
	mc_interp.cc \
	mc_interp.h \
	messages.cc \
	messages.h \
	methods.cc \
	methods.h \
	methods_aux.cc \
	montecarlo.cc \
	montecarlo.h \
	optproperties.cc \
	optproperties.h \
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
	rng.cc \
	rng.h \
	rte.cc \
	rte.h \
	sensor.cc \
	sensor.h \
	scatrte.cc \
	scatrte.h \
	sorting.h \
	special_interp.cc \
	special_interp.h \
	species_data.cc \
	sstream.h \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	workspace_ng.cc \
	workspace_ng.h \
	wsv_aux.h \
	xml_io.h \
	xml_io_instantiation.h \
	xml_io_private.h \
	xml_io_basic_types.cc \
	xml_io_basic_types.h \
	xml_io_compound_types.cc \
	xml_io_compound_types.h \
	xml_io_array_types.cc \
	xml_io_array_types.h \
	xml_io.cc \
	agenda_class.h \
	agenda_class.cc \
	agenda_record.h \
	agenda_record.cc \
	agendas.cc \
	lin_alg.h\
	lin_alg.cc\
	supergeneric.h \
	test_binaryio.cc

#test_geomag_calc_SOURCES = \
#	absorption.cc \
#	absorption.h \
#	arts.cc \
#	arts.h \
#	auto_md.cc \
#	auto_md.h \
#	auto_wsv.h \
#	auto_wsv_groups.h \
#	bifstream.h \
#	bifstream.cc \
#	bofstream.h \
#	bofstream.cc \
#	binio.h \
#	binio.cc \
#	check_input.cc \
#	check_input.h \
#	cloudbox.h \
#	cloudbox.cc \
#	constants.cc \
#	fastem.cc \
#	fastem.h \
#	file.cc \
#	file.h \
#	getopt.c \
#	getopt.h \
#	getopt1.c \
#	globals_2.cc \
#	gridded_fields.cc\
#	gridded_fields.h\
#	groups.cc \
#	interpolation.cc \
#	interpolation.h \
#	iotraits.h \
#	jacobian.cc \
#	jacobian.h \
#	legendre.cc \
#	legendre.h \
#	m_agenda.cc \
#	m_atmosphere.cc \
#	m_basic_types.cc \
#	m_batch.cc \
#	m_cloudbox.cc \
#	m_general.cc \
#	m_general.h \
#	m_jacobian.cc \
#	m_montecarlo.cc \
#	m_optproperties.cc \
#	m_physics.cc \
#	m_ppath.cc \
#	m_refraction.cc \
#	m_rte.cc \
#	m_scatrte.cc \
#	m_sensor.cc \
#	m_abs_lookup.cc \
#	m_copy.h \
#	m_ignore.h \
#	m_xml.cc \
#	m_xml.h \
#	mc_interp.cc \
#	mc_interp.h \
#	messages.cc \
#	messages.h \
#	methods.cc \
#	methods.h \
#	methods_aux.cc \
#	montecarlo.cc \
#	montecarlo.h \
#	optproperties.cc \
#	optproperties.h \
#	parameters.cc \
#	parameters.h \
#	parser.cc \
#	parser.h \
#	partition_function_data.cc \
#	ppath.h \
#	ppath.cc \
#	poly_roots.h \
#	poly_roots.cc \
#	physics_funcs.cc \
#	physics_funcs.h \
#	refraction.h \
#	refraction.cc \
#	rng.cc \
#	rng.h \
#	rte.cc \
#	rte.h \
#	sensor.cc \
#	sensor.h \
#	scatrte.cc \
#	scatrte.h \
#	special_interp.cc \
#	special_interp.h \
#	species_data.cc \
#	sstream.h \
#	token.cc \
#	token.h \
#	version.cc \
#	workspace.cc \
#	workspace_aux.cc \
#	workspace_ng.cc \
#	workspace_ng.h \
#	wsv_aux.h \
#	xml_io.h \
#	xml_io_instantiation.h \
#	xml_io_private.h \
#	xml_io_basic_types.cc \
#	xml_io_basic_types.h \
#	xml_io_compound_types.cc \
#	xml_io_compound_types.h \
#	xml_io_array_types.cc \
#	xml_io_array_types.h \
#	xml_io.cc \
#	agenda_class.h \
#	agenda_class.cc \
#	agenda_record.h \
#	agenda_record.cc \
#	agendas.cc \
#	lin_alg.h\
#	lin_alg.cc\
#	supergeneric.h \
#	geomag_calc.cc \
#        geomag_calc.h \
#	gas_abs_lookup.h \
#	gas_abs_lookup.cc \
#	test_geomag_calc.cc

test_interpolation_LDADD = libmatpack.a $(disort_ldflags)
test_interpolation_SOURCES = \
	absorption.h \
	absorption.cc \
	abs_species_tags.h \
	abs_species_tags.cc \
	gas_abs_lookup.h \
	gas_abs_lookup.cc \
	lineshapes.cc \
	arts.h \
	arts.cc \
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	binio.h \
	binio.cc \
	check_input.cc \
	check_input.h \
	cloudbox.h \
	cloudbox.cc \
	constants.cc \
	continua.cc \
	continua.h \
	disort.cc \
	disort.h \
	disort_DISORT.h \
	fastem.cc \
	fastem.h \
	file.cc \
	file.h \
	geomag_calc.cc \
	geomag_calc.h \
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
	jacobian.cc \
	jacobian.h \
	legendre.cc \
	legendre.h \
	m_abs.cc \
	m_agenda.cc \
        m_atmosphere.cc \
	m_basic_types.cc \
	m_basic_types.h \
	m_cloudbox.cc \
	m_disort.cc \
	m_general.cc \
	m_general.h \
	m_optproperties.cc \
	m_physics.cc \
	m_ppath.cc \
	m_refraction.cc \
	m_rte.cc \
	m_scatrte.cc \
	m_abs_lookup.cc \
	m_copy.h \
	m_ignore.h \
	m_jacobian.cc \
	m_montecarlo.cc \
	m_sensor.cc \
	m_xml.cc \
	m_xml.h \
	m_batch.cc \
	mc_antenna.cc \
	mc_antenna.h \
	mc_interp.cc \
	mc_interp.h \
	messages.cc \
	messages.h \
	methods.cc \
	methods.h \
	methods_aux.cc \
	montecarlo.cc \
	montecarlo.h \
	optproperties.cc \
	optproperties.h \
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
	rng.cc \
	rng.h \
	rte.cc \
	rte.h \
	sensor.cc \
	sensor.h \
	scatrte.cc \
	scatrte.h \
	sorting.h \
	special_interp.cc \
	special_interp.h \
	species_data.cc \
	sstream.h \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	workspace_ng.cc \
	workspace_ng.h \
	wsv_aux.h \
	xml_io.h \
	xml_io_instantiation.h \
	xml_io_private.h \
	xml_io_basic_types.cc \
	xml_io_basic_types.h \
	xml_io_compound_types.cc \
	xml_io_compound_types.h \
	xml_io_array_types.cc \
	xml_io_array_types.h \
	xml_io.cc \
	agenda_class.h \
	agenda_class.cc \
	agenda_record.h \
	agenda_record.cc \
	agendas.cc \
	lin_alg.h\
	lin_alg.cc\
	supergeneric.h \
	test_interpolation.cc

test_legendre_LDADD = libmatpack.a
test_legendre_SOURCES = \
	constants.cc \
	legendre.h \
	legendre.cc \
	test_legendre.cc

test_linalg_LDADD = libmatpack.a
test_linalg_SOURCES = \
	constants.cc \
	make_vector.h\
	make_vector.cc\
	lin_alg.cc \
	lin_alg.h \
	test_linalg.cc

test_poly_LDADD = libmatpack.a
test_poly_SOURCES = \
	poly_roots.h \
	poly_roots.cc \
	test_poly.cc

test_tensor_LDADD = libmatpack.a
test_tensor_SOURCES = \
	test_tensor.cc

test_matpack_LDADD = libmatpack.a
test_matpack_SOURCES = \
	describe.h \
	describe.cc \
	test_matpack.cc

test_sparse_LDADD = libmatpack.a $(disort_ldflags)
test_sparse_SOURCES = \
	absorption.h \
	absorption.cc \
	abs_species_tags.h \
	abs_species_tags.cc \
	gas_abs_lookup.h \
	gas_abs_lookup.cc \
	lineshapes.cc \
	arts.cc \
	arts.h \
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	binio.h \
	binio.cc \
	check_input.cc \
	check_input.h \
	cloudbox.h \
	cloudbox.cc \
	constants.cc \
	continua.cc \
	continua.h \
	disort.cc \
	disort.h \
	disort_DISORT.h \
	fastem.cc \
	fastem.h \
	file.cc \
	file.h \
	geomag_calc.cc \
	geomag_calc.h \
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
	jacobian.cc \
	jacobian.h \
	legendre.cc \
	legendre.h \
	m_abs.cc \
	m_agenda.cc \
        m_atmosphere.cc \
	m_basic_types.cc \
	m_basic_types.h \
	m_cloudbox.cc \
	m_disort.cc \
	m_general.cc \
	m_general.h \
	m_optproperties.cc \
	m_physics.cc \
	m_ppath.cc \
	m_refraction.cc \
	m_rte.cc \
	m_scatrte.cc \
	m_abs_lookup.cc \
	m_copy.h \
	m_ignore.h \
	m_jacobian.cc \
	m_montecarlo.cc \
	m_sensor.cc \
	m_xml.cc \
	m_xml.h \
	m_batch.cc \
	mc_interp.cc \
	mc_interp.h \
	mc_antenna.cc \
	mc_antenna.h \
	messages.cc \
	messages.h \
	methods.cc \
	methods.h \
	methods_aux.cc \
	montecarlo.cc \
	montecarlo.h \
	optproperties.cc \
	optproperties.h \
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
	rng.cc \
	rng.h \
	rte.cc \
	rte.h \
	sensor.cc \
	sensor.h \
	scatrte.cc \
	scatrte.h \
	sorting.h \
	special_interp.cc \
	special_interp.h \
	species_data.cc \
	sstream.h \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	workspace_ng.cc \
	workspace_ng.h \
	wsv_aux.h \
	xml_io.h \
	xml_io_instantiation.h \
	xml_io_private.h \
	xml_io_basic_types.cc \
	xml_io_basic_types.h \
	xml_io_compound_types.cc \
	xml_io_compound_types.h \
	xml_io_array_types.cc \
	xml_io_array_types.h \
	xml_io.cc \
	agenda_class.h \
	agenda_class.cc \
	agenda_record.h \
	agenda_record.cc \
	agendas.cc \
	lin_alg.h\
	lin_alg.cc\
	supergeneric.h \
	test_sparse.cc

test_xml_LDADD = libmatpack.a $(disort_ldflags)
test_xml_SOURCES = \
	absorption.h \
	absorption.cc \
	abs_species_tags.h \
	abs_species_tags.cc \
	gas_abs_lookup.h \
	gas_abs_lookup.cc \
	lineshapes.cc \
	arts.cc \
	arts.h \
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	binio.h \
	binio.cc \
	check_input.cc \
	check_input.h \
	cloudbox.h \
	cloudbox.cc \
	constants.cc \
	continua.cc \
	continua.h \
	disort.cc \
	disort.h \
	disort_DISORT.h \
	fastem.cc \
	fastem.h \
	file.cc \
	file.h \
	geomag_calc.cc \
	geomag_calc.h \
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
	jacobian.cc \
	jacobian.h \
	legendre.cc \
	legendre.h \
	m_abs.cc \
	m_agenda.cc \
        m_atmosphere.cc \
	m_basic_types.cc \
	m_basic_types.h \
	m_cloudbox.cc \
	m_disort.cc \
	m_general.cc \
	m_general.h \
	m_optproperties.cc \
	m_physics.cc \
	m_ppath.cc \
	m_refraction.cc \
	m_rte.cc \
	m_scatrte.cc \
	m_abs_lookup.cc \
	m_copy.h \
	m_ignore.h \
	m_jacobian.cc \
	m_montecarlo.cc \
	m_sensor.cc \
	m_xml.cc \
	m_xml.h \
	m_batch.cc \
	mc_antenna.cc \
	mc_antenna.h \
	mc_interp.cc \
	mc_interp.h \
	messages.cc \
	messages.h \
	methods.cc \
	methods.h \
	methods_aux.cc \
	montecarlo.cc \
	montecarlo.h \
	optproperties.cc \
	optproperties.h \
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
	rng.cc \
	rng.h \
	rte.cc \
	rte.h \
	sensor.cc \
	sensor.h \
	scatrte.cc \
	scatrte.h \
	sorting.h \
	special_interp.cc \
	special_interp.h \
	species_data.cc \
	sstream.h \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	workspace_ng.cc \
	workspace_ng.h \
	wsv_aux.h \
	xml_io.h \
	xml_io_instantiation.h \
	xml_io_private.h \
	xml_io_basic_types.cc \
	xml_io_basic_types.h \
	xml_io_compound_types.cc \
	xml_io_compound_types.h \
	xml_io_array_types.cc \
	xml_io_array_types.h \
	xml_io.cc \
	agenda_class.h \
	agenda_class.cc \
	agenda_record.h \
	agenda_record.cc \
	agendas.cc \
	lin_alg.h \
	lin_alg.cc \
	supergeneric.h \
	test_xml.cc

test_complex_LDADD = libmatpack.a
test_complex_SOURCES = \
	test_complex.cc

test_sorting_LDADD = libmatpack.a
test_sorting_SOURCES = \
	sorting.h \
	test_sorting.cc

test_sensor_LDADD = libmatpack.a $(disort_ldflags)
test_sensor_SOURCES = \
	absorption.h \
	absorption.cc \
	abs_species_tags.h \
	abs_species_tags.cc \
	gas_abs_lookup.h \
	gas_abs_lookup.cc \
	lineshapes.cc \
	arts.cc \
	arts.h \
	auto_md.cc \
	auto_md.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	bifstream.h \
	bifstream.cc \
	bofstream.h \
	bofstream.cc \
	binio.h \
	binio.cc \
	check_input.cc \
	check_input.h \
	cloudbox.h \
	cloudbox.cc \
	constants.cc \
	continua.cc \
	continua.h \
	disort.cc \
	disort.h \
	disort_DISORT.h \
	fastem.cc \
	fastem.h \
	file.cc \
	file.h \
	geomag_calc.cc \
	geomag_calc.h \
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
	jacobian.cc \
	jacobian.h \
	legendre.cc \
	legendre.h \
	m_abs.cc \
	m_agenda.cc \
        m_atmosphere.cc \
	m_basic_types.cc \
	m_basic_types.h \
	m_cloudbox.cc \
	m_disort.cc \
	m_general.cc \
	m_general.h \
	m_optproperties.cc \
	m_physics.cc \
	m_ppath.cc \
	m_refraction.cc \
	m_rte.cc \
	m_scatrte.cc \
	m_abs_lookup.cc \
	m_copy.h \
	m_ignore.h \
	m_jacobian.cc \
	m_montecarlo.cc \
	m_sensor.cc \
	m_xml.cc \
	m_xml.h \
	m_batch.cc \
	mc_antenna.cc \
	mc_antenna.h \
	mc_interp.cc \
	mc_interp.h \
	messages.cc \
	messages.h \
	methods.cc \
	methods.h \
	methods_aux.cc \
	montecarlo.cc \
	montecarlo.h \
	optproperties.cc \
	optproperties.h \
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
	rng.cc \
	rng.h \
	rte.cc \
	rte.h \
	sensor.cc \
	sensor.h \
	scatrte.cc \
	scatrte.h \
	sorting.h \
	special_interp.cc \
	special_interp.h \
	species_data.cc \
	sstream.h \
	token.cc \
	token.h \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	workspace_ng.cc \
	workspace_ng.h \
	wsv_aux.h \
	xml_io.h \
	xml_io_instantiation.h \
	xml_io_private.h \
	xml_io_basic_types.cc \
	xml_io_basic_types.h \
	xml_io_compound_types.cc \
	xml_io_compound_types.h \
	xml_io_array_types.cc \
	xml_io_array_types.h \
	xml_io.cc \
	agenda_class.h \
	agenda_class.cc \
	agenda_record.h \
	agenda_record.cc \
	agendas.cc \
	lin_alg.h \
	lin_alg.cc \
	supergeneric.h \
	test_sensor.cc

test_integration_LDADD = libmatpack.a
test_integration_SOURCES = \
	constants.cc \
	test_integration.cc 

test_readpp_LDADD = libmatpack.a
test_readpp_SOURCES = \
	bifstream.cc \
	bifstream.h \
	binio.cc \
	binio.h \
	test_readpp.cc

make_auto_wsv_groups_h_LDADD = libmatpack.a
make_auto_wsv_groups_h_SOURCES = \
	arts.cc \
	arts.h \
	file.cc \
	file.h \
	groups.cc \
	m_general.h \
	make_auto_wsv_groups_h.cc \
	messages.cc \
	version.cc

make_auto_wsv_h_LDADD = libmatpack.a
make_auto_wsv_h_SOURCES = \
	agenda_class.h \
	arts.cc \
	arts.h \
	auto_wsv_groups.h \
	binio.cc \
	binio.h \
	file.cc \
	file.h \
	groups.cc \
	m_general.h \
	make_auto_wsv_h.cc \
	messages.cc \
	version.cc \
	workspace.cc \
	workspace_aux.cc \
	wsv_aux.h

make_auto_md_h_LDADD = libmatpack.a
make_auto_md_h_SOURCES = \
	agenda_record.cc \
	agenda_record.h \
	agendas.cc \
	arts.cc \
	arts.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	binio.cc \
	binio.h \
	file.cc \
	file.h \
	globals_2.cc \
	groups.cc \
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
	workspace_ng.cc \
	workspace_ng.h \
	wsv_aux.h

make_auto_md_cc_LDADD = libmatpack.a
make_auto_md_cc_SOURCES = \
	agenda_record.cc \
	agenda_record.h \
	agendas.cc \
	arts.cc \
	arts.h \
	auto_wsv.h \
	auto_wsv_groups.h \
	binio.cc \
	binio.h \
	file.cc \
	file.h \
	globals_2.cc \
	groups.cc \
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
	workspace_ng.cc \
	workspace_ng.h \
	wsv_aux.h

auto_wsv.h: make_auto_wsv_h$(EXEEXT)
	$(top_builddir)/src/make_auto_wsv_h$(EXEEXT)

auto_wsv_groups.h: make_auto_wsv_groups_h$(EXEEXT)
	$(top_builddir)/src/make_auto_wsv_groups_h$(EXEEXT)

auto_md.h: make_auto_md_h$(EXEEXT)
	$(top_builddir)/src/make_auto_md_h$(EXEEXT)

auto_md.cc: make_auto_md_cc$(EXEEXT)
	$(top_builddir)/src/make_auto_md_cc$(EXEEXT)

BUILT_SOURCES = \
	auto_wsv_groups.h \
	auto_wsv.h \
	auto_md.h \
	auto_md.cc

CLEANFILES = \
	auto_wsv.h \
	auto_wsv_groups.h \
	auto_md.h \
	auto_md.cc \
	check-code-header-consistency.log \
	test_include \
	test_include.cc

TESTS_ENVIRONMENT = SRCDIR=$(srcdir) COMPILE="$(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) $(CXXNOSTRICTALIASING)" sh -c

TESTS = \
	check-code-cheaders.pl \
	check-code-tab-indent.sh \
	check-code-header-consistency.sh

