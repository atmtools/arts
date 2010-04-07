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
	check-code-header-consistency.sh \
	nc_io.h \
	nc_io.cc \
	nc_io_array_types.cc \
	nc_io_compound_types.cc \
	nc_io_basic_types.cc \
	nc_io_instantiation.h \
	nc_io_types.h

bin_PROGRAMS = arts

arts_LDADD = libartscore.a libmethods.a libspecies_data.a libpartition_function_data.a libcontinua.a libmatpack.a $(libgzstream_a_ldadd) $(disort_a_ldadd) $(netcdf_ldadd)
arts_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
arts_SOURCES = auto_version.h main.cc

libartscore_a_CFLAGS = $(CFLAGS_OPT) $(WARNINGCFLAGS)
libartscore_a_CXXFLAGS = $(CXXFLAGS_OPT) $(CXXNOSTRICTALIASING) $(WARNINGCFLAGS)
libartscore_a_ldflags = 
libartscore_a_SOURCES = \
	absorption.h \
	absorption.cc \
	abs_species_tags.h \
	abs_species_tags.cc \
	arts_omp.h \
	arts_omp.cc \
	gas_abs_lookup.h \
	gas_abs_lookup.cc \
	lineshapes.cc \
	arts.cc \
	arts.h \
	auto_md.cc \
	auto_md.h \
	auto_workspace.h \
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
	interpolation.cc \
	interpolation.h \
	interpolation_poly.cc \
	interpolation_poly.h \
	iotraits.h \
	jacobian.cc \
	jacobian.h \
	legendre.cc \
	legendre.h \
	math_funcs.h \
	math_funcs.cc \
	m_abs.cc \
	m_agenda.cc \
	m_append.h \
	m_atmosphere.cc \
	m_create.cc \
	m_basic_types.cc \
	m_basic_types.h \
	m_cloudbox.cc \
	m_disort.cc \
	m_general.cc \
	m_general.h \
	m_nc.h \
	m_optproperties.cc \
	m_physics.cc \
	m_ppath.cc \
	m_refraction.cc \
	m_rte.cc \
	m_doit.cc \
	m_surface.cc \
	m_abs_lookup.cc \
	m_copy.h \
	m_delete.h \
	m_ignore.h \
	m_jacobian.cc \
	m_montecarlo.cc \
	m_select.h \
	m_sensor.cc \
	m_xml.h \
	m_batch.cc \
	m_extract.h \
	mc_antenna.cc \
	mc_antenna.h \
	mc_interp.cc \
	mc_interp.h \
	messages.cc \
	messages.h \
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
	doit.cc \
	doit.h \
	sorting.h \
	sourcetext.cc \
	sourcetext.h \
	special_interp.cc \
	special_interp.h \
	sstream.h \
	token.cc \
	token.h \
	workspace.cc \
	workspace_ng.cc \
	workspace_ng.h \
	xml_io.h \
	xml_io_instantiation.h \
	xml_io_private.h \
	xml_io_basic_types.cc \
	xml_io_compound_types.cc \
	xml_io_array_types.cc \
	xml_io_types.h \
	xml_io.cc \
	agenda_class.h \
	agenda_class.cc \
	agenda_record.h \
	agenda_record.cc \
	agendas.cc \
	lin_alg.h \
	lin_alg.cc \
	supergeneric.h

if ENABLE_NETCDF
netcdf_ldadd = -lnetcdf -lnetcdf_c++
libartscore_a_SOURCES += \
	nc_io.h \
	nc_io.cc \
	nc_io_array_types.cc \
	nc_io_compound_types.cc \
	nc_io_basic_types.cc \
	nc_io_instantiation.h \
	nc_io_types.h
endif

libcontinua_a_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS) $(CXXNOSTRICTALIASING)
libcontinua_a_SOURCES = continua.cc continua.h

libmatpack_a_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
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

libmethods_a_CXXFLAGS = $(CXXFLAGS_NOOPT) $(WARNINGCFLAGS)
libmethods_a_SOURCES = \
	array.h \
	arts.h \
	make_array.h \
	methods.h \
	methods.cc \
	wsv_aux.h \
	groups.cc

libspecies_data_a_CXXFLAGS = $(CXXFLAGS_NOOPT) $(WARNINGCFLAGS)
libspecies_data_a_SOURCES = \
	absorption.h \
	arts.h \
	make_array.h \
	species_data.cc

libpartition_function_data_a_CXXFLAGS = $(CXXFLAGS_NOOPT) $(WARNINGCFLAGS)
libpartition_function_data_a_SOURCES = \
	absorption.h \
	arts.h \
	make_array.h \
	partition_function_data.cc

noinst_LIBRARIES = \
	libartscore.a \
	libcontinua.a \
	libmatpack.a \
	libmethods.a \
	libpartition_function_data.a \
	libspecies_data.a

if ENABLE_ZLIB

noinst_LIBRARIES += libgzstream.a

libgzstream_a_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
libgzstream_a_ldadd = libgzstream.a -lz
libgzstream_a_SOURCES = gzstream.cc gzstream.h

endif


if ENABLE_DISORT

noinst_LIBRARIES += libdisort.a

disort_a_ldadd = -lf2c libdisort.a

libartscore_a_SOURCES += disort_DISORT.h

libdisort_a_CFLAGS = $(CFLAGS_OPT) $(subst -Werror,,$(WARNINGCFLAGS))
libdisort_a_SOURCES = \
	disort_D1MACH.c \
	disort_DISORT.c \
	disort_ErrPack.c \
	disort_FAKEMAIN.c \
	disort_LINPAK.c \
	disort_R1MACH.c

endif

noinst_PROGRAMS = \
	make_auto_md_cc \
	make_auto_md_h \
	make_auto_workspace_h

make_auto_workspace_h_LDADD = libmatpack.a
make_auto_workspace_h_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
make_auto_workspace_h_SOURCES = \
	arts.cc \
	arts.h \
	arts_omp.cc \
	arts_omp.h \
	file.cc \
	file.h \
	groups.cc \
	m_general.h \
	make_auto_workspace_h.cc \
	messages.h \
	messages.cc \
	parameters.cc \
	parameters.h

make_auto_md_h_LDADD = libmethods.a libmatpack.a
make_auto_md_h_CXXFLAGS = $(CXXFLAGS_OPT) $(CXXNOSTRICTALIASING) $(WARNINGCFLAGS)
make_auto_md_h_SOURCES = \
	agenda_class.h \
	agenda_record.cc \
	agenda_record.h \
	agendas.cc \
	arts.cc \
	arts.h \
	arts_omp.cc \
	arts_omp.h \
	auto_workspace.h \
	binio.cc \
	binio.h \
	file.cc \
	file.h \
	globals_2.cc \
	gridded_fields.cc \
	make_auto_md_h.cc \
	messages.cc \
	messages.h \
	methods.h \
	methods_aux.cc \
	parameters.cc \
	token.cc \
	token.h \
	workspace.cc \
	workspace_ng.cc \
	workspace_ng.h

make_auto_md_cc_LDADD = libmethods.a libmatpack.a
make_auto_md_cc_CXXFLAGS = $(CXXFLAGS_OPT) $(CXXNOSTRICTALIASING) $(WARNINGCFLAGS)
make_auto_md_cc_SOURCES = \
	agenda_class.h \
	agenda_record.cc \
	agenda_record.h \
	agendas.cc \
	arts.cc \
	arts.h \
	arts_omp.cc \
	arts_omp.h \
	auto_workspace.h \
	binio.cc \
	binio.h \
	file.cc \
	file.h \
	globals_2.cc \
	gridded_fields.cc \
	make_auto_md_cc.cc \
	messages.cc \
	messages.h \
	methods.h \
	methods_aux.cc \
	parameters.cc \
	token.cc \
	token.h \
	workspace.cc \
	workspace_ng.cc \
	workspace_ng.h

auto_workspace.h: make_auto_workspace_h$(EXEEXT)
	$(top_builddir)/src/make_auto_workspace_h$(EXEEXT)

auto_md.h: make_auto_md_h$(EXEEXT)
	$(top_builddir)/src/make_auto_md_h$(EXEEXT)

auto_md.cc: make_auto_md_cc$(EXEEXT)
	$(top_builddir)/src/make_auto_md_cc$(EXEEXT)

BUILT_SOURCES = \
	auto_workspace.h \
	auto_md.h \
	auto_md.cc

CLEANFILES = \
	auto_wsv.h \
	make_auto_wsv_h$(EXEEXT) \
	make_auto_wsv_groups_h$(EXEEXT) \
	auto_wsv.txt \
	auto_wsv_groups.h \
	auto_workspace.h \
	auto_md.h \
	auto_md.cc \
	check-code-header-consistency.log \
	test_include \
	test_include.cc

if MAINTAINER_MODE

TESTS_ENVIRONMENT = SRCDIR=$(srcdir) COMPILE="$(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) $(CXXNOSTRICTALIASING)" sh -c

TESTS = \
	check-code-cheaders.pl \
	check-code-tab-indent.sh \
	check-code-header-consistency.sh

endif

noinst_PROGRAMS += \
	test_binaryio \
	test_matpack \
	test_sparse \
	test_poly \
	test_readpp \
	test_tensor \
	test_legendre \
	test_linalg \
	test_interpolation \
	test_xml \
	test_complex \
	test_sorting \
	test_integration

if ENABLE_DISORT

noinst_PROGRAMS += test_disort

test_disort_a_ldadd = $(disort_a_ldadd)
test_disort_a_CFLAGS = $(CFLAGS_OPT) $(subst -Werror,,$(WARNINGCFLAGS))
test_disort_SOURCES = \
	disort_DISOTEST.c

endif

test_binaryio_LDADD = libartscore.a libmethods.a libspecies_data.a libpartition_function_data.a libcontinua.a libmatpack.a $(libgzstream_a_ldadd) $(disort_a_ldadd) $(netcdf_ldadd)
test_binaryio_LDFLAGS = $(libartscore_a_ldflags)
test_binaryio_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_binaryio_SOURCES = test_binaryio.cc

#test_geomag_calc_LDADD = libartscore.a libmethods.a libspecies_data.a libpartition_function_data.a libcontinua.a libmatpack.a $(libgzstream_a_ldadd) $(disort_a_ldadd) $(netcdf_ldadd)
#test_geomag_calc_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
#test_geomag_calc_SOURCES = test_geomag_calc.cc

test_interpolation_LDADD = libartscore.a libmethods.a libspecies_data.a libpartition_function_data.a libcontinua.a libmatpack.a $(libgzstream_a_ldadd) $(disort_a_ldadd) $(netcdf_ldadd)
test_interpolation_LDFLAGS = $(libartscore_a_ldflags)
test_interpolation_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_interpolation_SOURCES = test_interpolation.cc

test_legendre_LDADD = libmatpack.a
test_legendre_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_legendre_SOURCES = \
	constants.cc \
	math_funcs.h \
	math_funcs.cc \
	legendre.h \
	legendre.cc \
	test_legendre.cc

test_linalg_LDADD = libmatpack.a
test_linalg_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_linalg_SOURCES = \
	constants.cc \
	make_vector.h\
	make_vector.cc\
	lin_alg.cc \
	lin_alg.h \
	test_linalg.cc

test_poly_LDADD = libmatpack.a
test_poly_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_poly_SOURCES = \
	poly_roots.h \
	poly_roots.cc \
	test_poly.cc

test_tensor_LDADD = libmatpack.a
test_tensor_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_tensor_SOURCES = \
	test_tensor.cc

test_matpack_LDADD = libmatpack.a
test_matpack_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_matpack_SOURCES = \
	describe.h \
	describe.cc \
	test_matpack.cc

test_sparse_LDADD = libartscore.a libmethods.a libspecies_data.a libpartition_function_data.a libcontinua.a libmatpack.a $(libgzstream_a_ldadd) $(disort_a_ldadd) $(netcdf_ldadd)
test_sparse_LDFLAGS = $(libartscore_a_ldflags)
test_sparse_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_sparse_SOURCES = test_sparse.cc

test_xml_LDADD = libartscore.a libmethods.a libspecies_data.a libpartition_function_data.a libcontinua.a libmatpack.a $(libgzstream_a_ldadd) $(disort_a_ldadd) $(netcdf_ldadd)
test_xml_LDFLAGS =  $(libartscore_a_ldflags)
test_xml_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_xml_SOURCES = test_xml.cc

test_complex_LDADD = libmatpack.a
test_complex_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_complex_SOURCES = test_complex.cc

test_sorting_LDADD = libmatpack.a
test_sorting_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_sorting_SOURCES = \
	sorting.h \
	test_sorting.cc

test_integration_LDADD = libmatpack.a
test_integration_CXXFLAGS = $(CXXFLAGS_OPT) $(WARNINGCFLAGS)
test_integration_SOURCES = \
	constants.cc \
	math_funcs.h \
	math_funcs.cc \
	test_integration.cc 

test_readpp_LDADD = libmatpack.a
test_readpp_CXXFLAGS = $(CXXFLAGS_OPT) $(CXXNOSTRICTALIASING) $(WARNINGCFLAGS)
test_readpp_SOURCES = \
	bifstream.cc \
	bifstream.h \
	binio.cc \
	binio.h \
	test_readpp.cc


