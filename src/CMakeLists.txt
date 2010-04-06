configure_file (${CMAKE_CURRENT_SOURCE_DIR}/auto_version.h.in
                ${CMAKE_CURRENT_BINARY_DIR}/auto_version.h @ONLY)

include_directories (${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
                     ${ZLIB_INCLUDE_DIR})

if (NETCDF_FOUND)
  include_directories (${NETCDF_INCLUDE_DIR})
endif (NETCDF_FOUND)

set (ALL_ARTS_LIBRARIES
  artscore
  methods
  species_data
  partition_function_data
  continua
  matpack
  gzstream
  )

if (ENABLE_DISORT)
  list (APPEND ALL_ARTS_LIBRARIES disort)
endif (ENABLE_DISORT)

########### next target ###############

add_executable (arts auto_version.h main.cc)

target_link_libraries (arts ${ALL_ARTS_LIBRARIES})

install (TARGETS arts RUNTIME DESTINATION bin)

if (CMAKE_COMPILER_IS_GNUCXX)
  set_source_files_properties (continua.cc PROPERTIES
                               COMPILE_FLAGS "-fno-strict-aliasing")
  set_source_files_properties (binio.cc PROPERTIES
                               COMPILE_FLAGS "-fno-strict-aliasing")
  set_source_files_properties (methods.cc PROPERTIES COMPILE_FLAGS "-O0")
  set_source_files_properties (partition_function_data.cc PROPERTIES
                               COMPILE_FLAGS "-O0")
  set_source_files_properties (species_data.cc PROPERTIES COMPILE_FLAGS "-O0")
endif(CMAKE_COMPILER_IS_GNUCXX)

########### next target ###############

add_executable (make_auto_workspace_h
	arts.cc
	arts.h
	arts_omp.cc
	arts_omp.h
	file.cc
	file.h
	groups.cc
	m_general.h
	make_auto_workspace_h.cc
	messages.h
	messages.cc
	parameters.cc
	parameters.h
	)

target_link_libraries (make_auto_workspace_h matpack)

add_custom_command (
	OUTPUT auto_workspace.h
	COMMAND make_auto_workspace_h
	DEPENDS make_auto_workspace_h
	)

add_custom_target (BuildAutoSources DEPENDS auto_workspace.h auto_md.h
                   auto_md.cc)

########### next target ###############

add_executable (make_auto_md_h
	agenda_class.h
	agenda_record.cc
	agenda_record.h
	agendas.cc
	arts.cc
	arts.h
	arts_omp.cc
	arts_omp.h
	auto_workspace.h
	binio.cc
	binio.h
	file.cc
	file.h
	globals_2.cc
	gridded_fields.cc
	make_auto_md_h.cc
	messages.cc
	messages.h
	methods.h
	methods_aux.cc
	parameters.cc
	token.cc
	token.h
	workspace.cc
	workspace_ng.cc
	workspace_ng.h
	)

target_link_libraries (make_auto_md_h methods matpack)

add_custom_command (
	OUTPUT auto_md.h
	COMMAND make_auto_md_h
	DEPENDS make_auto_md_h
)

########### next target ###############

add_executable (make_auto_md_cc
	agenda_class.h
	agenda_record.cc
	agenda_record.h
	agendas.cc
	arts.cc
	arts.h
	arts_omp.cc
	arts_omp.h
	auto_workspace.h
	binio.cc
	binio.h
	file.cc
	file.h
	globals_2.cc
	gridded_fields.cc
	make_auto_md_cc.cc
	messages.cc
	messages.h
	methods.h
	methods_aux.cc
	parameters.cc
	token.cc
	token.h
	workspace.cc
	workspace_ng.cc
	workspace_ng.h
	)

target_link_libraries (make_auto_md_cc methods matpack)

add_custom_command (
	OUTPUT auto_md.cc
	COMMAND make_auto_md_cc
	DEPENDS make_auto_md_cc
	)

########### next target ###############

add_library (artscore STATIC
	absorption.h
	absorption.cc
	abs_species_tags.h
	abs_species_tags.cc
	arts_omp.h
	arts_omp.cc
	gas_abs_lookup.h
	gas_abs_lookup.cc
	lineshapes.cc
	arts.cc
	arts.h
	auto_md.cc
	auto_md.h
	auto_workspace.h
	bifstream.h
	bifstream.cc
	bofstream.h
	bofstream.cc
	binio.h
	binio.cc
	check_input.cc
	check_input.h
	cloudbox.h
	cloudbox.cc
	constants.cc
	disort.cc
	disort.h
	fastem.cc
	fastem.h
	file.cc
	file.h
	geomag_calc.cc
	geomag_calc.h
	getopt.c
	getopt.h
	getopt1.c
	globals_2.cc
	gridded_fields.cc
	gridded_fields.h
	interpolation.cc
	interpolation.h
	interpolation_poly.cc
	interpolation_poly.h
	iotraits.h
	jacobian.cc
	jacobian.h
	legendre.cc
	legendre.h
	math_funcs.h
	math_funcs.cc
	m_abs.cc
	m_agenda.cc
	m_append.h
	m_atmosphere.cc
	m_create.cc
	m_basic_types.cc
	m_basic_types.h
	m_cloudbox.cc
	m_disort.cc
	m_general.cc
	m_general.h
	m_nc.h
	m_optproperties.cc
	m_physics.cc
	m_ppath.cc
	m_refraction.cc
	m_rte.cc
	m_scatrte.cc
	m_surface.cc
	m_abs_lookup.cc
	m_copy.h
	m_delete.h
	m_ignore.h
	m_jacobian.cc
	m_montecarlo.cc
	m_select.h
	m_sensor.cc
	m_xml.h
	m_batch.cc
	m_extract.h
	mc_antenna.cc
	mc_antenna.h
	mc_interp.cc
	mc_interp.h
	messages.cc
	messages.h
	methods.h
	methods_aux.cc
	montecarlo.cc
	montecarlo.h
  nc_io.h
  nc_io.cc
  nc_io_array_types.cc
	nc_io_compound_types.cc
	nc_io_basic_types.cc
	nc_io_instantiation.h
	nc_io_types.h
	optproperties.cc
	optproperties.h
	parameters.cc
	parameters.h
	parser.cc
	parser.h
	ppath.h
	ppath.cc
	poly_roots.h
	poly_roots.cc
	physics_funcs.cc
	physics_funcs.h
	refraction.h
	refraction.cc
	rng.cc
	rng.h
	rte.cc
	rte.h
	sensor.cc
	sensor.h
	scatrte.cc
	scatrte.h
	sorting.h
	sourcetext.cc
	sourcetext.h
	special_interp.cc
	special_interp.h
	sstream.h
	token.cc
	token.h
	workspace.cc
	workspace_ng.cc
	workspace_ng.h
	xml_io.h
	xml_io_instantiation.h
	xml_io_private.h
	xml_io_basic_types.cc
	xml_io_compound_types.cc
	xml_io_array_types.cc
	xml_io_types.h
	xml_io.cc
	agenda_class.h
	agenda_class.cc
	agenda_record.h
	agenda_record.cc
	agendas.cc
	lin_alg.h
	lin_alg.cc
	supergeneric.h
	)

if (NETCDF_FOUND)
  target_link_libraries (artscore ${NETCDF_LIBRARIES})
endif (NETCDF_FOUND)

########### next target ###############

add_library (continua STATIC continua.cc continua.h)

########### next target ###############

add_library (matpack STATIC
	array.h
	complex.h
	complex.cc
	exceptions.h
	logic.h
	logic.cc
	make_array.h
	make_vector.h
	make_vector.cc
	matpack.h
	matpackI.h
	matpackI.cc
	matpackII.h
	matpackII.cc
	matpackIII.h
	matpackIII.cc
	matpackIV.h
	matpackIV.cc
	matpackV.h
	matpackV.cc
	matpackVI.h
	matpackVI.cc
	matpackVII.h
	matpackVII.cc
	mystring.h
	)

########### next target ###############

add_library (methods STATIC
	array.h
	arts.h
	make_array.h
	methods.h
	methods.cc
	wsv_aux.h
	groups.cc
	)

########### next target ###############

add_library (species_data STATIC
	absorption.h
	arts.h
	make_array.h
	species_data.cc
	)

########### next target ###############

add_library (partition_function_data STATIC
	absorption.h
	arts.h
	make_array.h
	partition_function_data.cc
	)

########### next target ###############

add_library (gzstream STATIC gzstream.cc gzstream.h)
target_link_libraries (gzstream ${ZLIB_LIBRARIES})

########### next target ###############

if (ENABLE_DISORT)
  add_library (disort STATIC
    disort_D1MACH.c
    disort_DISORT.c
    disort_ErrPack.c
    disort_FAKEMAIN.c
    disort_LINPAK.c
    disort_R1MACH.c
    )
  target_link_libraries (disort f2c)

  ########### testcase ###############

  add_executable (test_disort disort_DISOTEST.c)

  target_link_libraries (test_disort disort)
endif (ENABLE_DISORT)

########### testcases ###############

add_executable (test_binaryio test_binaryio.cc)

target_link_libraries (test_binaryio ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_interpolation test_interpolation.cc)
target_link_libraries (test_interpolation ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_legendre
  constants.cc
  math_funcs.h
  math_funcs.cc
  legendre.h
  legendre.cc
  test_legendre.cc)

target_link_libraries (test_legendre matpack)

########### next testcase ###############

add_executable (test_linalg
  constants.cc
  make_vector.h
  make_vector.cc
  lin_alg.cc
  lin_alg.h
  test_linalg.cc)

target_link_libraries (test_linalg ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_integration
  constants.cc
  math_funcs.h
  math_funcs.cc
  test_integration.cc)

target_link_libraries (test_integration matpack)

########### next testcase ###############

add_executable (test_poly poly_roots.h poly_roots.cc test_poly.cc)
target_link_libraries (test_poly matpack)

########### next testcase ###############

add_executable (test_tensor test_tensor.cc)
target_link_libraries (test_tensor matpack)

########### next testcase ###############

add_executable (test_matpack describe.h describe.cc test_matpack.cc)
target_link_libraries (test_matpack matpack)

########### next testcase ###############

add_executable (test_sparse test_sparse.cc)
target_link_libraries (test_sparse ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_xml test_xml.cc)
target_link_libraries (test_xml ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_complex test_complex.cc)
target_link_libraries (test_complex matpack)

########### next testcase ###############

add_executable (test_sorting sorting.h test_sorting.cc)
target_link_libraries (test_sorting matpack)

########### next testcase ###############

add_executable (test_readpp test_readpp.cc)
target_link_libraries (test_readpp ${ALL_ARTS_LIBRARIES})

