
include_directories (${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
                     ${ZLIB_INCLUDE_DIR}
                     ${BLAS_INCLUDE_DIR}
                     /usr/local/MATLAB/R2015a/extern/include)




if (NETCDF_FOUND)
  include_directories (${NETCDF_INCLUDE_DIR})
endif (NETCDF_FOUND)

set (ALL_ARTS_LIBRARIES
  artscore
  methods
  matpack
  )

file (GLOB HEADERFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
add_custom_target(UtilityHeadersArts SOURCES ${HEADERFILES})

if (ENABLE_GUI)
  include (${QT_USE_FILE})
  add_definitions (${QT_DEFINITIONS})

  list (APPEND ALL_ARTS_LIBRARIES artsgui ${QT_LIBRARIES})
  include_directories ("gui")
endif (ENABLE_GUI)

if (ENABLE_DISORT)
  list (APPEND ALL_ARTS_LIBRARIES disort)
endif (ENABLE_DISORT)

if (ENABLE_FASTEM)
  list (APPEND ALL_ARTS_LIBRARIES fastem)
endif (ENABLE_FASTEM)

if (ENABLE_TMATRIX)
  list (APPEND ALL_ARTS_LIBRARIES tmatrix)
endif (ENABLE_TMATRIX)

if (ENABLE_REFICE)
  list (APPEND ALL_ARTS_LIBRARIES refice)
endif (ENABLE_REFICE)


########### target to extract version number from ChangeLog ###############

add_custom_target (auto_version_h DEPENDS auto_version.h)
add_dependencies (auto_version_h UpdateAutoVersion)

add_custom_command (
  OUTPUT auto_version.h
  COMMAND ${CMAKE_COMMAND} -D "ARTS_BINARY_DIR:STRING=\"${CMAKE_BINARY_DIR}\"" -D "ARTS_SOURCE_DIR:STRING=\"${CMAKE_SOURCE_DIR}\"" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/update_auto_version_h.cmake
  DEPENDS ${CMAKE_BINARY_DIR}/auto_version.txt
)

########### next target ###############

add_executable (arts main.cc)
add_dependencies (arts auto_version_h)

target_link_libraries (arts ${ALL_ARTS_LIBRARIES})

install (TARGETS arts RUNTIME DESTINATION bin)

set_source_files_properties (continua.cc PROPERTIES
                             COMPILE_FLAGS "-fno-strict-aliasing")
set_source_files_properties (binio.cc PROPERTIES
                             COMPILE_FLAGS "-fno-strict-aliasing")
set_source_files_properties (methods.cc PROPERTIES COMPILE_FLAGS "-O0")
set_source_files_properties (partition_function_data.cc PROPERTIES
                             COMPILE_FLAGS "-O0")
set_source_files_properties (species_data.cc PROPERTIES COMPILE_FLAGS "-O0")

########### next target ###############

add_executable (make_auto_workspace_h
	arts.cc
	file.cc
  make_auto_workspace_h.cc
	messages.cc
	parameters.cc
	)

target_link_libraries (make_auto_workspace_h matpack methods)

add_custom_command (
	OUTPUT auto_workspace.h
	COMMAND make_auto_workspace_h
	DEPENDS make_auto_workspace_h
	)

add_custom_target (BuildAutoSources DEPENDS auto_workspace.h auto_md.h
                   auto_md.cc)

########### next target ###############

add_executable (make_auto_md_h
	agenda_record.cc
	arts.cc
  auto_workspace.h
	file.cc
	make_auto_md_h.cc
	messages.cc
	parameters.cc
	workspace.cc
	workspace_ng.cc
	)

target_link_libraries (make_auto_md_h methods matpack)

add_custom_command (
	OUTPUT auto_md.h
	COMMAND make_auto_md_h
	DEPENDS make_auto_md_h
)

########### next target ###############

add_executable (make_auto_md_cc
	agenda_record.cc
	arts.cc
  auto_workspace.h
	file.cc
	make_auto_md_cc.cc
	messages.cc
	parameters.cc
	workspace.cc
	workspace_ng.cc
	)

target_link_libraries (make_auto_md_cc methods matpack)

add_custom_command (
	OUTPUT auto_md.cc
	COMMAND make_auto_md_cc
	DEPENDS make_auto_md_cc
	)

########### next target ###############

set (ARTSCORE_EXTRA_SOURCES)
if (NOT HAVE_GETOPT_H)
  list (APPEND ARTSCORE_EXTRA_SOURCES
    arts_getopt.c arts_getopt1.c)
endif (NOT HAVE_GETOPT_H)

add_library (artscore STATIC
	auto_md.h
	auto_md.cc
	auto_workspace.h
  ${ARTSCORE_EXTRA_SOURCES}
  abs_species_tags.cc
  absorption.cc
  agenda_class.cc
  agenda_record.cc
  arts.cc
  arts_omp.cc
  bifstream.cc
  binio.cc
  bofstream.cc
  check_input.cc
  cia.cc
  cloudbox.cc
  constants.cc
  continua.cc
  disort.cc
  docserver.cc
  doit.cc
  Faddeeva.cc
  fastem.cc
  file.cc
  gas_abs_lookup.cc
  geomag_calc.cc
  geodetic.cc
  gridded_fields.cc
  gzstream.cc
  interpolation.cc
  interpolation_poly.cc
  jacobian.cc
  legendre.cc
  lin_alg.cc
  linemixingrecord.cc
  linemixingdata.cc
  linerecord.cc
  lineshapes.cc
  m_abs.cc
  m_abs_lookup.cc
  m_agenda.cc
  m_atmosphere.cc
  m_basic_types.cc
  m_batch.cc
  m_checked.cc
  m_cia.cc
  m_cloudbox.cc
  m_cloudradar.cc
  m_disort.cc
  m_doit.cc
  m_fos.cc
  m_general.cc
  m_geodetic.cc
  m_jacobian.cc
  m_linemixing.cc
  m_linerecord.cc
  m_montecarlo.cc
  m_optproperties.cc
  m_physics.cc
  m_planets.cc
  m_ppath.cc
  m_refraction.cc
  m_rte.cc
  m_sensor.cc
  m_surface.cc
  m_tmatrix.cc
  m_transmitter.cc
  m_xml.cc
  m_zeeman.cc
  math_funcs.cc
  mc_antenna.cc
  mc_interp.cc
  messages.cc
  montecarlo.cc
  nc_io.cc
  nc_io_array_types.cc
  nc_io_basic_types.cc
  nc_io_compound_types.cc
  oem.cc
  optproperties.cc
  parameters.cc
  parser.cc
  partition_function_data.cc
  linerecord_partitionfunctiondata.cc
  physics_funcs.cc
  poly_roots.cc
  ppath.cc
  pressurebroadeningdata.cc
  quantum.cc
  quantum_parser_hitran.cc
  rational.cc
  refraction.cc
  rng.cc
  rte.cc
  sensor.cc
  sourcetext.cc
  special_interp.cc
  species_data.cc
  surface.cc
  tmatrix.cc
  token.cc
  workspace.cc
  workspace_ng.cc
  xml_io.cc
  xml_io_array_types.cc
  xml_io_basic_types.cc
  xml_io_compound_types.cc
  zeeman.cc
	)
add_dependencies (artscore auto_version_h)

target_link_libraries (artscore ${ZLIB_LIBRARIES})

if (NETCDF_FOUND)
  target_link_libraries (artscore ${NETCDF_LIBRARIES})
endif (NETCDF_FOUND)

if (ENABLE_DOCSERVER)
  target_link_libraries (artscore microhttpd)
endif (ENABLE_DOCSERVER)

########### next target ###############

add_library (matpack STATIC
	complex.cc
	lin_alg.cc
	logic.cc
	make_vector.cc
	rational.cc
	matpackI.cc
	matpackII.cc
	matpackIII.cc
	matpackIV.cc
	matpackV.cc
	matpackVI.cc
	matpackVII.cc
	)

target_link_libraries( matpack ${BLAS_LIBRARIES} )

########### next target ###############

add_library (methods STATIC
	agendas.cc
	agenda_record.cc
	methods.cc
	methods_aux.cc
	groups.cc
	)

########### next target ###############

add_library (test_utils STATIC
        test_utils.cc
	)

########### testcases ###############

add_executable (test_binaryio test_binaryio.cc)

target_link_libraries (test_binaryio ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_gridded_fields
  gridded_fields.cc
  test_gridded_fields.cc)

target_link_libraries (test_gridded_fields matpack)

########### next testcase ###############

add_executable (test_interpolation test_interpolation.cc)
target_link_libraries (test_interpolation ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_legendre
  constants.cc
  math_funcs.cc
  legendre.cc
  test_legendre.cc)

target_link_libraries (test_legendre matpack)

########### next testcase ###############

add_executable (test_linalg
  constants.cc
  make_vector.cc
  lin_alg.cc
  test_linalg.cc)

target_link_libraries (test_linalg ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_integration
  constants.cc
  math_funcs.cc
  test_integration.cc)

target_link_libraries (test_integration matpack)

########### next testcase ###############

set ( MATLABDIR OFF CACHE PATH
  "Path to the Matlab install folder on the system." )
set ( ATMLABDIR OFF CACHE PATH
  "Path to the Atmlab install folder on the system" )

find_program( CSHFOUND csh )

if ( MATLABDIR AND ATMLABDIR AND CSHFOUND )
  add_executable (test_oem
    oem.cc
    test_oem.cc)

  add_definitions(-DSOURCEDIR="${CMAKE_BINARY_DIR}/src")
  add_definitions(-DATMLABDIR="${ATMLABDIR}")
  FIND_LIBRARY(ENGVAR libeng.so ${MATLABDIR}/R2015a/bin/glnxa64)
  FIND_LIBRARY(MXVAR libmx.so ${MATLABDIR}/R2015a/bin/glnxa64)
  target_link_libraries (test_oem matpack test_utils ${MXVAR} ${ENGVAR})
  FILE(MAKE_DIRECTORY test_oem_files)
  FILE(COPY test_oem/test_oem.m test_oem/linear_forward_model.m test_oem/make_plot.m
  test_oem/test_mult.m test_oem/test_inv.m
  DESTINATION test_oem_files)

endif ( MATLABDIR AND ATMLABDIR AND CSHFOUND )
########### next testcase ###############

add_executable (test_omp
  test_omp.cc)

target_link_libraries (test_omp  ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_poly poly_roots.h poly_roots.cc test_poly.cc)
target_link_libraries (test_poly matpack)

########### next testcase ###############

add_executable (test_tensor test_tensor.cc)
target_link_libraries (test_tensor matpack)

########### next testcase ###############

add_executable (test_matpack describe.h describe.cc
                wigner_functions.cc test_matpack.cc)
target_link_libraries (test_matpack matpack test_utils)

########### next testcase ###############

add_executable (test_quantum test_quantum.cc)
target_link_libraries (test_quantum ${ALL_ARTS_LIBRARIES})

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

########### next testcase ###############

add_executable (test_cia test_cia.cc cia.h cia.cc m_cia.cc)
target_link_libraries (test_cia ${ALL_ARTS_LIBRARIES})

########### subdirs ###############

add_subdirectory (libmicrohttpd)

if (WITH_GUI)
  add_subdirectory (gui)
endif (WITH_GUI)

########### tests ###############

arts_test_cmdline("describe" -d Copy)
arts_test_cmdline("groups"   -g)
arts_test_cmdline("help" -h)
arts_test_cmdline("input" -i Index)
arts_test_cmdline("methods" -m all)
arts_test_cmdline("version" -v)
arts_test_cmdline("workspacevariables" -w all)

