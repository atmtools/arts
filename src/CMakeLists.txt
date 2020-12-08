include_directories ( ${CMAKE_CURRENT_BINARY_DIR}
                      ${CMAKE_CURRENT_SOURCE_DIR}
                      ${ZLIB_INCLUDE_DIR}
                      ${ARTS_SOURCE_DIR}/3rdparty/invlib/src
                      ${ARTS_SOURCE_DIR}/3rdparty/cdisort
                    )

if (NETCDF_FOUND)
  include_directories (${NETCDF_INCLUDE_DIR})
  set (NETCDF_CC_FILES
    nc_io.cc
    nc_io_array_types.cc
    nc_io_basic_types.cc
    nc_io_compound_types.cc
    )

endif (NETCDF_FOUND)

if (ENABLE_MPI)
  include_directories (${MPI_CXX_INCLUDE_PATH})
endif (ENABLE_MPI)

set (ALL_ARTS_LIBRARIES
  artscore
  methods
  matpack
  cdisort
  )

file (GLOB HEADERFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
add_custom_target(UtilityHeadersArts SOURCES ${HEADERFILES})

if (OpenMP_CXX_FOUND)
  list (APPEND ALL_ARTS_LIBRARIES OpenMP::OpenMP_CXX)
endif ()

if (ENABLE_RT4)
  list (APPEND ALL_ARTS_LIBRARIES rt4)
endif (ENABLE_RT4)

if (ENABLE_FASTEM)
  list (APPEND ALL_ARTS_LIBRARIES fastem)
endif (ENABLE_FASTEM)

if (ENABLE_TMATRIX)
  list (APPEND ALL_ARTS_LIBRARIES tmatrix)
endif (ENABLE_TMATRIX)

if (ENABLE_REFICE)
  list (APPEND ALL_ARTS_LIBRARIES refice)
endif (ENABLE_REFICE)

list (APPEND ALL_ARTS_LIBRARIES Faddeeva)

########### target to extract version number from ChangeLog ###############

add_custom_target (auto_version_h DEPENDS auto_version.h)
add_dependencies (auto_version_h UpdateAutoVersion)

add_custom_command (
  OUTPUT auto_version.h
  COMMAND ${CMAKE_COMMAND} -D "ARTS_BINARY_DIR:STRING=\"${ARTS_BINARY_DIR}\"" -D "ARTS_SOURCE_DIR:STRING=\"${ARTS_SOURCE_DIR}\"" -P ${ARTS_SOURCE_DIR}/cmake/scripts/update_auto_version_h.cmake
  DEPENDS ${ARTS_BINARY_DIR}/auto_version.txt
)

########### next target ###############

add_executable (arts main.cc)
add_dependencies (arts auto_version_h)

target_link_libraries (arts ${ALL_ARTS_LIBRARIES})

install (TARGETS arts RUNTIME DESTINATION bin)

set_source_files_properties (legacy_continua.cc PROPERTIES
                             COMPILE_FLAGS "-fno-strict-aliasing")
set_source_files_properties (binio.cc PROPERTIES
                             COMPILE_FLAGS "-fno-strict-aliasing")

# Reduce compile time for very large initialization functions
set_source_files_properties (methods.cc PROPERTIES COMPILE_FLAGS "-O0")
set_source_files_properties (partition_function_data.cc PROPERTIES
                             COMPILE_FLAGS "-O0")
set_source_files_properties (species_data.cc PROPERTIES COMPILE_FLAGS "-O0")

# Silence invlib warning spam
set_source_files_properties (m_oem.cc PROPERTIES
                             COMPILE_FLAGS "-Wno-conversion")

########### next target ###############

if (C_API)
add_library (arts_api SHARED arts_api.cc interactive_workspace.cc arts_api_classes.cc)
add_dependencies (arts arts_api)
set_target_properties(arts_api PROPERTIES SUFFIX .so)
target_link_libraries (arts_api ${ALL_ARTS_LIBRARIES})
endif (C_API)

########### next target ###############

add_executable (make_workspace_memory_handler_cc
                make_workspace_memory_handler_cc.cc
                arts.cc
                groups.cc
                file.cc
                parameters.cc
                messages.cc)
add_custom_command (
        OUTPUT workspace_memory_handler.cc
        COMMAND make_workspace_memory_handler_cc
        DEPENDS make_workspace_memory_handler_cc
        )
add_custom_target (BuildAutoSources DEPENDS auto_md.h
                   auto_md.cc)

########### next target ###############

add_executable (make_auto_md_h
        agenda_record.cc
        arts.cc
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
  ${ARTSCORE_EXTRA_SOURCES}
  abs_species_tags.cc
  absorption.cc
  absorptionlines.cc
  agenda_class.cc
  agenda_record.cc
  arts.cc
  arts_omp.cc
  artstime.cc
  bifstream.cc
  binio.cc
  bofstream.cc
  check_input.cc
  cia.cc
  cloudbox.cc
  constants.cc
  covariance_matrix.cc
  disort.cc
  docserver.cc
  doit.cc
  energylevelmap.cc
  fastem.cc
  predefined_absorption_models.cc
  file.cc
  gas_abs_lookup.cc
  geodetic.cc
  gridded_fields.cc
  gzstream.cc
  hitran_xsec.cc
  interpolation.cc
  interpolation_poly.cc
  jacobian.cc
  legacy_continua.cc
  legendre.cc
  lin_alg.cc
  linemixing.cc
  linemixing_hitran.cc
  linescaling.cc
  lineshapemodel.cc
  linefunctions.cc
  m_abs.cc
  m_absorptionlines.cc
  m_abs_lookup.cc
  m_agenda.cc
  m_array.cc
  m_artstime.cc
  m_atmosphere.cc
  m_basic_types.cc
  m_batch.cc
  m_checked.cc
  m_cia.cc
  m_cloudbox.cc
  m_cloudradar.cc
  m_disort.cc
  m_doit.cc
  m_fluxes.cc
  m_predefined_absorption_models.cc
  m_fos.cc
  m_general.cc
  m_geodetic.cc
  m_hitran_xsec.cc
  m_jacobian.cc
  m_linemixing.cc
  m_microphysics.cc
  m_montecarlo.cc
  m_nlte.cc
  m_oem.cc
  m_optproperties.cc
  m_physics.cc
  m_planets.cc
  m_ppath.cc
  m_psd.cc
  m_quantum.cc
  m_radiation_field.cc
  m_raw.cc
  m_refraction.cc
  m_retrieval.cc
  m_rt4.cc
  m_rte.cc
  m_retrieval.cc
  m_sensor.cc
  m_surface.cc
  m_telsem.cc
  m_tessem.cc
  m_tmatrix.cc
  m_transmitter.cc
  m_wigner.cc
  m_xml.cc
  m_zeeman.cc
  math_funcs.cc
  mc_antenna.cc
  mc_interp.cc
  messages.cc
  microphysics.cc
  montecarlo.cc
  nlte.cc
  optproperties.cc
  parameters.cc
  parser.cc
  partition_function_data.cc
  physics_funcs.cc
  poly_roots.cc
  ppath.cc
  propagationmatrix.cc
  propmat_field.cc
  psd.cc
  quantum.cc
  quantum_parser_hitran.cc
  radiation_field.cc
  rational.cc
  raw.cc
  refraction.cc
  rng.cc
  rt4.cc
  rte.cc
  sensor.cc
  sourcetext.cc
  special_interp.cc
  species_data.cc
  species_info.cc
  surface.cc
  telsem.cc
  tessem.cc
  timings.cc
  tmatrix.cc
  token.cc
  transmissionmatrix.cc
  wigner_functions.cc
  workspace.cc
  workspace_ng.cc
  workspace_memory_handler.cc
  xml_io.cc
  xml_io_array_types.cc
  xml_io_basic_types.cc
  xml_io_compound_types.cc
  zeeman.cc
  zeemandata.cc
  ${NETCDF_CC_FILES}
        m_fluxes.cc)
add_dependencies (artscore auto_version_h)

if (ENABLE_GUI)
  add_subdirectory(gui)
  target_include_directories(artscore PUBLIC ${ARTS_SOURCE_DIR}/src
                                             ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui/
                                             ${ARTS_SOURCE_DIR}/3rdparty/gui/implot/
                                             ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui-filebrowser/
                                             ${ARTS_SOURCE_DIR}/3rdparty/)
  target_link_libraries (artscore artsgui)
endif()

target_link_libraries (artscore ${ZLIB_LIBRARIES} wigner)

if (FFTW_FOUND)
  include_directories (${FFTW_INCLUDE_DIR})
  target_link_libraries (artscore ${FFTW_LIBRARIES})
endif()

if (NETCDF_FOUND)
  target_link_libraries (artscore ${NETCDF_LIBRARIES})
endif (NETCDF_FOUND)

if (ENABLE_MPI)
  target_link_libraries (artscore ${MPI_CXX_LIBRARIES})
endif (ENABLE_MPI)

if (ENABLE_DOCSERVER)
  target_link_libraries (artscore microhttpd ${CMAKE_THREAD_LIBS_INIT})
endif (ENABLE_DOCSERVER)

########### next target ###############

add_library (matpack STATIC
        complex.cc
        lin_alg.cc
        logic.cc
        rational.cc
        matpackI.cc
        matpackII.cc
        matpackIII.cc
        matpackIV.cc
        matpackV.cc
        matpackVI.cc
        matpackVII.cc
        )
        
include_directories ( SYSTEM ${ARTS_SOURCE_DIR}/3rdparty )
target_link_libraries( matpack ${LAPACK_LIBRARIES} )

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
target_link_libraries( test_utils matpack )

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
  lin_alg.cc
  test_linalg.cc)

target_link_libraries (test_linalg ${ALL_ARTS_LIBRARIES}
  ${LAPACK_LIBRARIES} test_utils)

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

if (MATLABDIR AND ATMLABDIR AND CSHFOUND AND OEM_SUPPORT)

  if (OpenMP_CXX_FOUND)
    add_definitions(-DOMP)
  endif ()

  add_executable (test_oem
    timings.cc
    test_oem.cc)

  include_directories(${MATLABDIR}/extern/include)

  # Macros for Matlab interface.
  add_definitions(-DSOURCEDIR="${ARTS_BINARY_DIR}/src")
  add_definitions(-DATMLABDIR="${ATMLABDIR}")

  FIND_LIBRARY(ENGVAR libeng.so ${MATLABDIR}/bin/glnxa64)
  FIND_LIBRARY(MXVAR libmx.so ${MATLABDIR}/bin/glnxa64)

  target_link_libraries (test_oem matpack test_utils ${MXVAR} ${ENGVAR})

  # Copy test scripts.
  FILE(MAKE_DIRECTORY test_oem_files)
  FILE(COPY test_oem/test_oem.m test_oem/linear_forward_model.m test_oem/make_plot.m
  test_oem/test_mult.m test_oem/forward_model.m test_oem/test_oem_gauss_newton.m
  test_oem/test_oem_levenberg_marquardt.m test_oem/test_inv.m
  DESTINATION test_oem_files)

endif()
########### next testcase ###############

add_executable (test_poly poly_roots.h poly_roots.cc test_poly.cc)
target_link_libraries (test_poly matpack)

########### next testcase ###############

add_executable (test_tensor test_tensor.cc)
target_link_libraries (test_tensor matpack)

########### next testcase ###############

add_executable (test_matpack describe.h describe.cc
                wigner_functions.cc test_matpack.cc)
target_link_libraries (test_matpack matpack artscore test_utils)

########### next testcase ###############

add_executable (test_matpack_simple test_matpack_simple.cc)
target_link_libraries (test_matpack_simple test_utils)

########### next testcase ###############

add_executable (test_sparse test_sparse.cc)
target_link_libraries (test_sparse ${ALL_ARTS_LIBRARIES} test_utils)

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

########### next testcase ###############

add_executable (test_propagationmatrix test_propagationmatrix.cc absorption.h jacobian.h linescaling.h propagationmatrix.h global_data.h)
target_link_libraries (test_propagationmatrix ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_time test_time.cc artstime.h)
target_link_libraries (test_time ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_covariance_matrix test_covariance_matrix.cc)
target_link_libraries(test_covariance_matrix test_utils ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_telsem test_telsem.cc)
target_link_libraries(test_telsem ${ALL_ARTS_LIBRARIES})

########### subdirs ###############

add_subdirectory (libmicrohttpd)

########### tests ###############

arts_test_cmdline("describe" -d Copy)
arts_test_cmdline("groups"   -g)
arts_test_cmdline("help" -h)
arts_test_cmdline("input" -i Index)
arts_test_cmdline("methods" -m all)
arts_test_cmdline("version" -v)
arts_test_cmdline("workspacevariables" -w all)
arts_test_cmdline("check-docs" -C)

########### ARTS Interface ###############

########################################################################################
# ARTS Library for including in C++ (must link to autoarts.h to get file)
add_library (public_arts_interface INTERFACE)
target_link_libraries(public_arts_interface INTERFACE ${ALL_ARTS_LIBRARIES})
target_include_directories(public_arts_interface INTERFACE ${ARTS_BINARY_DIR}
                                                           ${ARTS_BINARY_DIR}/src/
                                                           ${ARTS_SOURCE_DIR}/src/
                                                           ${ARTS_SOURCE_DIR}/3rdparty/
                                                           ${ARTS_SOURCE_DIR}/3rdparty/invlib/src)
add_dependencies(public_arts_interface make_autoarts_h_target)
########################################################################################

########################################################################################
# Generate the internal structure for the ARTS C++ interface
add_executable(make_autoarts_h make_autoarts_h.cc)
target_link_libraries(make_autoarts_h PUBLIC ${ALL_ARTS_LIBRARIES})
########################################################################################

########################################################################################
# Generate autoarts.h for the ARTS interface
add_custom_command(OUTPUT autoarts.h
                   COMMAND make_autoarts_h > ${ARTS_BINARY_DIR}/src/autoarts.h
                   DEPENDS make_autoarts_h)
add_custom_target(make_autoarts_h_target SOURCES ${ARTS_BINARY_DIR}/src/autoarts.h)
########################################################################################

########################################################################################
# Generate autoarts.h for the ARTS interface
add_executable(test_cpp_api test_cpp_api.cc)
target_link_libraries(test_cpp_api public_arts_interface)
add_test(NAME "arts.cpp_api.fast.silly_groundbased_water_test" COMMAND test_cpp_api)
########################################################################################

# Build a test for the plotting tool
if (ENABLE_GUI)
  add_executable(test_gui test_gui.cc)
  target_link_libraries(test_gui PUBLIC ${ALL_ARTS_LIBRARIES})
endif()
