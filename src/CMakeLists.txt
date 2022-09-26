include_directories ( ${CMAKE_CURRENT_BINARY_DIR}
                      ${CMAKE_CURRENT_SOURCE_DIR}
                      ${CMAKE_CURRENT_SOURCE_DIR}/matpack
                      ${ZLIB_INCLUDE_DIR}
                      ${ARTS_SOURCE_DIR}/3rdparty/invlib/src
                      ${ARTS_SOURCE_DIR}/3rdparty/cdisort
                      ${ARTS_SOURCE_DIR}/3rdparty/Eigen
                      ${ARTS_SOURCE_DIR}/3rdparty/mdspan/include()
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
  species
  )
  
if (ENABLE_GUI)
  list(PREPEND ALL_ARTS_LIBRARIES artsgui)
endif()

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

add_custom_target (auto_version_cc DEPENDS UpdateAutoVersion)

add_custom_command (
  OUTPUT auto_version.cc
  COMMAND ${CMAKE_COMMAND} -D "ARTS_BINARY_DIR:STRING=\"${ARTS_BINARY_DIR}\"" -D "ARTS_SOURCE_DIR:STRING=\"${ARTS_SOURCE_DIR}\"" -P ${ARTS_SOURCE_DIR}/cmake/scripts/update_auto_version_cc.cmake
  DEPENDS ${ARTS_BINARY_DIR}/auto_version.txt
)

########### next target ###############

add_executable (arts main.cc)
add_dependencies(check-deps arts)

target_link_libraries (arts ${ALL_ARTS_LIBRARIES})

install (TARGETS arts RUNTIME DESTINATION bin)

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
                   auto_md.cc tokval.h)

########### next target ###############

add_executable (make_auto_md_h
        agenda_record.cc
        arts.cc
        file.cc
        make_auto_md_h.cc
        messages.cc
        parameters.cc
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
        )

target_link_libraries (make_auto_md_cc methods matpack)

add_custom_command (
        OUTPUT auto_md.cc
        COMMAND make_auto_md_cc
        DEPENDS make_auto_md_cc
        )

########### next target ###############

add_executable (make_default_gins
        make_default_gins.cc
        )

target_link_libraries (make_default_gins methods)

add_custom_command (
        OUTPUT default_gins.cc default_gins.h
        COMMAND make_default_gins
        DEPENDS make_default_gins
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
  absorption.cc
  absorptionlines.cc
  agenda_class.cc
  agenda_record.cc
  agenda_set.cc
  arts.cc
  arts_omp.cc
  artstime.cc
  atm.cc
  auto_version.cc
  callback.cc
  check_input.cc
  cia.cc
  cloudbox.cc
  compile_time_tests.cc
  covariance_matrix.cc
  default_gins.cc
  disort.cc
  docserver.cc
  doit.cc
  energylevelmap.cc
  fastem.cc
  predefined_absorption_models.cc
  file.cc
  gas_abs_lookup.cc
  gas_scattering.cc
  geodetic.cc
  geodeticZZZ.cc
  gridded_fields.cc
  gui_propmat.cc
  hitran_species.cc
  xsec_fit.cc
  igrf13.cc
  interpolation.cc
  jacobian.cc
  jpl_species.cc
  legendre.cc
  linemixing.cc
  linemixing_hitran.cc
  linescaling.cc
  lineshape.cc
  lineshapemodel.cc
  m_abs.cc
  m_absorptionlines.cc
  m_abs_lookup.cc
  m_agenda.cc
  m_agenda_set.cc
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
  m_propmat.cc
  m_gas_scattering.cc
  m_general.cc
  m_geodetic.cc
  m_geodeticZZZ.cc
  m_xsec_fit.cc
  m_jacobian.cc
  m_linemixing.cc
  m_microphysics.cc
  m_montecarlo.cc
  m_nlte.cc
  m_oem.cc
  m_optproperties.cc
  m_partfun.cc
  m_physics.cc
  m_planets.cc
  m_ppath.cc
  m_ppathZZZ.cc
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
        m_sun.cc
  m_surface.cc
  m_telsem.cc
  m_tessem.cc
  m_tmatrix.cc
  m_transmitter.cc
  m_wigner.cc
  m_xml.cc
  m_zeeman.cc
  m_callback.cc
  math_funcs.cc
  mc_antenna.cc
  mc_interp.cc
  messages.cc
  microphysics.cc
  minimize.cc
  montecarlo.cc
  nlte.cc
  optproperties.cc
  parameters.cc
  parser.cc
  partfun.cc
  physics_funcs.cc
  poly_roots.cc
  ppath.cc
  ppath_struct.cc
  ppathZZZ.cc
  propmat_field.cc
  psd.cc
  quantum_numbers.cc
  quantum_term_symbol.cc
  radiation_field.cc
  raw.cc
  refraction.cc
  rng.cc
  rt4.cc
  rte.cc
  sensor.cc
  sourcetext.cc
  special_interp.cc
  species_info.cc
  species_tags.cc
        sun.cc
  surface.cc
  telsem.cc
  tessem.cc
  timings.cc
  tmatrix.cc
  tokval_io.cc
  transmissionmatrix.cc
  wigner_functions.cc
  workspace_ng.cc
  workspace_memory_handler.cc
  wsv_aux_operator.cc
  xml_io.cc
  xml_io_array_types.cc
  xml_io_basic_types.cc
  xml_io_compound_types.cc
  zeeman.cc
  zeemandata.cc
  ${NETCDF_CC_FILES}
        m_fluxes.cc)
add_dependencies (artscore auto_version_cc species)
target_link_libraries(artscore auto_partfun)

if (ENABLE_PCH)
  target_precompile_headers(artscore PRIVATE pch_artscore.h)
endif()

if (ENABLE_GUI)
  add_subdirectory(gui)
  
  target_include_directories(artscore PUBLIC ${ARTS_SOURCE_DIR}/src
                                             ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui/
                                             ${ARTS_SOURCE_DIR}/3rdparty/gui/implot/
                                             ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui-filebrowser/
                                             ${ARTS_SOURCE_DIR}/3rdparty/)
  
  target_link_libraries (artscore artsgui)

  target_compile_definitions(artscore PUBLIC ARTS_GUI_ENABLED)
endif()


add_subdirectory(predefined)
target_link_libraries (artscore predef)

target_link_libraries (artscore xmliobase wigner)

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
        
include_directories ( SYSTEM ${ARTS_SOURCE_DIR}/3rdparty )

########### next target ###############

add_library (xmliobase STATIC
        messages.cc
        xml_io_base.cc
        xml_io_general_types.cc
        )

target_link_libraries (xmliobase binio matpack)

########### next target ###############

add_library (binio STATIC
        bifstream.cc
        bofstream.cc
        binio.cc
        gzstream.cc 
        )

set_source_files_properties (binio.cc PROPERTIES
                             COMPILE_FLAGS "-fno-strict-aliasing")
target_link_libraries (binio ${ZLIB_LIBRARIES})

########### next target ###############

add_library (methods STATIC
        methods.cc
        tokval.cc

        agendas.cc
        agenda_record.cc
        methods_aux.cc
        groups.cc
        workspace.cc
        wsv_aux.cc
        )

set_source_files_properties (methods.cc PROPERTIES COMPILE_FLAGS "-O0")
target_link_libraries(methods PUBLIC predef species)

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
  artstime.cc
  gridded_fields.cc
  test_gridded_fields.cc)

target_link_libraries (test_gridded_fields matpack)

########### next testcase ###############

add_executable (test_interpolation test_interpolation.cc)
target_link_libraries (test_interpolation ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_atm test_atm.cc)
target_link_libraries (test_atm ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_linalg
  test_linalg.cc)

target_link_libraries (test_linalg ${ALL_ARTS_LIBRARIES}
  ${LAPACK_LIBRARIES} test_utils)

########### next testcase ###############

add_executable (test_integration
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

add_executable (test_matpack
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
target_link_libraries (test_xml xmliobase)

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

add_executable (test_propagationmatrix test_propagationmatrix.cc)
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

########### next testcase ###############

add_executable (test_hitran test_hitran.cc)
target_link_libraries(test_hitran ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

if (ENABLE_GUI)
  add_executable (test_lineshape test_lineshape.cc)
  target_link_libraries (test_lineshape ${ALL_ARTS_LIBRARIES})
endif ()

########### next testcase ###############

add_executable (test_magfield test_magfield.cc)
target_link_libraries (test_magfield ${ALL_ARTS_LIBRARIES})

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

# Build a test for the plotting tool
if (ENABLE_GUI)
  add_executable(test_gui test_gui.cc)
  target_link_libraries(test_gui PUBLIC ${ALL_ARTS_LIBRARIES})
endif()
########################################################################################

########################################################################################
add_subdirectory(partfun)

add_library(species STATIC
            isotopologues.cc
            )
target_link_libraries(species PUBLIC auto_partfun)
########################################################################################

########################################################################################
add_executable(test_species_tags test_species_tags.cc)
target_link_libraries (test_species_tags ${ALL_ARTS_LIBRARIES})
########################################################################################

########################################################################################
add_executable(test_isotope test_isotope.cc)
target_link_libraries (test_isotope species matpack)
########################################################################################

########################################################################################
# Test that the builtin partition functions are complete
add_executable(test_partfun test_partfun.cc partfun.cc)
target_link_libraries(test_partfun ${ALL_ARTS_LIBRARIES})
add_dependencies(check-deps test_partfun)
add_test(NAME "arts.cpp.fast.test_partfun" COMMAND test_partfun)
########################################################################################

########################################################################################
# Test that the implentation for predefined models are complete
add_executable(test_predefined test_predefined.cc)
target_link_libraries(test_predefined ${ALL_ARTS_LIBRARIES})
add_dependencies(check-deps test_predefined)
add_test(NAME "arts.cpp.fast.test_predefined" COMMAND test_predefined)
########################################################################################

########### next testcase ###############

add_executable (linemixing_jmh linemixing_jmh.cc)
target_link_libraries(linemixing_jmh ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

if (ENABLE_GUI)
add_executable (calc_predef predefined.cc)
target_link_libraries(calc_predef ${ALL_ARTS_LIBRARIES})
endif()

add_executable (arts_qn arts_quantum_numbers.cc quantum_numbers.cc)
target_link_libraries(arts_qn species xmliobase)

########### python interface
add_subdirectory(python_interface)

########### variant tokval
add_executable (make_tokval make_tokval.cc groups.cc)
add_custom_command(OUTPUT tokval.h tokval.cc
                   COMMAND make_tokval
                   DEPENDS make_tokval)

########################################################################################
pybind11_add_module(test_pyversion test_pyversion.cc)
add_custom_target(check-pyversion
  COMMAND ${Python3_EXECUTABLE} -c "import test_pyversion\\; test_pyversion.check_pyversion\\(\\)"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS test_pyversion
  COMMENT "Check Python version")
########################################################################################

########### IPO / LTO support
if (IPO_SUPPORTED)
  set_property(TARGET arts PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET artscore PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET binio PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET make_auto_md_cc PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET make_auto_md_h PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET matpack PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET methods PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET species PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET xmliobase PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

########################################################################################
# Test OpenMP
add_executable(test_omp EXCLUDE_FROM_ALL test_omp.cc)
target_link_libraries(test_omp PUBLIC ${ALL_ARTS_LIBRARIES})
########################################################################################
# Test Eigen
add_executable(test_matpack_eigen test_matpack_eigen.cc)
target_link_libraries(test_matpack_eigen PUBLIC matpack)
add_test(NAME "arts.cpp.fast.test_matpack_eigen" COMMAND test_matpack_eigen)
add_dependencies(check-deps test_matpack_eigen)
########################################################################################

add_executable(run_partfun run_partfun.cc)
target_link_libraries(run_partfun PRIVATE artscore)

#####
add_subdirectory(matpack)
add_executable(test_mdspan test_mdspan.cc)
target_link_libraries(test_mdspan PUBLIC ${LAPACK_LIBRARIES} matpack)
add_test(NAME "arts.cpp.fast.test_mdspan" COMMAND test_mdspan)
add_dependencies(check-deps test_mdspan)

#####
add_executable(test_matpack_perf test_matpack_perf.cc artstime.cc)
target_link_libraries(test_matpack_perf PUBLIC matpack)

#####
add_executable(test_interp_perf test_interp_perf.cc artstime.cc)
target_link_libraries(test_interp_perf PUBLIC matpack artscore)
