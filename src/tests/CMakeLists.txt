add_library (test_utils STATIC
        test_utils.cc
        )
target_link_libraries( test_utils matpack )

########### testcases ###############

add_executable (test_binaryio test_binaryio.cc)

target_link_libraries (test_binaryio ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_gridded_fields
  ../gridded_fields.cc
  test_gridded_fields.cc)

target_link_libraries (test_gridded_fields matpack)

########### next testcase ###############

add_executable (test_interpolation test_interpolation.cc)
target_link_libraries (test_interpolation ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_linalg
  test_linalg.cc)

target_link_libraries (test_linalg ${ALL_ARTS_LIBRARIES}
  ${LAPACK_LIBRARIES} test_utils)

########### next testcase ###############

add_executable (test_integration
  ../math_funcs.cc
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

add_executable (test_poly ../poly_roots.h ../poly_roots.cc test_poly.cc)
target_link_libraries (test_poly matpack)

########### next testcase ###############

add_executable (test_tensor test_tensor.cc)
target_link_libraries (test_tensor matpack)

########### next testcase ###############

add_executable (test_matpack
                ../wigner_functions.cc test_matpack.cc)
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

add_executable (test_sorting test_sorting.cc)
target_link_libraries (test_sorting matpack)

########### next testcase ###############

add_executable (test_readpp test_readpp.cc)
target_link_libraries (test_readpp ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_cia test_cia.cc ../cia.cc ../m_cia.cc)
target_link_libraries (test_cia ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_propagationmatrix test_propagationmatrix.cc)
target_link_libraries (test_propagationmatrix ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

add_executable (test_time test_time.cc)
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

########### tests ###############

arts_test_cmdline("describe" -d Copy)
arts_test_cmdline("groups"   -g)
arts_test_cmdline("help" -h)
arts_test_cmdline("input" -i Index)
arts_test_cmdline("methods" -m all)
arts_test_cmdline("version" -v)
arts_test_cmdline("workspacevariables" -w all)
if (ENABLE_DOCSERVER)
  arts_test_cmdline("check-docs" -C)
endif ()

# Build a test for the plotting tool
if (ENABLE_GUI)
  add_executable(test_gui test_gui.cc)
  target_link_libraries(test_gui PUBLIC ${ALL_ARTS_LIBRARIES})
endif()

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
add_executable(test_partfun test_partfun.cc ../partfun.cc)
target_link_libraries(test_partfun ${ALL_ARTS_LIBRARIES})
add_dependencies(check-deps test_partfun)
add_test(NAME "cpp.fast.test_partfun" COMMAND test_partfun)
########################################################################################

########################################################################################
# Test that the implentation for predefined models are complete
add_executable(test_predefined test_predefined.cc)
target_link_libraries(test_predefined ${ALL_ARTS_LIBRARIES})
add_dependencies(check-deps test_predefined)
add_test(NAME "cpp.fast.test_predefined" COMMAND test_predefined)
########################################################################################

########### next testcase ###############

add_executable (linemixing_jmh linemixing_jmh.cc)
target_link_libraries(linemixing_jmh ${ALL_ARTS_LIBRARIES})

########### next testcase ###############

if (ENABLE_GUI)
add_executable (calc_predef predefined.cc)
target_link_libraries(calc_predef ${ALL_ARTS_LIBRARIES})
endif()

add_executable (arts_qn arts_quantum_numbers.cc ../quantum_numbers.cc)
target_link_libraries(arts_qn species xmliobase)

########################################################################################
# Test OpenMP
add_executable(test_omp EXCLUDE_FROM_ALL test_omp.cc)
target_link_libraries(test_omp PUBLIC ${ALL_ARTS_LIBRARIES})
########################################################################################
# Test Eigen
add_executable(test_matpack_eigen test_matpack_eigen.cc)
target_link_libraries(test_matpack_eigen PUBLIC matpack)
add_test(NAME "cpp.fast.test_matpack_eigen" COMMAND test_matpack_eigen)
add_dependencies(check-deps test_matpack_eigen)
########################################################################################

add_executable(run_partfun run_partfun.cc)
target_link_libraries(run_partfun PRIVATE artscore)

#####
add_executable(test_mdspan test_mdspan.cc)
target_link_libraries(test_mdspan PUBLIC ${LAPACK_LIBRARIES} matpack)
add_test(NAME "cpp.fast.test_mdspan" COMMAND test_mdspan)
add_dependencies(check-deps test_mdspan)

########################################################################################
pybind11_add_module(test_pyversion test_pyversion.cc)
add_custom_target(check-pyversion
  COMMAND ${Python3_EXECUTABLE} -c "import test_pyversion\\; test_pyversion.check_pyversion\\(\\)"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS test_pyversion
  COMMENT "Check Python version")
########################################################################################

#####
add_executable(test_matpack_perf test_matpack_perf.cc ../artstime.cc)
target_link_libraries(test_matpack_perf PUBLIC matpack)

#####
add_executable(test_interp_perf test_interp_perf.cc ../artstime.cc)
target_link_libraries(test_interp_perf PUBLIC matpack artscore)

#####
add_executable(test_rng test_rng.cc artstime.cc)
target_link_libraries(test_rng PUBLIC matpack)
