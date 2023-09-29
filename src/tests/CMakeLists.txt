add_library (test_utils STATIC
        test_utils.cc
        )
target_link_libraries( test_utils matpack )

########### testcases ###############

add_executable (test_binaryio test_binaryio.cc)

target_link_libraries (test_binaryio artsworkspace)

########### next testcase ###############

add_executable (test_gridded_fields
  test_gridded_fields.cc)

target_link_libraries (test_gridded_fields artscore)

########### next testcase ###############

add_executable (test_interpolation test_interpolation.cc)
target_link_libraries (test_interpolation artsworkspace)

########### next testcase ###############

add_executable (test_linalg
  test_linalg.cc)

target_link_libraries (test_linalg artsworkspace
  ${LAPACK_LIBRARIES} test_utils)

########### next testcase ###############

add_executable (test_integration
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

add_executable (test_poly test_poly.cc)
target_link_libraries (test_poly matpack)

########### next testcase ###############

add_executable (test_tensor test_tensor.cc)
target_link_libraries (test_tensor artsworkspace)

########### next testcase ###############

add_executable (test_matpack test_matpack.cc)
target_link_libraries (test_matpack artscore test_utils)

########### next testcase ###############

add_executable (test_matpack_simple test_matpack_simple.cc)
target_link_libraries (test_matpack_simple test_utils)

########### next testcase ###############

add_executable (test_sparse test_sparse.cc)
target_link_libraries (test_sparse artsworkspace test_utils)

########### next testcase ###############

add_executable (test_xml test_xml.cc)
target_link_libraries (test_xml xmliobase)
target_include_directories(test_xml PRIVATE ${ARTS_SOURCE_DIR}/src)

########### next testcase ###############

add_executable (test_complex test_complex.cc)
target_link_libraries (test_complex matpack)

########### next testcase ###############

add_executable (test_sorting test_sorting.cc)
target_link_libraries (test_sorting artsworkspace)

########### next testcase ###############

add_executable (test_readpp test_readpp.cc)
target_link_libraries (test_readpp artsworkspace)

########### next testcase ###############

add_executable (test_cia test_cia.cc)
target_link_libraries (test_cia artsworkspace)

########### next testcase ###############

add_executable (test_time test_time.cc)
target_link_libraries (test_time artsworkspace)

########### next testcase ###############

add_executable (test_covariance_matrix test_covariance_matrix.cc)
target_link_libraries(test_covariance_matrix test_utils artsworkspace)

########### next testcase ###############

add_executable (test_telsem test_telsem.cc)
target_link_libraries(test_telsem artsworkspace)

########### next testcase ###############

add_executable (test_hitran test_hitran.cc)
target_link_libraries(test_hitran artsworkspace)

########### next testcase ###############

add_executable (test_magfield test_magfield.cc)
target_link_libraries (test_magfield artsworkspace)

########### tests ###############

########################################################################################
add_executable(test_species_tags test_species_tags.cc)
target_link_libraries (test_species_tags artsworkspace)
########################################################################################

########################################################################################
add_executable(test_isotope test_isotope.cc)
target_link_libraries (test_isotope species matpack)
########################################################################################

########################################################################################
# Test that the builtin partition functions are complete
add_executable(test_partfun test_partfun.cc)
target_link_libraries(test_partfun artsworkspace partfun)
add_dependencies(check-deps test_partfun)
add_test(NAME "cpp.fast.test_partfun" COMMAND test_partfun)
########################################################################################

########################################################################################
# Test that the implentation for predefined models are complete
add_executable(test_predefined test_predefined.cc)
target_link_libraries(test_predefined artsworkspace)
add_dependencies(check-deps test_predefined)
add_test(NAME "cpp.fast.test_predefined" COMMAND test_predefined)
########################################################################################

########### next testcase ###############

add_executable (linemixing_jmh linemixing_jmh.cc)
target_link_libraries(linemixing_jmh artsworkspace)

########### next testcase ###############

add_executable (arts_qn arts_quantum_numbers.cc)
target_link_libraries(arts_qn species xmliobase quantum)

########################################################################################
# Test OpenMP
add_executable(test_omp EXCLUDE_FROM_ALL test_omp.cc)
target_link_libraries(test_omp PUBLIC artsworkspace)
########################################################################################
# Test Eigen
add_executable(test_matpack_eigen test_matpack_eigen.cc)
target_link_libraries(test_matpack_eigen PUBLIC matpack)
add_test(NAME "cpp.fast.test_matpack_eigen" COMMAND test_matpack_eigen)
add_dependencies(check-deps test_matpack_eigen)
########################################################################################

add_executable(run_partfun run_partfun.cc)
target_link_libraries(run_partfun PRIVATE artsworkspace)

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
add_executable(test_matpack_perf test_matpack_perf.cc)
target_link_libraries(test_matpack_perf PUBLIC artscore)
target_include_directories(test_matpack_perf PRIVATE ${ARTS_SOURCE_DIR}/src)

#####
add_executable(test_interp_perf test_interp_perf.cc)
target_link_libraries(test_interp_perf PUBLIC matpack artsworkspace)

#####
add_executable(test_rng test_rng.cc)
target_link_libraries(test_rng PUBLIC artscore)
target_include_directories(test_rng PRIVATE ${ARTS_SOURCE_DIR}/src)

#####
add_executable(test_rtepack test_rtepack.cc)
target_link_libraries(test_rtepack PUBLIC artstime rtepack)
