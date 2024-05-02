add_library(test_utils STATIC
  test_utils.cc
)
target_link_libraries(test_utils matpack)

# ########## testcases ###############
add_executable(test_binaryio test_binaryio.cc)

target_link_libraries(test_binaryio artsworkspace)

# ########## next testcase ###############
add_executable(test_gridded_fields
  test_gridded_fields.cc)

target_link_libraries(test_gridded_fields matpack)

# ########## next testcase ###############
add_executable(test_interpolation test_interpolation.cc)
target_link_libraries(test_interpolation artsworkspace)

# ########## next testcase ###############
add_executable(test_linalg
  test_linalg.cc)

target_link_libraries(test_linalg artsworkspace
  ${LAPACK_LIBRARIES} test_utils)

# ########## next testcase ###############
add_executable(test_integration
  test_integration.cc)

target_link_libraries(test_integration PUBLIC artstime matpack)

# ########## next testcase ###############
set(MATLABDIR OFF CACHE PATH
  "Path to the Matlab install folder on the system.")
set(ATMLABDIR OFF CACHE PATH
  "Path to the Atmlab install folder on the system")

find_program(CSHFOUND csh)

if(MATLABDIR AND ATMLABDIR AND CSHFOUND AND OEM_SUPPORT)
  if(OpenMP_CXX_FOUND)
    add_definitions(-DOMP)
  endif()

  add_executable(test_oem
    timings.cc
    test_oem.cc)

  include_directories(${MATLABDIR}/extern/include)

  # Macros for Matlab interface.
  add_definitions(-DSOURCEDIR="${ARTS_BINARY_DIR}/src")
  add_definitions(-DATMLABDIR="${ATMLABDIR}")

  FIND_LIBRARY(ENGVAR libeng.so ${MATLABDIR}/bin/glnxa64)
  FIND_LIBRARY(MXVAR libmx.so ${MATLABDIR}/bin/glnxa64)

  target_link_libraries(test_oem matpack test_utils ${MXVAR} ${ENGVAR})

  # Copy test scripts.
  FILE(MAKE_DIRECTORY test_oem_files)
  FILE(COPY test_oem/test_oem.m test_oem/linear_forward_model.m test_oem/make_plot.m
    test_oem/test_mult.m test_oem/forward_model.m test_oem/test_oem_gauss_newton.m
    test_oem/test_oem_levenberg_marquardt.m test_oem/test_inv.m
    DESTINATION test_oem_files)
endif()

# ########## next testcase ###############
add_executable(test_poly test_poly.cc)
target_link_libraries(test_poly matpack)

# ########## next testcase ###############
add_executable(test_tensor test_tensor.cc)
target_link_libraries(test_tensor artsworkspace)

# ########## next testcase ###############
add_executable(test_matpack test_matpack.cc)
target_link_libraries(test_matpack artscore test_utils)

# ########## next testcase ###############
add_executable(test_matpack_simple test_matpack_simple.cc)
target_link_libraries(test_matpack_simple test_utils)

# ########## next testcase ###############
add_executable(test_sparse test_sparse.cc)
target_link_libraries(test_sparse artsworkspace test_utils)

# ########## next testcase ###############
add_executable(test_xml test_xml.cc)
target_link_libraries(test_xml xmliobase)
target_include_directories(test_xml PRIVATE ${ARTS_SOURCE_DIR}/src)

# ########## next testcase ###############
add_executable(test_complex test_complex.cc)
target_link_libraries(test_complex matpack)

# ########## next testcase ###############
add_executable(test_sorting test_sorting.cc)
target_link_libraries(test_sorting artsworkspace)

# ########## next testcase ###############
add_executable(test_readpp test_readpp.cc)
target_link_libraries(test_readpp artsworkspace)

# ########## next testcase ###############
add_executable(test_cia test_cia.cc)
target_link_libraries(test_cia artsworkspace)

# ########## next testcase ###############
add_executable(test_time test_time.cc)
target_link_libraries(test_time artsworkspace)

# ########## next testcase ###############
add_executable(test_telsem test_telsem.cc)
target_link_libraries(test_telsem artsworkspace)

# ########## next testcase ###############
add_executable(test_hitran test_hitran.cc)
target_link_libraries(test_hitran artsworkspace)

# ########## next testcase ###############
add_executable(test_magfield test_magfield.cc)
target_link_libraries(test_magfield artsworkspace)

# ########## tests ###############

# #######################################################################################
add_executable(test_species_tags test_species_tags.cc)
target_link_libraries(test_species_tags artsworkspace)

# #######################################################################################

# #######################################################################################
add_executable(test_isotope test_isotope.cc)
target_link_libraries(test_isotope species matpack)

# #######################################################################################

# #######################################################################################
# Test that the builtin partition functions are complete
add_executable(test_partfun test_partfun.cc)
target_link_libraries(test_partfun artsworkspace partfun)
add_dependencies(check-deps test_partfun)
add_test(NAME "cpp.fast.test_partfun" COMMAND test_partfun)

# #######################################################################################

# #######################################################################################
# Test that the implentation for predefined models are complete
add_executable(test_predefined test_predefined.cc)
target_link_libraries(test_predefined artsworkspace)
add_dependencies(check-deps test_predefined)
add_test(NAME "cpp.fast.test_predefined" COMMAND test_predefined)

# #######################################################################################

# ########## next testcase ###############
add_executable(linemixing_jmh linemixing_jmh.cc)
target_link_libraries(linemixing_jmh artsworkspace)

# ########## next testcase ###############
add_executable(arts_qn arts_quantum_numbers.cc)
target_link_libraries(arts_qn species xmliobase quantum)

# #######################################################################################
# Test OpenMP
add_executable(test_omp EXCLUDE_FROM_ALL test_omp.cc)
target_link_libraries(test_omp PUBLIC artsworkspace)

# #######################################################################################
# Test Eigen
add_executable(test_matpack_eigen test_matpack_eigen.cc)
target_link_libraries(test_matpack_eigen PUBLIC matpack)
add_test(NAME "cpp.fast.test_matpack_eigen" COMMAND test_matpack_eigen)
add_dependencies(check-deps test_matpack_eigen)

# #######################################################################################
add_executable(run_partfun run_partfun.cc)
target_link_libraries(run_partfun PRIVATE artsworkspace)

# ####
add_executable(test_mdspan test_mdspan.cc)
target_link_libraries(test_mdspan PUBLIC ${LAPACK_LIBRARIES} matpack artstime)
add_test(NAME "cpp.fast.test_mdspan" COMMAND test_mdspan)
add_dependencies(check-deps test_mdspan)

# #######################################################################################
pybind11_add_module(test_pyversion test_pyversion.cc)
add_custom_target(check-pyversion
  COMMAND ${Python3_EXECUTABLE} -c "import test_pyversion\\; test_pyversion.check_pyversion\\(\\)"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS test_pyversion
  COMMENT "Check Python version")

# #######################################################################################

# ####
add_executable(test_matpack_perf test_matpack_perf.cc)
target_link_libraries(test_matpack_perf PUBLIC artscore)
target_include_directories(test_matpack_perf PRIVATE ${ARTS_SOURCE_DIR}/src)

add_custom_target(
  run_matpack_perf
  COMMAND test_matpack_perf 10 100000000 100000000 8000 2000 100000000 9000 100000000 6000 > matpack_perf.txt
  DEPENDS test_matpack_perf
  BYPRODUCTS matpack_perf.txt
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Running performance test for matpack"
)

# ####
add_executable(test_interp_perf test_interp_perf.cc)
target_link_libraries(test_interp_perf PUBLIC matpack artsworkspace)

add_custom_target(
  run_interp_perf
  COMMAND test_interp_perf 10 1000000 1000000 3000000 20000 > interp_perf.txt
  DEPENDS test_interp_perf
  BYPRODUCTS interp_perf.txt
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Running performance test for interpolation"
)

# ####
add_executable(test_rng test_rng.cc)
target_link_libraries(test_rng PUBLIC artscore)
target_include_directories(test_rng PRIVATE ${ARTS_SOURCE_DIR}/src)

# ####
add_executable(test_rtepack test_rtepack.cc)
target_link_libraries(test_rtepack PUBLIC artstime rtepack)

# ####
add_executable(test_path_point test_path_point.cc)
target_link_libraries(test_path_point PUBLIC artstime path)
add_test(NAME "cpp.fast.test_path_point" COMMAND test_path_point)
add_dependencies(check-deps test_path_point)

# ####
add_executable(test_fwd test_fwd.cc)
target_link_libraries(test_fwd PUBLIC fwd)
add_test(NAME "cpp.fast.test_fwd" COMMAND test_fwd)
add_dependencies(check-deps test_fwd)

# ###  Set up a bunch of performance tests
# ###  NOTE: New tests should be added as dependencies to the run_perf target,
# ###        but also to one-another so the tests are not run at the same time
# ###        (affecting performance, which is what we want to test, so we want to avoid that)
add_dependencies(run_interp_perf run_matpack_perf)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/perf_results.py perf_results.py COPYONLY)
add_custom_target(run_perf
  COMMAND ${Python3_EXECUTABLE} perf_results.py matpack_perf.txt interp_perf.txt > perf_report.rst
  DEPENDS run_matpack_perf run_interp_perf
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Creating performance test report"
)

# ####
add_executable(test_legendre test_legendre.cc)
target_link_libraries(test_legendre PUBLIC legendre)
add_test(NAME "cpp.fast.test_legendre" COMMAND test_legendre)
add_dependencies(check-deps test_legendre)


# ####
add_executable(test_einsum_perf test_einsum_perf.cc)
target_link_libraries(test_einsum_perf PUBLIC matpack artstime)
add_test(NAME "cpp.fast.test_einsum_perf" COMMAND test_einsum_perf)
add_dependencies(check-deps test_einsum_perf)


# ####
add_executable(test_band_matrix_solver test_band_matrix_solver.cc)
target_link_libraries(test_band_matrix_solver PUBLIC matpack lapack)
add_test(NAME "cpp.fast.test_band_matrix_solver" COMMAND test_band_matrix_solver)
add_dependencies(check-deps test_band_matrix_solver)
