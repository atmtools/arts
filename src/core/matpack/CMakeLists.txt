add_library(matpack STATIC
  double_imanip.cc
  lin_alg.cc
  matpack_sparse.cc
  rational.cc
  poly_roots.cc
  matpack_mdspan_math.cc
  matpack_mdspan_helpers.cc
  math_funcs.cc
  matpack_mdspan_helpers_band_matrix.cc
  matpack_mdspan_helpers_eigen.cc
)

find_package (Eigen3 3.4 REQUIRED NO_MODULE)
message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")

target_link_libraries(matpack PRIVATE ${LAPACK_LIBRARIES})
target_link_libraries(matpack PUBLIC coretypes util arts_enum_options Eigen3::Eigen)
target_include_directories(matpack PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(matpack PUBLIC "${ARTS_SOURCE_DIR}/3rdparty/mdspan/include")

if (IPO_SUPPORTED)
  set_property(TARGET matpack PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
