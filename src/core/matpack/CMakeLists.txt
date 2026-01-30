add_library(matpack STATIC
  xml_io_stream_matpack_mdspan.cc
  xml_io_stream_matpack_mdspan_helpers.cc
  matpack_mdspan_helpers_gridded_data_t.cc
  xml_io_stream_matpack_range.cc
  lin_alg.cc
  matpack_sparse.cc
  rational.cc
  matpack_mdspan_math.cc
  matpack_mdspan_helpers.cc
  math_funcs.cc
  matpack_mdspan_helpers_band_matrix.cc
  matpack_mdspan_helpers_eigen.cc
  matpack_mdspan_view_t.cc
  interpolation.cc
  minimize.cc
  lagrange_interp.cc
  matpack_mdspan_algorithm.cc
)

find_package (Eigen3 3.4 REQUIRED NO_MODULE)
message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")

target_link_libraries(matpack PRIVATE ${LAPACK_LIBRARIES})
target_link_libraries(matpack PUBLIC coretypes util binio arts_enum_options time_report Eigen3::Eigen xml_io)
target_include_directories(matpack PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(matpack PUBLIC "${ARTS_SOURCE_DIR}/3rdparty/mdspan/include")

if (IPO_SUPPORTED)
  set_property(TARGET matpack PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_subdirectory(test)
