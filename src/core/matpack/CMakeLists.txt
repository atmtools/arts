add_library(matpack STATIC
  double_imanip.cc
  lin_alg.cc
  logic.cc
  matpack_math.cc
  matpack_sparse.cc
  rational.cc
  math_funcs.cc
  poly_roots.cc
)

target_link_libraries(matpack PRIVATE ${LAPACK_LIBRARIES})
target_link_libraries(matpack PUBLIC coretypes util arts_enum_options)
target_include_directories(matpack PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(matpack PUBLIC "${ARTS_SOURCE_DIR}/3rdparty/Eigen")
target_include_directories(matpack PUBLIC "${ARTS_SOURCE_DIR}/3rdparty/mdspan/include")

if (IPO_SUPPORTED)
  set_property(TARGET matpack PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
