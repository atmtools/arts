add_library(rtepack STATIC
  rtepack_mueller_matrix.cc
  rtepack_multitype.cc
  rtepack_propagation_matrix.cc
  rtepack_rtestep.cc
  rtepack_scattering.cc
  rtepack_source.cc
  rtepack_stokes_vector.cc
  rtepack_transmission.cc
)
target_include_directories(rtepack PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(rtepack PUBLIC matpack physics)
