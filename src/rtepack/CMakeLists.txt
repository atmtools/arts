add_library(rtepack STATIC rtepack_transmission.cc rtepack_multitype.cc rtepack_mueller_matrix.cc rtepack_scattering.cc)
target_include_directories(rtepack PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(rtepack PUBLIC matpack)
