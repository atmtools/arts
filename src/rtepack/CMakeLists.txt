add_library(rtepack STATIC rtepack_transmission.cc rtepack_multitype.cc)
target_include_directories(rtepack PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
