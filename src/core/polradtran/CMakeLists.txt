if(ENABLE_RT3 OR ENABLE_RT4)
  # The C++ that the ports of RT3 (rt3/) and RT4 (rt4/) share, in the
  # namespace polradtran
  add_library(arts_polradtran STATIC radintg.cc radutil.cc)
  target_link_libraries(arts_polradtran PUBLIC matpack PRIVATE physics rtepack)
  target_include_directories(arts_polradtran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endif()

add_subdirectory(rt3)
add_subdirectory(rt4)
