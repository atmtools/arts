# What RT3 (rt3/) and RT4 (rt4/) share, in the namespace polradtran: their
# streams and grounds (polradtran.h) and their layers from an ARTS path
# (polradtran_arts.h), always, and the C++ that their ports share
# (radintg.h, radutil.h, polradtran_workdata.h) when either is built
add_library(arts_polradtran STATIC polradtran.cc polradtran_arts.cc)
target_link_libraries(arts_polradtran PUBLIC matpack atm path rtepack PRIVATE scattering)
target_include_directories(arts_polradtran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(ENABLE_RT3 OR ENABLE_RT4)
  target_sources(arts_polradtran PRIVATE radintg.cc radutil.cc)
  target_link_libraries(arts_polradtran PRIVATE physics rtepack)
endif()

add_subdirectory(rt3)
add_subdirectory(rt4)
