# What RT3 (rt3/) and RT4 (rt4/) share, in the namespace polradtran: their
# streams and grounds (polradtran.h), their layers from an ARTS path
# (polradtran_arts.h), and the C++ that their ports share (radintg.h,
# radutil.h, polradtran_workdata.h)
add_library(arts_polradtran STATIC polradtran.cc polradtran_arts.cc radintg.cc radutil.cc)
target_link_libraries(arts_polradtran PUBLIC matpack atm path rtepack PRIVATE scattering physics)
target_include_directories(arts_polradtran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(rt3)
add_subdirectory(rt4)
