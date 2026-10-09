# Evans' RT4: the C++ port of RADTRANO (radtran4.cc) and of the routines of
# radintg4.f and radutil4.f that RT3 does not share (radintg4.cc,
# radutil4.cc), behind rt4::solve (rt4.cc) and the inputs from ARTS data
# (rt4_arts.cc)
add_library(arts_rt4 STATIC rt4.cc rt4_arts.cc radtran4.cc radintg4.cc radutil4.cc)
target_link_libraries(arts_rt4 PUBLIC arts_polradtran matpack atm path rtepack scattering PRIVATE physics)
target_include_directories(arts_rt4 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(test)
