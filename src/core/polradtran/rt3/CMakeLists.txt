# Evans' RT3: the C++ port of RADTRAN (radtran3.cc) and of every routine it
# calls, those of radscat3.f (radscat3.cc, and Evans' FFT in rt3_fft.cc),
# radutil3.f (radutil3.cc) and radintg3.f (radintg3.cc) that RT4 does not
# share, behind rt3::solve (rt3.cc) and the inputs from ARTS data (rt3_arts.cc)
add_library(arts_rt3 STATIC rt3.cc rt3_arts.cc radtran3.cc radintg3.cc radscat3.cc rt3_fft.cc radutil3.cc)
target_link_libraries(arts_rt3 PUBLIC arts_polradtran matpack atm path rtepack scattering PRIVATE legendre physics)
target_include_directories(arts_rt3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(test)
