add_library(arts_rt3 STATIC rt3.cc rt3_arts.cc)
target_link_libraries(arts_rt3 PUBLIC matpack atm path rtepack scattering PRIVATE legendre)
target_include_directories(arts_rt3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(ENABLE_RT3)
  # The C++ port of RADTRAN (radtran3.cc) and of every routine it calls:
  # radscat3.f (radscat3.cc, and Evans' FFT in rt3_fft.cc, for which a build
  # may choose FFTW as a compile-time option), radutil3.f (radutil3.cc) and
  # radintg3.f (radintg3.cc).  It calls no Fortran; the Fortran library rt3
  # is the reference the tests compare the port with
  target_sources(arts_rt3 PRIVATE radtran3.cc radintg3.cc radscat3.cc rt3_fft.cc radutil3.cc)
  target_link_libraries(arts_rt3 PUBLIC arts_polradtran PRIVATE physics)
  target_compile_definitions(arts_rt3 PRIVATE ARTS_HAS_RT3=1)
  add_subdirectory(test)
endif()
