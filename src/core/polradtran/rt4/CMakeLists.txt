add_library(arts_rt4 STATIC rt4.cc rt4_arts.cc)
target_link_libraries(arts_rt4 PUBLIC arts_polradtran matpack atm path rtepack scattering)
target_include_directories(arts_rt4 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(ENABLE_RT4)
  # The C++ port of RADTRANO (radtran4.cc) and of the radintg4.f and
  # radutil4.f routines (radintg4.cc, radutil4.cc); the Fortran library
  # rt4 is the reference the tests compare the port with
  target_sources(arts_rt4 PRIVATE radtran4.cc radintg4.cc radutil4.cc)
  target_link_libraries(arts_rt4 PRIVATE rt4 physics)
  target_compile_definitions(arts_rt4 PRIVATE ARTS_HAS_RT4=1)
  add_subdirectory(test)
endif()
