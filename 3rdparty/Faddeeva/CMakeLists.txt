add_library (Faddeeva STATIC Faddeeva.cc)
arts_add_bundled_license(Faddeeva SPDX MIT FILES LICENSE)
target_include_directories(Faddeeva PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if (IPO_SUPPORTED)
  set_property(TARGET Faddeeva PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
