add_library (Faddeeva STATIC Faddeeva.cc)
target_include_directories(Faddeeva PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if (IPO_SUPPORTED)
  set_property(TARGET Faddeeva PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
