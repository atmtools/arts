add_library (Faddeeva STATIC Faddeeva.cc)

if (IPO_SUPPORTED)
  set_property(TARGET Faddeeva PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
