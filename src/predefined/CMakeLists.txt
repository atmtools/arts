add_library(predef STATIC
    MPM89.cc
    MPM2020.cc
    CKDMT350.cc
    hitran_water_mtckd.cc
    predef_data.cc
)

if (ENABLE_PCH)
    target_precompile_headers(predef PRIVATE <pch_matpack.h>)
endif()
target_link_libraries(predef PUBLIC matpack)
target_include_directories(predef PRIVATE ${ARTS_SOURCE_DIR}/3rdparty/)

if (IPO_SUPPORTED)
  set_property(TARGET predef PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
