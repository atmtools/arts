add_library(predef STATIC
    ELL07.cc
    MPM89.cc
    MPM2020.cc
    MT_CKD100.cc
    MT_CKD252.cc
    CKDMT350.cc
    PWR98.cc
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
