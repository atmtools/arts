add_library(predef STATIC
    MPM2020.cc
    CKDMT350.cc
    hitran_water_mtckd.cc
    predef_data.cc
)

if (ENABLE_PCH)
    target_precompile_headers(predef PRIVATE <pch_matpack.h>)
endif()
target_link_libraries(predef PUBLIC matpack)
