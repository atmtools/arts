add_library(predef STATIC
    MPM2020.cc
    PWR2021.cc
    CKDMT350.cc
)

target_link_libraries(predef PUBLIC matpack Faddeeva)
