add_library(predef STATIC
    ELL07.cc
    MPM89.cc
    MPM2020.cc
    PWR20xx.cc
    PWR98.cc
    TRE05.cc
    predef_data.cc
    standard.cc
)

if (NOT ENABLE_ARTS_LGPL)
  add_library(predef_aer STATIC
    MT_CKD100.cc
    MT_CKD252.cc
    CKDMT350.cc
    MT_CKD400.cc)
  target_link_libraries(predef_aer PUBLIC matpack)
  target_include_directories(predef_aer PRIVATE ${ARTS_SOURCE_DIR}/3rdparty/)
  target_include_directories(predef_aer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
else()
  add_library(predef_aer STATIC no_aer.cc)
endif()

target_link_libraries(predef PUBLIC matpack Faddeeva predef_aer)
target_include_directories(predef PRIVATE ${ARTS_SOURCE_DIR}/3rdparty/)
target_include_directories(predef PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if (IPO_SUPPORTED)
  set_property(TARGET predef PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
