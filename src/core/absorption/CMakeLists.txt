add_library(absorption STATIC
  cia.cc
  xsec_fit.cc
  zeemandata.cc
  predefined_absorption_models.cc
)
target_link_libraries(absorption PUBLIC path predef physics lbl)
target_include_directories(absorption PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
