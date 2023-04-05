if(ENABLE_ARTS_LGPL)
  add_library(artscore_gpl_selector STATIC gsl_gauss_legendre_lgpl.cc)
else()
  add_library(artscore_gpl_selector STATIC gsl_gauss_legendre_gpl.cc)
endif()

target_include_directories(artscore_gpl_selector PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
