if(ENABLE_ARTS_LGPL)
  add_library(artsworkspace_gpl_selector STATIC gsl_gauss_legendre_lgpl.cc)
else()
  add_library(artsworkspace_gpl_selector STATIC gsl_gauss_legendre_gpl.cc)
endif()

target_include_directories(artsworkspace_gpl_selector PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(artsworkspace_gpl_selector PUBLIC matpack)
