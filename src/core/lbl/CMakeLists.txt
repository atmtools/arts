add_library(lbl STATIC
  lbl_data.cpp
  lbl_fwd.cpp
  lbl_hitran.cpp
  lbl_lineshape.cpp
  lbl_lineshape_linemixing.cpp
  lbl_lineshape_model.cpp
  lbl_lineshape_voigt_ecs.cpp
  lbl_lineshape_voigt_ecs_hartmann.cpp
  lbl_lineshape_voigt_ecs_makarov.cpp
  lbl_lineshape_voigt_lte.cpp
  lbl_lineshape_voigt_lte_mirrored.cpp
  lbl_lineshape_voigt_nlte.cpp
  lbl_temperature_model.cpp
  lbl_zeeman.cpp
)

target_link_libraries(lbl PUBLIC matpack atm Faddeeva predef quantum surface wigner util arts_options)
target_include_directories(lbl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
