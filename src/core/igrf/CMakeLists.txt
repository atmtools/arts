# The glob notices when the file appears in, or disappears from, arts-cat-data
file(GLOB igrf_data_file CONFIGURE_DEPENDS ${ARTS_CAT_DATA_DIR}/mag-field/igrf.txt)
if (NOT igrf_data_file)
  message(WARNING "No IGRF coefficient file found at ${ARTS_CAT_DATA_DIR}/mag-field/igrf.txt
ARTS is compiled without the IGRF model, and using it is an error.  Update your copy of arts-cat-data to use it.")
endif()

add_executable(make_igrf_data make_igrf_data.cc)

arts_generate(make_igrf_data
  ARGS ${igrf_data_file}
  OUTPUT auto_igrf_data.h auto_igrf_data.cpp
  DEPENDS ${igrf_data_file}
  COMMENT "Generating automatic IGRF coefficients"
)

add_library(igrf STATIC
  igrf.cpp
  auto_igrf_data.cpp
)

target_link_libraries(igrf PUBLIC matpack geodesy artstime legendre)
target_include_directories(igrf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
