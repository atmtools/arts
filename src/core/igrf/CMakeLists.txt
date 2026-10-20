set(igrf_data_file ${ARTS_CAT_DATA_DIR}/mag-field/igrf.txt)
if (NOT EXISTS ${igrf_data_file})
  message(FATAL_ERROR "No IGRF coefficient file found at ${igrf_data_file}
You need to pass -DARTS_CAT_DATA_PATH=<path> to cmake to point to the directory containing your copy of arts-cat-data.")
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
