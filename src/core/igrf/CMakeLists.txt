add_library(igrf STATIC
  igrf13.cpp
  geodetic.cpp
)

target_link_libraries(igrf PUBLIC matpack artstime legendre)
target_include_directories(igrf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
