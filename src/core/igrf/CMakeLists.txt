add_library(igrf STATIC
  igrf13.cpp
)

target_link_libraries(igrf PUBLIC matpack geodesy artstime legendre)
target_include_directories(igrf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
