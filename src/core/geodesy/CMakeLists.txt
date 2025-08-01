add_library(geodesy STATIC
  geodetic.cpp
)

target_link_libraries(geodesy PUBLIC matpack)
target_include_directories(geodesy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
