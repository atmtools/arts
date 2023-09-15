add_library(griddedfield STATIC
  gridded_fields.cc
)
target_link_libraries(griddedfield PUBLIC matpack artstime)
target_include_directories(griddedfield PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
