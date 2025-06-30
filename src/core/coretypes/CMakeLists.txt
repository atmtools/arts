# #### Give Index, Numeric, String, and Array<> types
set(CORETYPE_EXTRA)

add_library(coretypes STATIC
  ${CORETYPE_EXTRA}
  file.cc
  parameters.cc
  double_imanip.cc
)
target_include_directories(coretypes PUBLIC "${ARTS_SOURCE_DIR}/3rdparty")
target_include_directories(coretypes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(coretypes PUBLIC "${ARTS_BINARY_DIR}/src")
target_link_libraries(coretypes PUBLIC util strings)
