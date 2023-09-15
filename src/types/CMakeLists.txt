##### Give Index, Numeric, String, and Array<> types

add_library(coretypes STATIC
  array.cpp
  mystring.cpp
)
target_include_directories(coretypes PUBLIC "${ARTS_SOURCE_DIR}/3rdparty")
target_include_directories(coretypes PUBLIC "${ARTS_BINARY_DIR}/src")
target_include_directories(coretypes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

######
