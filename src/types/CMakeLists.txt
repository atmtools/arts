##### Give Index, Numeric, String, and Array<> types

add_library(coretypes STATIC
  array.cpp
  mystring.cpp
)
target_include_directories(coretypes PUBLIC "${ARTS_SOURCE_DIR}/3rdparty")
target_include_directories(coretypes PUBLIC "${ARTS_BINARY_DIR}/src")
target_include_directories(coretypes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

######

add_subdirectory(binio)

add_subdirectory(matpack)
add_subdirectory(rtepack)
add_subdirectory(predefined)
add_subdirectory(operators)
add_subdirectory(fwd)
add_subdirectory(quantum)
add_subdirectory(spec)
add_subdirectory(griddedfield)
add_subdirectory(compound)
add_subdirectory(artstime)

add_library(artstypes INTERFACE)
target_link_libraries(artstypes INTERFACE
  artstime
  compoundtypes
  coretypes
  griddedfield
  fwd
  matpack
  operators
  predef
  quantum
  rtepack
  species
  species_tags
)
