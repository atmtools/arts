##### Give Index, Numeric, String, and Array<> types
set (CORETYPE_EXTRA)
if (NOT HAVE_GETOPT_H)
  list (APPEND CORETYPE_EXTRA
    arts_getopt.c arts_getopt1.c)
endif (NOT HAVE_GETOPT_H)

add_library(coretypes STATIC
  ${CORETYPE_EXTRA}
  array.cpp
  file.cc
  mystring.cpp
  parameters.cc
)
target_include_directories(coretypes PUBLIC "${ARTS_SOURCE_DIR}/3rdparty")
target_include_directories(coretypes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(coretypes PUBLIC "${ARTS_BINARY_DIR}/src")

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
