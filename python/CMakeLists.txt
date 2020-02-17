add_subdirectory(doc)

#
# Copy Python source files.
#
file(GLOB_RECURSE PYTHON_SRC_FILES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
set(ARTS_BINARY_DIR ${CMAKE_BINARY_DIR/src})
set(ARTS_SRC_DIR ${CMAKE_SOURCE_DIR})
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION ARTS_VERSION)
foreach(file ${PYTHON_SRC_FILES})
  configure_file(${file} ${file} @ONLY)
endforeach()

file(COPY test/reference DESTINATION test/)
file(COPY test/xml/reference DESTINATION test/xml/)
file(COPY test/test_files DESTINATION test/)
