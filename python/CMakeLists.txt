add_subdirectory(doc)

#
# Copy Python source files.
#

file(GLOB_RECURSE PYTHON_SRC_FILES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
set(ARTS_BINARY_DIR ${CMAKE_BINARY_DIR/src})
set(ARTS_SRC_DIR ${CMAKE_SOURCE_DIR})
if (ARTS_XML_DATA_PATH)
  set(ARTS_XML_DIR "\"${ARTS_XML_DATA_PATH}\"")
else()
  set(ARTS_XML_DIR "")
endif()
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION ARTS_VERSION)
foreach(file ${PYTHON_SRC_FILES})
  configure_file(${file} ${file} @ONLY)
endforeach()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/controlfiles/)

#
# Install pyarts
#

find_package(PythonInterp)
if (PYTHONINTERP_FOUND)
    add_custom_target(pyarts
      ALL
      COMMAND ${PYTHON_EXECUTABLE} setup.py build
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      DEPENDS arts_api
      COMMENT "Building ARTS python package.")
    set(CONTROLFILE_DIR ${CMAKE_SOURCE_DIR}/controlfiles)
    add_custom_target(python_tests
      ALL
      COMMAND ${PYTHON_EXECUTABLE} bin/arts_convert.py ${CONTROLFILE_DIR} -o controlfiles
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      DEPENDS pyarts
      COMMENT "Converting *.arts controlfiles to Python")
endif()

file(COPY test/reference DESTINATION test/)
file(COPY test/xml/reference DESTINATION test/xml/)
file(COPY test/test_files DESTINATION test/)
