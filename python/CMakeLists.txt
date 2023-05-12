add_subdirectory(doc)

#
# Copy Python source files.
#

file(GLOB_RECURSE PYTHON_SRC_FILES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
set(ARTS_SRC_DIR ${ARTS_SOURCE_DIR})
if (ARTS_XML_DATA_PATH)
  set(ARTS_XML_DIR "\"${ARTS_XML_DATA_PATH}\"")
else()
  set(ARTS_XML_DIR "")
endif()
file(STRINGS ${ARTS_SOURCE_DIR}/VERSION ARTS_VERSION)
foreach(file ${PYTHON_SRC_FILES})
  configure_file(${file} ${file} @ONLY)
endforeach()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/controlfiles/)

#
# Install pyarts
#

configure_file(MANIFEST.in MANIFEST.in)

add_custom_target(pyarts
  ALL
  COMMAND ${Python3_EXECUTABLE} setup.py build
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts_cpp
  COMMENT "Building ARTS python package.")
set(CONTROLFILE_DIR ${ARTS_SOURCE_DIR}/controlfiles)
add_custom_target(python_conversion_controlfiles
  COMMAND ${ARTS_PYTHON_INTERPRETER} bin/arts_convert.py --nofail ${CONTROLFILE_DIR} -o controlfiles
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Converting *.arts controlfiles to Python")
set(CONTROLFILE_DIR ${ARTS_SOURCE_DIR}/tests)
add_custom_target(python_conversion_tests
  COMMAND ${ARTS_PYTHON_INTERPRETER} bin/arts_convert.py ${CONTROLFILE_DIR} -o controlfiles
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Converting *.arts controlfiles to Python")
set(CONTROLFILE_DIR ${ARTS_SOURCE_DIR}/examples)
add_custom_target(python_conversion_examples
  COMMAND ${ARTS_PYTHON_INTERPRETER} bin/arts_convert.py ${CONTROLFILE_DIR} -o controlfiles
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Converting *.arts controlfiles to Python")
add_custom_target(python_conversion
  DEPENDS python_conversion_controlfiles python_conversion_tests python_conversion_examples
  COMMENT "Converting *.arts controlfiles to Python")

add_custom_target(pyarts_tests
  COMMAND ${Python3_EXECUTABLE} -m pytest -v test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Run PyARTS tests")

add_test(NAME pytest
  COMMAND ${Python3_EXECUTABLE} -m pytest -v test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(
  pytest PROPERTIES
  DEPENDS pyarts
)
file(COPY test/plots/reference DESTINATION test/plots/)
file(COPY test/reference DESTINATION test/)
file(COPY test/xml/reference DESTINATION test/xml/)
file(COPY test/test_files DESTINATION test/)

#
# Clean
#
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
  build pyarts/workspace/libarts_api.so)
