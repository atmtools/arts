add_subdirectory(doc)

#
# Copy Python source files.
#

set(ARTS_SRC_DIR ${ARTS_SOURCE_DIR})
if (ARTS_XML_DATA_PATH)
  set(ARTS_XML_DIR "\"${ARTS_XML_DATA_PATH}\"")
else()
  set(ARTS_XML_DIR "")
endif()

file(STRINGS ${ARTS_SOURCE_DIR}/VERSION ARTS_VERSION)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/pyarts/__init__.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/pyarts/__init__.py
  @ONLY
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/setup.py
  @ONLY
)

file(GLOB_RECURSE PYTHON_PYARTS_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/pyarts/*.py")
add_custom_target(
  pyarts_copy_python_source_files
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different pyarts ${CMAKE_CURRENT_BINARY_DIR}/pyarts
  COMMAND ${CMAKE_COMMAND} -E rm ${CMAKE_CURRENT_BINARY_DIR}/pyarts/__init__.py.in
  DEPENDS ${PYTHON_PYARTS_FILES}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copying python/pyarts *.py files"
)

file(GLOB_RECURSE PYTHON_PYARTS_TEST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/test/*")
add_custom_target(
  pyarts_copy_test_python_source_files
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different test ${CMAKE_CURRENT_BINARY_DIR}/test
  DEPENDS ${PYTHON_PYARTS_TEST_FILES}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copying python/test files"
)

file(GLOB_RECURSE PYTHON_PYARTS_BIN_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/bin/*")
add_custom_target(
  pyarts_copy_bin_python_source_files
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different bin ${CMAKE_CURRENT_BINARY_DIR}/bin
  DEPENDS ${PYTHON_PYARTS_BIN_FILES}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copying python/bin files"
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/controlfiles/)

#
# Install pyarts
#

configure_file(MANIFEST.in MANIFEST.in)

add_custom_target(pyarts
  ALL
  DEPENDS pyarts_cpp pyarts_cpp_stub pyarts_copy_python_source_files pyarts_copy_test_python_source_files pyarts_copy_bin_python_source_files
  COMMENT "Updating ARTS python package.")
add_custom_target(pyarts-package
  ALL
  COMMAND ${Python_EXECUTABLE} -m build
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts_cpp
  COMMENT "Building ARTS python package.")

set(CONTROLFILE_DIR ${ARTS_SOURCE_DIR}/tests)
add_custom_target(python_conversion_tests
  COMMAND ${ARTS_PYTHON_INTERPRETER} bin/arts_convert.py ${CONTROLFILE_DIR} -o controlfiles
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Converting ${CONTROLFILE_DIR}/.../*.arts controlfiles to Python")

set(CONTROLFILE_DIR ${ARTS_SOURCE_DIR}/examples)
add_custom_target(python_conversion_examples
  COMMAND ${ARTS_PYTHON_INTERPRETER} bin/arts_convert.py ${CONTROLFILE_DIR} -o controlfiles
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Converting ${CONTROLFILE_DIR}/.../*.arts controlfiles to Python")
add_custom_target(python_conversion
  DEPENDS python_conversion_tests python_conversion_examples
  COMMENT "Converting *.arts controlfiles to Python")

add_custom_target(pyarts_tests
  COMMAND ${Python_EXECUTABLE} -m pytest -v test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Run PyARTS tests")

add_test(NAME pytest
  COMMAND ${Python_EXECUTABLE} -m pytest -v test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(
  pytest PROPERTIES
  DEPENDS pyarts
)
