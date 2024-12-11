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

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/pyarts/__init__.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/src/pyarts/__init__.py
  @ONLY
)

configure_file(
  ${ARTS_SOURCE_DIR}/LICENSE.txt
  ${CMAKE_CURRENT_BINARY_DIR}/LICENSE.txt
  COPYONLY
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml.in
  ${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml
  @ONLY
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
  ${CMAKE_CURRENT_BINARY_DIR}/setup.py
  COPYONLY
)

file(GLOB_RECURSE PYTHON_PYARTS_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/src/pyarts/*.py")
add_custom_target(
  pyarts_copy_python_source_files
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different src ${CMAKE_CURRENT_BINARY_DIR}/src
  COMMAND ${CMAKE_COMMAND} -E rm ${CMAKE_CURRENT_BINARY_DIR}/src/pyarts/__init__.py.in
  DEPENDS ${PYTHON_PYARTS_FILES}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copying python/src *.py files"
)

file(GLOB_RECURSE PYTHON_PYARTS_TEST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/test/*")
add_custom_target(
  pyarts_copy_test_python_source_files
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different test ${CMAKE_CURRENT_BINARY_DIR}/test
  DEPENDS ${PYTHON_PYARTS_TEST_FILES}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copying python/test files"
)

#
# Install pyarts
#

configure_file(MANIFEST.in MANIFEST.in)

add_custom_target(pyarts
  ALL
  DEPENDS pyarts_cpp pyarts_cpp_stub pyarts_copy_python_source_files pyarts_copy_test_python_source_files 
  COMMENT "Updating ARTS Python package")
add_custom_target(pyarts-package
  ALL
  COMMAND ${Python_EXECUTABLE} -m build
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts_cpp
  COMMENT "Building ARTS Python package")

add_custom_target(pyarts_tests
  COMMAND ${Python_EXECUTABLE} -m pytest -v test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Running PyARTS pytest")

