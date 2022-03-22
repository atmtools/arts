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
# Controlfiles
#

file(GLOB_RECURSE ARTS_CONTROLFILES
  RELATIVE ${ARTS_SOURCE_DIR}
  ${ARTS_SOURCE_DIR}/controlfiles/general/*.arts)
file(GLOB_RECURSE ARTS_CONTROLFILES_INSTRUMENTS
  RELATIVE ${ARTS_SOURCE_DIR}
  ${ARTS_SOURCE_DIR}/controlfiles/instruments/*.arts)
list(APPEND ARTS_CONTROLFILES ${ARTS_CONTROLFILES_INSTRUMENTS})
foreach(file ${ARTS_CONTROLFILES})
  configure_file(${ARTS_SOURCE_DIR}/${file} pyarts/${file})
endforeach()

#
# Install pyarts
#

configure_file(MANIFEST.in MANIFEST.in)

add_custom_target(pyarts
  ALL
  COMMAND ${PYTHON_EXECUTABLE} setup.py build
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts_cpp
  COMMENT "Building ARTS python package.")
set(CONTROLFILE_DIR ${ARTS_SOURCE_DIR}/controlfiles)
add_custom_target(python_tests
  COMMAND PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR} ${PYTHON_EXECUTABLE} bin/arts_convert.py ${CONTROLFILE_DIR} -o controlfiles
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Converting *.arts controlfiles to Python")
#add_test(
#  NAME arts.pyarts.pytest
#  COMMAND ${PYTHON_EXECUTABLE} -m pytest test
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#)
add_custom_target(pyarts_agendatest
  COMMAND ${PYTHON_EXECUTABLE} -m test.workspace.test_agendas
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Run PyARTS Agenda test")

add_custom_target(pyarts_tests
  COMMAND ${PYTHON_EXECUTABLE} -m pytest -v test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS pyarts
  COMMENT "Run PyARTS tests")
add_dependencies(check-deps pyarts_tests)

file(COPY test/plots/reference DESTINATION test/plots/)
file(COPY test/reference DESTINATION test/)
file(COPY test/xml/reference DESTINATION test/xml/)
file(COPY test/test_files DESTINATION test/)

#
# Clean
#
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
  build pyarts/workspace/libarts_api.so)
