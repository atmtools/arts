macro (ARTS_TEST_RUN_CTLFILE TESTNAME CTLFILE)
  set(ARTS_TEST_INCLUDE_PATH "${ARTS_SOURCE_DIR}/controlfiles")

  if (ARTS_DATA_DIR)
    string(APPEND ARTS_TEST_INCLUDE_PATH ":${ARTS_DATA_DIR}")
  endif()

  if (ARTS_XML_DATA_DIR)
    string(APPEND ARTS_TEST_INCLUDE_PATH ";ARTS_XML_DATA_DIR=${ARTS_XML_DATA_DIR}")
  else()
    string(APPEND ARTS_TEST_INCLUDE_PATH ";ARTS_XML_DATA_DIR=${ARTS_BINARY_DIR}/tests/testdata/arts-xml-data")
  endif()

  if (ARTS_CAT_DATA_DIR)
    string(APPEND ARTS_TEST_INCLUDE_PATH ";ARTS_CAT_DATA_DIR=${ARTS_CAT_DATA_DIR}")
  else()
    string(APPEND ARTS_TEST_INCLUDE_PATH ";ARTS_CAT_DATA_DIR=${ARTS_BINARY_DIR}/tests/testdata/arts-cat-data")
  endif()

  string(REGEX REPLACE "/" "." TESTNAME_LONG ${CTLFILE})
  string(REGEX REPLACE ".arts$" "" TESTNAME_LONG ${TESTNAME_LONG})
  set(TESTNAME_LONG ctlfile.${TESTNAME}.${TESTNAME_LONG})
  set(ARTS arts -r002 -I${CMAKE_CURRENT_SOURCE_DIR})
  add_test(
    NAME ${TESTNAME_LONG}
    COMMAND ${ARTS} ${CMAKE_CURRENT_SOURCE_DIR}/${CTLFILE}
  )
  set_tests_properties(
    ${TESTNAME_LONG} PROPERTIES
    ENVIRONMENT "ARTS_HEADLESS=1;ARTS_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${CFILESUBDIR}:${ARTS_TEST_INCLUDE_PATH}"
  )

  string(REGEX REPLACE ".arts$" ".py" PYTHONCTLFILE ${CTLFILE})
  string(REGEX REPLACE "^ctlfile" "converted" TESTNAME_LONG ${TESTNAME_LONG})
  set(PYTHONCTLFILE "controlfiles/${PYTHONCTLFILE}")
  get_filename_component(CFILESUBDIR ${CTLFILE} DIRECTORY)
  add_test(
    NAME ${TESTNAME_LONG}
    COMMAND ${Python3_EXECUTABLE} ${PYTHONCTLFILE}
    WORKING_DIRECTORY ${ARTS_BINARY_DIR}/python
    )
  set_tests_properties(
    ${TESTNAME_LONG} PROPERTIES
    ENVIRONMENT "PYTHONPATH=${ARTS_BINARY_DIR}/python;ARTS_HEADLESS=1;ARTS_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${CFILESUBDIR}:${ARTS_TEST_INCLUDE_PATH}"
    DEPENDS python_tests
  )
endmacro (ARTS_TEST_RUN_CTLFILE)

macro (ARTS_TEST_RUN_PYFILE TESTNAME PYFILE)
  set(PYARTS_TEST_INCLUDE_PATH "${ARTS_SOURCE_DIR}/controlfiles")

  if (ARTS_DATA_DIR)
    string(APPEND PYARTS_TEST_INCLUDE_PATH ":${ARTS_DATA_DIR}")
  endif()

  if (ARTS_XML_DATA_DIR)
    string(APPEND PYARTS_TEST_INCLUDE_PATH ";ARTS_XML_DATA_DIR=${ARTS_XML_DATA_DIR}")
  else()
    string(APPEND PYARTS_TEST_INCLUDE_PATH ";ARTS_XML_DATA_DIR=${ARTS_BINARY_DIR}/tests/testdata/arts-xml-data")
  endif()

  if (ARTS_CAT_DATA_DIR)
    string(APPEND PYARTS_TEST_INCLUDE_PATH ";ARTS_CAT_DATA_DIR=${ARTS_CAT_DATA_DIR}")
  else()
    string(APPEND PYARTS_TEST_INCLUDE_PATH ";ARTS_CAT_DATA_DIR=${ARTS_BINARY_DIR}/tests/testdata/arts-cat-data")
  endif()

  string(REGEX REPLACE "/" "." TESTNAME_LONG ${PYFILE})
  string(REGEX REPLACE ".py$" "" TESTNAME_LONG ${TESTNAME_LONG})
  set(TESTNAME_LONG pyarts.${TESTNAME}.${TESTNAME_LONG})

  set(PYFILE "controlfiles/${PYTHONCTLFILE}")
  get_filename_component(CFILESUBDIR ${PYFILE} DIRECTORY)
  add_test(
    NAME ${TESTNAME_LONG}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${PYFILE}
    WORKING_DIRECTORY ${ARTS_BINARY_DIR}/controlfiles-python
    )
  set_tests_properties(
    ${TESTNAME_LONG} PROPERTIES
    ENVIRONMENT "PYTHONPATH=${ARTS_BINARY_DIR}/python;ARTS_HEADLESS=1;ARTS_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${CFILESUBDIR}:${PYARTS_TEST_INCLUDE_PATH}"
    DEPENDS python_tests
  )
endmacro()

macro (ARTS_TEST_CMDLINE TESTNAME OPTIONS)
  set(ARTS arts)
  add_test(
    NAME cmdline.${TESTNAME}
    COMMAND ${ARTS} ${OPTIONS} ${ARGN}
    )
endmacro (ARTS_TEST_CMDLINE TESTNAME OPTIONS)

macro (ARTS_TEST_CTLFILE_DEPENDS TESTNAME DEPENDNAME)
  set_tests_properties(
    ctlfile.${TESTNAME}
    PROPERTIES DEPENDS ctlfile.${DEPENDNAME}
    )
endmacro (ARTS_TEST_CTLFILE_DEPENDS)

macro (COLLECT_TEST_SUBDIR SUBDIR)
  file(GLOB_RECURSE PYFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${SUBDIR}/*.py)
  get_filename_component(CURRENTDIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  foreach(PYFILE ${PYFILES})
    arts_test_run_pyfile(${CURRENTDIR} ${PYFILE})
  endforeach()

  file(GLOB_RECURSE ARTSFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${SUBDIR}/*.arts)
  foreach(ARTSFILE ${ARTSFILES})
    arts_test_run_ctlfile(${CURRENTDIR} ${ARTSFILE})
  endforeach()

endmacro ()

macro (SETUP_ARTS_CHECKS)
  set(CTEST_ARGS ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS})

  # This will be removed in ARTS-3, only here to work around currently broken controlfiles
  set(EXCLUDED_FROM_CHECK '\(\\.nocheck\\.|\\.slow\\.|\\.xmldata\\.|\\.planettoolbox\\.\)')

  add_custom_target(check-deps)

  add_custom_target(check
    COMMAND ${CTEST_ARGS}
    -R '\(^ctlfile|^pytest|^pyarts|^doc|^cmdline|^cpp\)' -E ${EXCLUDED_FROM_CHECK}
    DEPENDS check-deps pyarts)

  add_custom_target(check-pyarts
    COMMAND ${CTEST_ARGS}
    -R '\(^pyarts\)' -E ${EXCLUDED_FROM_CHECK}
    DEPENDS check-deps pyarts)

  add_custom_target(check-controlfiles
    COMMAND ${CTEST_ARGS}
    -R '\(^ctlfile\)' -E ${EXCLUDED_FROM_CHECK}
    DEPENDS check-deps)

  add_custom_target(check-conversion
    COMMAND ${CTEST_ARGS}
    -R '\(^converted\)' -E ${EXCLUDED_FROM_CHECK}
    DEPENDS check-deps python_conversion)

  add_custom_target(check-examples
    COMMAND ${CTEST_ARGS}
    -R '\(\\.examples\\.\)' -E '\(^converted\)'
    DEPENDS check-deps pyarts)

  add_custom_target(check-tests
    COMMAND ${CTEST_ARGS}
    -R '\(\\.tests\\.\)' -E '\(^converted\)'
    DEPENDS check-deps pyarts)

  add_custom_target(check-doc
    COMMAND ${CTEST_ARGS}
    -R '\(^doc\)'
    DEPENDS check-deps)

  add_custom_target(check-pytest DEPENDS pyarts_tests)

  # Everything below this point will go away in ARTS-3
  add_custom_target(check-xmldata
    COMMAND ${CTEST_ARGS}
    -R '\(\\.xmldata\\.\)'
    DEPENDS check-deps pyarts)
  add_dependencies(check-xmldata mkdir-arts-results)

  add_custom_target(check-planettoolbox
    COMMAND ${CTEST_ARGS}
    -R '\(\\.planettoolbox\\.\)'
    DEPENDS check-deps)

  add_custom_target(check-all
    COMMAND ${CTEST_ARGS}
    -E '\\.nocheck\\.'
    DEPENDS check-deps python_conversion)
  add_dependencies(check-all mkdir-arts-results)

  add_custom_target(check-all-controlfiles
    COMMAND ${CTEST_ARGS}
    -R '\(^ctlfile\)' -E '\\.nocheck\\.'
    DEPENDS check-deps)

endmacro ()