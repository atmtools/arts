if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(DELIM "\;")
else()
  set(DELIM ":")
endif()

macro (ARTS_TEST_RUN_EXE TESTNAME EXE)
  if (ARTS_DATA_DIR)
    string(APPEND CPP_TEST_INCLUDE_PATH "${DELIM}${ARTS_DATA_DIR}")
  endif()

  if (ARTS_XML_DATA_DIR)
    string(APPEND CPP_TEST_INCLUDE_PATH ";ARTS_XML_DATA_DIR=${ARTS_XML_DATA_DIR}")
  else()
    string(APPEND CPP_TEST_INCLUDE_PATH ";ARTS_XML_DATA_DIR=${ARTS_BINARY_DIR}/tests/testdata/arts-xml-data")
  endif()

  if (ARTS_CAT_DATA_DIR)
    string(APPEND CPP_TEST_INCLUDE_PATH ";ARTS_CAT_DATA_DIR=${ARTS_CAT_DATA_DIR}")
  else()
    string(APPEND CPP_TEST_INCLUDE_PATH ";ARTS_CAT_DATA_DIR=${ARTS_BINARY_DIR}/tests/testdata/arts-cat-data")
  endif()

  set(TESTNAME_LONG cpp.${TESTNAME})

  add_test(
    NAME ${TESTNAME_LONG}
    COMMAND ${EXE}
    )
  set_tests_properties(
    ${TESTNAME_LONG} PROPERTIES
    ENVIRONMENT "ARTS_HEADLESS=1;ARTS_INCLUDE_PATH=${CPP_TEST_INCLUDE_PATH}"
  )
endmacro()

macro (ARTS_TEST_RUN_PYFILE TESTNAME PYFILE)
  if (ARTS_DATA_DIR)
    string(APPEND PYARTS_TEST_INCLUDE_PATH "${DELIM}${ARTS_DATA_DIR}")
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
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${PYFILE}
    )
  set_tests_properties(
    ${TESTNAME_LONG} PROPERTIES
    ENVIRONMENT "PYTHONPATH=${ARTS_BINARY_DIR}/python/src;ARTS_HEADLESS=1;ARTS_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${CFILESUBDIR}${DELIM}${PYARTS_TEST_INCLUDE_PATH}"
    DEPENDS python_tests
  )
endmacro()

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
  set(CTEST_ARGS ${CMAKE_CTEST_COMMAND} ${ARTS_CTEST_USER_OPTIONS} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS})

  add_custom_target(check-deps)

  add_custom_target(check
    COMMAND ${CTEST_ARGS}
    -R \"\(^ctlfile|^pytest|^pyarts|^doc||^cpp\)\"
    DEPENDS check-deps pyarts)

  add_custom_target(check-pyarts
    COMMAND ${CTEST_ARGS}
    -R \"\(^pyarts\)\"
    DEPENDS check-deps pyarts)

  add_custom_target(check-examples
    COMMAND ${CTEST_ARGS}
    -R \"\(\\.examples\\.\)\"
    DEPENDS check-deps pyarts)

  add_custom_target(check-tests
    COMMAND ${CTEST_ARGS}
    -R \"\(\\.tests\\.\)\"
    DEPENDS check-deps pyarts)

  add_custom_target(check-doc
    COMMAND ${CTEST_ARGS}
    -R \"\(^doc\)\"
    DEPENDS check-deps)

  add_custom_target(check-pytest DEPENDS pyarts_tests)

endmacro ()
