macro (ARTS_TEST_RUN_EXE TESTNAME EXE)

  set(TESTNAME_LONG cpp.${TESTNAME})

  add_test(
    NAME ${TESTNAME_LONG}
    COMMAND ${EXE}
    )
  set_tests_properties(
    ${TESTNAME_LONG} PROPERTIES
    ENVIRONMENT "ARTS_HEADLESS=1;ARTS_XML_DATA_DIR=${ARTS_XML_DATA_DIR}"
  )
endmacro()

macro (ARTS_TEST_RUN_PYFILE TESTNAME PYFILE)
  if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(DELIM "\;")
  else()
    set(DELIM ":")
  endif()
  string(REGEX REPLACE "/" "." TESTNAME_LONG ${PYFILE})
  string(REGEX REPLACE ".py$" "" TESTNAME_LONG ${TESTNAME_LONG})
  set(TESTNAME_LONG pyarts.${TESTNAME}.${TESTNAME_LONG})

  set(PYFILE "controlfiles/${PYTHONCTLFILE}")
  get_filename_component(CFILESUBDIR ${PYFILE} DIRECTORY)
  add_custom_target(
    make_dir_${TESTNAME}_${TESTNAME_LONG}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/${CFILESUBDIR}
  )
  add_dependencies(check-deps make_dir_${TESTNAME}_${TESTNAME_LONG})

  add_test(
    NAME ${TESTNAME_LONG}
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${PYFILE}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CFILESUBDIR}"
  )
  set_tests_properties(
    ${TESTNAME_LONG} PROPERTIES
    ENVIRONMENT "PYTHONPATH=${ARTS_BINARY_DIR}/python/src;ARTS_HEADLESS=1;ARTS_XML_DATA_DIR=${ARTS_XML_DATA_DIR};ARTS_CAT_DATA_DIR=${ARTS_CAT_DATA_DIR};ARTS_DATA_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${CFILESUBDIR}"
    DEPENDS python_tests
  )
endmacro()

macro (COLLECT_TEST_SUBDIR SUBDIR)
  file(GLOB_RECURSE PYFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${SUBDIR}/*.py)
  get_filename_component(CURRENTDIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  foreach(PYFILE ${PYFILES})
    set(SKIPFILE 0)
    if (ENABLE_ARTS_LGPL AND "${PYFILE}" MATCHES ".nolgpl.")
      set(SKIPFILE 1)
    endif()
    if (NOT SKIPFILE)
      arts_test_run_pyfile(${CURRENTDIR} ${PYFILE})
    else()
      message(STATUS "Build settings don't support ${PYFILE}, skipping.")
    endif()
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
