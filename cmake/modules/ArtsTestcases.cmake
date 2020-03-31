macro (ARTS_TEST_RUN_CTLFILE TESTNAME CTLFILE)
  string(REGEX REPLACE "/" "." TESTNAME_LONG ${CTLFILE})
  string(REGEX REPLACE ".arts$" "" TESTNAME_LONG ${TESTNAME_LONG})
  set(TESTNAME_LONG arts.ctlfile.${TESTNAME}.${TESTNAME_LONG})
  set(ARTS arts -r002 -I${CMAKE_CURRENT_SOURCE_DIR})
  if (ARTS_XML_DATA_DIR)
    set(ARTS ${ARTS} -D${ARTS_XML_DATA_DIR})
  endif()
  add_test(
    NAME ${TESTNAME_LONG}
    COMMAND ${ARTS} ${CMAKE_CURRENT_SOURCE_DIR}/${CTLFILE}
  )

  string(REGEX REPLACE ".arts$" ".py" PYTHONCTLFILE ${CTLFILE})
  set(PYTHONCTLFILE "controlfiles/${PYTHONCTLFILE}")
  get_filename_component(CFILESUBDIR ${CTLFILE} DIRECTORY)
  add_test(
    NAME python.${TESTNAME_LONG}
    COMMAND ${PYTHON_EXECUTABLE} ${PYTHONCTLFILE}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python
    )
  if (ARTS_XML_DATA_DIR)
    set_tests_properties(
      python.${TESTNAME_LONG} PROPERTIES
      ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;ARTS_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${CFILESUBDIR}:${CMAKE_SOURCE_DIR}/controlfiles:${ARTS_XML_DATA_DIR}"
      DEPENDS python_tests
    )
  else()
    set_tests_properties(
      python.${TESTNAME_LONG} PROPERTIES
      ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;ARTS_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${CFILESUBDIR}:${CMAKE_SOURCE_DIR}/controlfiles"
      DEPENDS python_tests
    )
  endif()
endmacro (ARTS_TEST_RUN_CTLFILE)

macro (ARTS_TEST_RUN_PYFILE TESTNAME PYFILE)
  string(REGEX REPLACE "/" "." TESTNAME_LONG ${PYFILE})
  string(REGEX REPLACE ".py$" "" TESTNAME_LONG ${TESTNAME_LONG})
  set(TESTNAME_LONG arts.pyarts.${TESTNAME}.${TESTNAME_LONG})

  set(PYFILE "controlfiles/${PYTHONCTLFILE}")
  get_filename_component(CFILESUBDIR ${PYFILE} DIRECTORY)
  add_test(
    NAME ${TESTNAME_LONG}
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${PYFILE}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/controlfiles-python
    )
  if (ARTS_XML_DATA_DIR)
    set_tests_properties(
      ${TESTNAME_LONG} PROPERTIES
      ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;ARTS_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${CFILESUBDIR}:${CMAKE_SOURCE_DIR}/controlfiles:${ARTS_XML_DATA_DIR};ARTS_XML_DATA_DIR=${ARTS_XML_DATA_DIR}"
      DEPENDS python_tests
    )
  else()
    set_tests_properties(
      ${TESTNAME_LONG} PROPERTIES
      ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;ARTS_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/${CFILESUBDIR}:${CMAKE_SOURCE_DIR}/controlfiles"
      DEPENDS python_tests
    )
  endif()
endmacro()

macro (ARTS_TEST_CMDLINE TESTNAME OPTIONS)
  set(ARTS arts)
  add_test(
    NAME arts.cmdline.${TESTNAME}
    COMMAND ${ARTS} ${OPTIONS} ${ARGN}
    )
endmacro (ARTS_TEST_CMDLINE TESTNAME OPTIONS)

macro (ARTS_TEST_CTLFILE_DEPENDS TESTNAME DEPENDNAME)
  set_tests_properties(
    arts.ctlfile.${TESTNAME}
    PROPERTIES DEPENDS arts.ctlfile.${DEPENDNAME}
    )
endmacro (ARTS_TEST_CTLFILE_DEPENDS)

