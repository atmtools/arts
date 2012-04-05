macro (ARTS_TEST_RUN_CTLFILE TESTNAME CTLFILE)
  set(ARTS arts -r002)
  add_test(
    NAME arts.ctlfile.${TESTNAME}
    COMMAND ${ARTS} ${CMAKE_CURRENT_SOURCE_DIR}/${CTLFILE}
    )
endmacro (ARTS_TEST_RUN_CTLFILE)

macro (ARTS_TEST_CTLFILE_DEPENDS TESTNAME DEPENDNAME)
  set_tests_properties(
    arts.ctlfile.${TESTNAME}
    PROPERTIES DEPENDS arts.ctlfile.${DEPENDNAME}
    )
endmacro (ARTS_TEST_CTLFILE_DEPENDS)

