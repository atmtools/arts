find_package(PythonInterp)

if (PYTHONINTERP_FOUND)

add_custom_command(
  OUTPUT PrepareTestScripts
  COMMAND ${CMAKE_COMMAND} -E
    ARGS copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/artsXML.py
                           ${CMAKE_CURRENT_BINARY_DIR}/artsXML.py
  COMMAND ${CMAKE_COMMAND} -E
    ARGS copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/testall.py
                           ${CMAKE_CURRENT_BINARY_DIR}/testall.py
)

set(TESTSLOWTARGETS "")
set(TESTFASTTARGETS "")

arts_add_prepare_testcases (PREPTARGETS OdinSMR)
arts_add_testcases (TARGETS "${PREPTARGETS}" OdinSMR)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS MonteCarlo)
arts_add_testcases (TARGETS "${PREPTARGETS}" MonteCarloDataPrepare)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

set (PREPTARGETS ${PREPTARGETS} TestMonteCarloDataPrepare)
arts_add_testcases (TARGETS "${PREPTARGETS}"
                    MonteCarloGeneral
                    MonteCarloGeneralGaussian
                    RteCalcMC
)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS DOIT)
arts_add_testcases (TARGETS "${PREPTARGETS}" DOIT)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS ClearSky)
arts_add_testcases (TARGETS "${PREPTARGETS}" ClearSky)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS Abs)
arts_add_testcases (TARGETS "${PREPTARGETS}" Abs)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS GroundBased)
arts_add_testcases (TARGETS "${PREPTARGETS}" GroundBased)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS AMSU)
arts_add_testcases (TARGETS "${PREPTARGETS}" AMSUB_fast)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

set (PREPTARGETS ${PREPTARGETS} TestAMSUB_fast)
arts_add_testcases (TARGETS "${PREPTARGETS}" AMSUB)
set(TESTSLOWTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS MHS)
arts_add_testcases (TARGETS "${PREPTARGETS}" MHS_fast)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

set (PREPTARGETS ${PREPTARGETS} TestMHS_fast)
arts_add_testcases (TARGETS "${PREPTARGETS}" MHS)
set(TESTSLOWTARGETS ${TESTFASTTARGETS} ${TARGETS})

add_custom_target (check-fast
                   DEPENDS
                     ${TESTFASTTARGETS}
)

add_custom_target (check
                   DEPENDS
                     ${TESTFASTTARGETS}
                     ${TESTSLOWTARGETS}
)

else (PYTHONINTERP_FOUND)
  message(STATUS "Test cases disabled due to missing Python interpreter")
endif (PYTHONINTERP_FOUND)

