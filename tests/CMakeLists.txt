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

arts_add_prepare_testcases (PREPTARGETS atmosphere_data OdinSMR)
arts_add_testcases (TARGETS "${PREPTARGETS}" OdinSMR)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS atmosphere_data MonteCarlo)
arts_add_testcases (TARGETS "${PREPTARGETS}" MonteCarloDataPrepare)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

set (PREPTARGETS ${PREPTARGETS} TestMonteCarloDataPrepare)
arts_add_testcases (TARGETS "${PREPTARGETS}"
                    MonteCarloGeneral
                    MonteCarloGeneralGaussian
                    RteCalcMC
)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS atmosphere_data DOIT)
arts_add_testcases (TARGETS "${PREPTARGETS}" DOIT)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS DOITBatch/scatData DOITBatch)
arts_add_testcases (TARGETS "${PREPTARGETS}" DOITBatch)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS atmosphere_data ClearSky)
arts_add_testcases (TARGETS "${PREPTARGETS}" ClearSky)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS atmosphere_data Abs)
arts_add_testcases (TARGETS "${PREPTARGETS}" Abs)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS atmosphere_data GroundBased)
arts_add_testcases (TARGETS "${PREPTARGETS}" GroundBased)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS AMSU)
arts_add_testcases (TARGETS "${PREPTARGETS}" AMSUB_fast)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

set (PREPTARGETS ${PREPTARGETS} TestAMSUB_fast)
arts_add_testcases (TARGETS "${PREPTARGETS}" AMSUB)
set(TESTSLOWTARGETS ${TESTSLOWTARGETS} ${TARGETS})

arts_add_prepare_testcases (PREPTARGETS MHS)
arts_add_testcases (TARGETS "${PREPTARGETS}" MHS_fast)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

set (PREPTARGETS ${PREPTARGETS} TestMHS_fast)
arts_add_testcases (TARGETS "${PREPTARGETS}" MHS)
set(TESTSLOWTARGETS ${TESTSLOWTARGETS} ${TARGETS})

if (NOT "${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  get_directory_property (TESTCLEANFILES ADDITIONAL_MAKE_CLEAN_FILES)
  set (TESTCLEANFILES ${TESTCLEANFILES}
       artsXML.py artsXML.pyc testall.py atmosphere_data
       AMSU Abs ClearSky DOIT DOITBatch GroundBased MHS MonteCarlo OdinSMR)
  set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                            "${TESTCLEANFILES}")
endif (NOT "${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")

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

