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

arts_add_prepare_testcases (PREPTARGETS Wfuns)
arts_add_testcases (TARGETS "${PREPTARGETS}" Wfuns)
set(TESTSLOWTARGETS ${TESTSLOWTARGETS} ${TARGETS})

# GH 2011-05-05: comment out for now, as this depends on PyARTS
# To be discussed (PyARTS makes making tests a lot easier, but would
# create a depedency on Python+PyARTS+numpy+scipy)
#
#set (PREPTARGETS ${PREPTARGETS} Testatm_fields_compactAddSpecies)
#arts_add_testcases (TARGETS "${PREPTARGETS}" atm_fields_compactAddSpecies)
#set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

if (NOT "${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  get_directory_property (TESTCLEANFILES ADDITIONAL_MAKE_CLEAN_FILES)
  set (TESTCLEANFILES ${TESTCLEANFILES}
       artsXML.py artsXML.pyc testall.py atmosphere_data
       AMSU Abs ClearSky DOIT DOITBatch GroundBased MHS MonteCarlo OdinSMR Wfuns)
  set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                            "${TESTCLEANFILES}")
endif (NOT "${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")

add_custom_target (check-fast DEPENDS ${TESTFASTTARGETS})

add_custom_target (check DEPENDS ${TESTFASTTARGETS} ${TESTSLOWTARGETS})

else (PYTHONINTERP_FOUND)
  message(STATUS "Test cases disabled due to missing Python interpreter")
endif (PYTHONINTERP_FOUND)

############ tests ##########

add_custom_command(
  OUTPUT clean-tests-cmd
  COMMAND rm -f "Test*.xml" "Test*.rep" "Test*.bin"
  )
add_custom_target (clean-tests DEPENDS clean-tests-cmd)

arts_test_run_ctlfile(fast.abs            Abs/TestAbs.arts)
arts_test_run_ctlfile(fast.ppath1d        Ppath/TestPpath1D.arts)
arts_test_run_ctlfile(fast.ppath2d        Ppath/TestPpath2D.arts)
arts_test_run_ctlfile(fast.ppath3d        Ppath/TestPpath3D.arts)
arts_test_run_ctlfile(fast.pencilbeam     PencilBeam/TestPencilBeam.arts)
arts_test_run_ctlfile(fast.clearsky       ClearSky/TestClearSky.arts)
arts_test_run_ctlfile(fast.groundbased    GroundBased/TestGbased.arts)
arts_test_run_ctlfile(fast.transmission   Transmission/TestTransmission.arts)
arts_test_run_ctlfile(fast.radiolink      RadioLink/TestRadioLink.arts)
arts_test_run_ctlfile(fast.amsub          AMSU/TestAMSUB_fast.arts)
arts_test_run_ctlfile(fast.doit           DOIT/TestDOIT.arts)

arts_test_run_ctlfile(fast.montecarlo.dataprepare
                      MonteCarlo/TestMonteCarloDataPrepare.arts)

arts_test_run_ctlfile(slow.montecarlo.general
                      MonteCarlo/TestMonteCarloGeneral.arts)
arts_test_ctlfile_depends(slow.montecarlo.general
                          fast.montecarlo.dataprepare)

arts_test_run_ctlfile(slow.montecarlo.generalgaussian
                      MonteCarlo/TestMonteCarloGeneralGaussian.arts)
arts_test_ctlfile_depends(slow.montecarlo.generalgaussian
                          fast.montecarlo.dataprepare)

arts_test_run_ctlfile(fast.montecarlo.rtecalc MonteCarlo/TestRteCalcMC.arts)
arts_test_ctlfile_depends(fast.montecarlo.rtecalc fast.montecarlo.dataprepare)

arts_test_run_ctlfile(slow.amsub AMSU/TestAMSUB.arts)

arts_test_run_ctlfile(slow.doitbatch DOITBatch/TestDOITBatch.arts)

arts_test_run_ctlfile(fast.mhs MHS/TestMHS_fast.arts)
arts_test_run_ctlfile(slow.mhs MHS/TestMHS.arts)

