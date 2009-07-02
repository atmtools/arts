
add_custom_command(
  OUTPUT PrepareTests
  COMMAND cmake -E
    ARGS copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/artsXML.py
                           ${CMAKE_CURRENT_BINARY_DIR}/artsXML.py
  COMMAND cmake -E
    ARGS copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/testall.py
                           ${CMAKE_CURRENT_BINARY_DIR}/testall.py
)

arts_add_prepare_testcases (TARGETS
                            Abs
                            AMSU
                            ClearSky
                            DISORT
                            DOIT
                            MonteCarlo
                            OdinSMR
)
set (PREPARETESTTARGETS PrepareTests ${TARGETS})

arts_add_testcases (TARGETS ""
                    MonteCarloDataPrepare
                    OdinSMR
                    DOIT
                    ClearSky
                    Abs
                    AMSUB_fast
)
set(TESTFASTTARGETS ${TARGETS})

arts_add_testcases (TARGETS TestMonteCarloDataPrepare
                    MonteCarloGeneral
                    MonteCarloGeneralGaussian
                    RteCalcMC
)
set(TESTFASTTARGETS ${TESTFASTTARGETS} ${TARGETS})

arts_add_testcases (TARGETS TestAMSUB_fast
                    AMSUB
)
set(TESTTARGETS ${TARGETS})

add_custom_target (check-fast
                   DEPENDS
                     ${PREPARETESTTARGETS}
                     ${TESTFASTTARGETS}
)

add_custom_target (check
                   DEPENDS
                     ${PREPARETESTTARGETS}
                     ${TESTFASTTARGETS}
                     ${TESTTARGETS}
)

