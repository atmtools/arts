if (ARTS_XML_DATA_DIR)
  arts_test_run_pyfile(slow artscomponents/surfaceblackbody/TestSurfaceBlackbody.py)
  arts_test_run_pyfile(fast artscomponents/linemixing/TestMolOxyAdaptation.py)
endif()

arts_test_run_pyfile(fast artscomponents/cutoff/ycalc.py)
arts_test_run_pyfile(fast artscomponents/zeeman/ycalc.py)
arts_test_run_pyfile(fast artscomponents/predefined/TestCKDMT350.py)
