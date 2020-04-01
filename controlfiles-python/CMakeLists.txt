
###########################
### Python Test Example ###
arts_test_run_pyfile(fast example/TestExample.py)
arts_test_run_pyfile(fast classes/TestEnergyLevelMap.py)
arts_test_run_pyfile(fast classes/TestBasicTypes.py)
arts_test_run_pyfile(fast classes/TestAgenda.py)
arts_test_run_pyfile(fast classes/TestGriddedField.py)
arts_test_run_pyfile(fast classes/TestGasAbsLookup.py)
arts_test_run_pyfile(fast classes/TestGridPos.py)
arts_test_run_pyfile(fast classes/TestMatpackTypes.py)
arts_test_run_pyfile(fast classes/TestMCAntenna.py)
arts_test_run_pyfile(fast classes/TestPpath.py)
arts_test_run_pyfile(fast classes/TestPropagationTypes.py)
arts_test_run_pyfile(fast classes/TestQuantum.py)
arts_test_run_pyfile(fast classes/TestRetrievalQuantity.py)
arts_test_run_pyfile(fast classes/TestRational.py)

if (ARTS_XML_DATA_DIR)
  arts_test_run_pyfile(fast classes/TestAbsorptionLines.py)
  arts_test_run_pyfile(fast classes/TestCIARecord.py)
endif()
