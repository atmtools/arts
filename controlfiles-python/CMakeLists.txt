
###########################
### Python Test Example ###
arts_test_run_pyfile(fast example/TestExample.py)

if (ARTS_XML_DATA_DIR)
  arts_test_run_pyfile(fast classes/TestClasses.py)
endif()
