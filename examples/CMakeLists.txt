
###########################
### Python Examples ###
if (ENABLE_GUI)
  if (ARTS_XML_DATA_DIR)
    arts_test_run_pyfile(examples classroom/rotational-spectra/start_gui.py)
  endif()
endif()
