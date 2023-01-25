
###########################
### Python Examples ###
if (ENABLE_GUI)
  if (ARTS_CAT_DATA_DIR)
    arts_test_run_pyfile(examples classroom/rotational-spectra/start_gui.py)
    arts_test_run_pyfile(examples classroom/vibrational-spectra/start_gui.py)
  endif()
endif()

arts_test_run_pyfile(examples arts-cat-data/cia/cia.py)
