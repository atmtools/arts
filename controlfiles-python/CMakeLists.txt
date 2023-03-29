
###########################
### Python Test Example ###

arts_test_run_pyfile(fast raw/calib.py)
arts_test_run_pyfile(fast raw/corr.py)

if (ARTS_XML_DATA_DIR)
  arts_test_run_pyfile(slow artscomponents/surfaceblackbody/TestSurfaceBlackbody.py)
endif()

# Data downloaded as part of the automatic downloader
arts_test_run_pyfile(fast artscomponents/linemixing/TestMolOxyAdaptation.py)

arts_test_run_pyfile(fast artscomponents/cloudbox/Testcloudbox_fieldInterp2Azimuth.py)
arts_test_run_pyfile(fast artscomponents/cutoff/ycalc.py)
arts_test_run_pyfile(fast artscomponents/disort/TestDisortCalcIrradiance.py)
arts_test_run_pyfile(fast artscomponents/gasscattering/TestClearsky_StarGasScatteringIsotropic.py)
arts_test_run_pyfile(fast artscomponents/gasscattering/TestClearsky_StarGasScatteringRayleigh.py)
arts_test_run_pyfile(fast artscomponents/predefined/TestCKDMT350.py)
arts_test_run_pyfile(fast artscomponents/refractiveindex/Test_refr_index_waterVisibleNIR.py)
arts_test_run_pyfile(fast artscomponents/surface/Test_iySurfaceFlatReflectivity.py)
arts_test_run_pyfile(fast artscomponents/surface/Test_iySurfaceFlatRefractiveIndex.py)
arts_test_run_pyfile(fast artscomponents/surface/Test_iySurfaceLambertian.py)
arts_test_run_pyfile(fast artscomponents/zeeman/ycalc.py)
