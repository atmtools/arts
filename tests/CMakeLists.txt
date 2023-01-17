add_subdirectory(testdata)

# core/geodetic
arts_test_run_ctlfile(fast core/geodetic/TestBasicGeodeticAccuracy.arts)
arts_test_run_ctlfile(fast core/geodetic/TestIntersections.arts)

# core/propmat
arts_test_run_pyfile(fast core/propmat/rescale_species.py)

# core/ppath
arts_test_run_ctlfile(fast core/ppath/TestPpathAddGridCrossings.arts)
arts_test_run_ctlfile(fast core/ppath/TestPpathGeometric.arts)
arts_test_run_ctlfile(fast core/ppath/TestPpathRefracted.arts)

# core/spectroscopy
arts_test_run_pyfile(fast core/spectroscopy/shape_derivatives.py)
arts_test_run_pyfile(fast core/spectroscopy/vmr_derivatives.py)
arts_test_run_pyfile(fast core/spectroscopy/ell07.py)
arts_test_run_pyfile(fast core/spectroscopy/mpm89.py)
arts_test_run_pyfile(fast core/spectroscopy/ckdmt400.py)
arts_test_run_pyfile(fast core/spectroscopy/ckdmt252.py)
arts_test_run_pyfile(fast core/spectroscopy/ckdmt100.py)
arts_test_run_pyfile(fast core/spectroscopy/pwr98.py)
arts_test_run_pyfile(fast core/spectroscopy/tre05.py)
arts_test_run_pyfile(fast core/spectroscopy/standard.py)
arts_test_run_pyfile(fast core/star/solar_irradiance.py)

# core/star

# core/surface
arts_test_run_ctlfile(fast core/surface/TestSpecularLos.arts)
arts_test_run_ctlfile(fast core/surface/TestSurfRtpropAgenda.arts)
