add_subdirectory(testdata)

arts_test_run_pyfile(fast core/spectroscopy/shape_derivatives.py)
arts_test_run_pyfile(fast core/spectroscopy/ell07.py)
arts_test_run_pyfile(fast core/spectroscopy/mpm89.py)
arts_test_run_pyfile(fast core/spectroscopy/ckdmt252.py)
arts_test_run_pyfile(fast core/spectroscopy/ckdmt100.py)
arts_test_run_pyfile(fast core/spectroscopy/pwr98.py)
