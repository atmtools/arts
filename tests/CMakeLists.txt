add_subdirectory(testdata)

arts_test_run_pyfile(fast core/spectroscopy/shape_derivatives.py)
arts_test_run_pyfile(fast core/spectroscopy/mpm89.py)
