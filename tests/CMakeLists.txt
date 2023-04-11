add_subdirectory(testdata)

# CORE SPECTROSCOPY
arts_test_run_pyfile(fast core/spectroscopy/shape_derivatives.py)
arts_test_run_pyfile(fast core/spectroscopy/vmr_derivatives.py)
arts_test_run_pyfile(fast core/spectroscopy/ell07.py)
arts_test_run_pyfile(fast core/spectroscopy/mpm89.py)
arts_test_run_pyfile(fast core/spectroscopy/ckdmt400.py)
arts_test_run_pyfile(fast core/spectroscopy/ckdmt252.py)
arts_test_run_pyfile(fast core/spectroscopy/ckdmt100.py)
arts_test_run_pyfile(fast core/spectroscopy/pwr98.py)
arts_test_run_pyfile(fast core/spectroscopy/pwr2021.py)
arts_test_run_pyfile(fast core/spectroscopy/pwr2022.py)
arts_test_run_pyfile(fast core/spectroscopy/tre05.py)
arts_test_run_pyfile(fast core/spectroscopy/standard.py)
arts_test_run_pyfile(fast core/star/solar_irradiance.py)

# CORE PROPMAT
arts_test_run_pyfile(fast core/propmat/rescale_species.py)