add_subdirectory(testdata)

collect_test_subdir(core)

if (NOT ENABLE_ARTS_LGPL)
  collect_test_subdir(aer)
  collect_test_subdir(disort)
endif()
