add_subdirectory(testdata)

collect_test_subdir(core)
collect_test_subdir(hitran)
collect_test_subdir(python)

if (NOT ENABLE_ARTS_LGPL)
  collect_test_subdir(aer)
endif()
