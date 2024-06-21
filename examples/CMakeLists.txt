
###########################
### Python Examples ###
if (ENABLE_GUI)
  if (ARTS_CAT_DATA_DIR)
    collect_test_subdir(classroom)
  endif()
endif()

collect_test_subdir(arts-catalogue-data)
collect_test_subdir(getting-started)
collect_test_subdir(recipes)
