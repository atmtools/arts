
###########################
### Python Examples ###
if (ENABLE_GUI)
  if (ARTS_CAT_DATA_DIR)
    collect_test_subdir(classroom)
  endif()
endif()

collect_test_subdir(arts-cat-data)
