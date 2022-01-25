add_executable(gen_auto_py gen_auto_py.cpp)
target_link_libraries(gen_auto_py PUBLIC ${ALL_ARTS_LIBRARIES})
target_include_directories(gen_auto_py PUBLIC ${ARTS_SOURCE_DIR}/src)

########### pybind11 class interface
pybind11_add_module(pyarts_cpp py_module.cpp py_matpack.cpp py_basic.cpp py_ppath.cpp py_griddedfield.cpp py_time.cpp py_tessem.cpp py_rte.cpp py_telsem.cpp py_species.cpp)
target_link_libraries(pyarts_cpp PUBLIC ${ALL_ARTS_LIBRARIES})
target_include_directories(pyarts_cpp PUBLIC ${ARTS_SOURCE_DIR}/src)
########################################################################################
