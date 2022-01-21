########### pybind11 class interface
pybind11_add_module(pyarts_cpp py_module.cpp py_matpack.cpp py_basic.cpp)
target_link_libraries(pyarts_cpp PUBLIC ${ALL_ARTS_LIBRARIES})
target_include_directories(pyarts_cpp PUBLIC ${ARTS_SOURCE_DIR}/src)
########################################################################################
