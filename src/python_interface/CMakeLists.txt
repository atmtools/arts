add_executable(gen_auto_py gen_auto_py.cpp)
target_link_libraries(gen_auto_py PUBLIC ${ALL_ARTS_LIBRARIES})
target_include_directories(gen_auto_py PUBLIC ${ARTS_SOURCE_DIR}/src)

########################################################################################
# Generate automatic python files
add_custom_command(OUTPUT py_auto_workspace.cc
    py_auto_workspace_split_0.cc
    py_auto_workspace_split_1.cc
    py_auto_workspace_split_2.cc
    py_auto_workspace_split_3.cc
    py_auto_workspace_split_4.cc
    py_auto_workspace_split_5.cc
    py_auto_workspace_split_6.cc
    py_auto_workspace_split_7.cc
    py_auto_workspace_split_8.cc
    py_auto_workspace_split_9.cc
    py_auto_workspace_split_10.cc
    py_auto_workspace_split_11.cc
    py_auto_workspace_split_12.cc
    py_auto_workspace_split_13.cc
    py_auto_workspace_split_14.cc
    py_auto_workspace_split_15.cc
    py_auto_workspace_split_16.cc
    py_auto_workspace_split_17.cc
    py_auto_workspace_split_18.cc
    py_auto_workspace_split_19.cc
    py_auto_workspace_split_20.cc
    py_auto_workspace_split_21.cc
    py_auto_workspace_split_22.cc
    py_auto_workspace_split_23.cc
    py_auto_workspace_split_24.cc
    py_auto_workspace_split_25.cc
    py_auto_workspace_split_26.cc
    py_auto_workspace_split_27.cc
    py_auto_workspace_split_28.cc
    py_auto_workspace_split_29.cc
    py_auto_workspace_split_30.cc
    py_auto_workspace_split_31.cc
                   COMMAND gen_auto_py
                   DEPENDS gen_auto_py)
########################################################################################

########### pybind11 class interface
pybind11_add_module(pyarts_cpp
    py_module.cpp
    py_matpack.cpp
    py_basic.cpp
    py_ppath.cpp
    py_griddedfield.cpp
    py_time.cpp
    py_tessem.cpp
    py_rte.cpp
    py_telsem.cpp
    py_species.cpp
    py_sparse.cpp
    py_mcantenna.cpp
    py_scattering.cpp
    py_quantum.cpp
    py_spectroscopy.cpp
    py_jac.cpp
    py_agenda.cpp
    py_workspace_references.cpp
    py_workspace.cpp
    py_auto_workspace.cc
    py_auto_workspace_split_0.cc
    py_auto_workspace_split_1.cc
    py_auto_workspace_split_2.cc
    py_auto_workspace_split_3.cc
    py_auto_workspace_split_4.cc
    py_auto_workspace_split_5.cc
    py_auto_workspace_split_6.cc
    py_auto_workspace_split_7.cc
    py_auto_workspace_split_8.cc
    py_auto_workspace_split_9.cc
    py_auto_workspace_split_10.cc
    py_auto_workspace_split_11.cc
    py_auto_workspace_split_12.cc
    py_auto_workspace_split_13.cc
    py_auto_workspace_split_14.cc
    py_auto_workspace_split_15.cc
    py_auto_workspace_split_16.cc
    py_auto_workspace_split_17.cc
    py_auto_workspace_split_18.cc
    py_auto_workspace_split_19.cc
    py_auto_workspace_split_20.cc
    py_auto_workspace_split_21.cc
    py_auto_workspace_split_22.cc
    py_auto_workspace_split_23.cc
    py_auto_workspace_split_24.cc
    py_auto_workspace_split_25.cc
    py_auto_workspace_split_26.cc
    py_auto_workspace_split_27.cc
    py_auto_workspace_split_28.cc
    py_auto_workspace_split_29.cc
    py_auto_workspace_split_30.cc
    py_auto_workspace_split_31.cc
)

target_link_libraries(pyarts_cpp PUBLIC ${ALL_ARTS_LIBRARIES})
target_include_directories(pyarts_cpp PUBLIC ${ARTS_SOURCE_DIR}/src)
target_include_directories(pyarts_cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pyarts_cpp PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
########################################################################################
