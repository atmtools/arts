########################################################################################
# Generate list of files to be generated
if (NOT NUM_PYARTS_WSM)
    set(NUM_PYARTS_WSM 5)
endif()
message(STATUS "Using ${NUM_PYARTS_WSM} python interface automatic workspace method files (-DNUM_PYARTS_WSM=N to change)")

if (NOT NUM_PYARTS_WSV)
    set(NUM_PYARTS_WSV 3)
endif()
message(STATUS "Using ${NUM_PYARTS_WSV} python interface automatic workspace variable files (-DNUM_PYARTS_WSV=N to change)")

if (NOT NUM_PYARTS_WSC)
    set(NUM_PYARTS_WSC 3)
endif()
message(STATUS "Using ${NUM_PYARTS_WSC} python interface automatic workspace create method files (-DNUM_PYARTS_WSC=N to change)")

if (NOT NUM_PYARTS_WSG)
    set(NUM_PYARTS_WSG 3)
endif()
message(STATUS "Using ${NUM_PYARTS_WSG} python interface automatic workspace generic method files (-DNUM_PYARTS_WSG=N to change)")

set(PYARTS_WORKSPACE_METHODS "")
list(LENGTH PYARTS_WORKSPACE_METHODS N)
while (N LESS NUM_PYARTS_WSM)
    list(APPEND PYARTS_WORKSPACE_METHODS "py_auto_workspace_split_methods_${N}.cc")
    list(LENGTH PYARTS_WORKSPACE_METHODS N)
endwhile()

set(PYARTS_WORKSPACE_VARS "")
list(LENGTH PYARTS_WORKSPACE_VARS N)
while (N LESS NUM_PYARTS_WSV)
    list(APPEND PYARTS_WORKSPACE_VARS "py_auto_workspace_split_vars_${N}.cc")
    list(LENGTH PYARTS_WORKSPACE_VARS N)
endwhile()

set(PYARTS_WORKSPACE_CREATE "")
list(LENGTH PYARTS_WORKSPACE_CREATE N)
while (N LESS NUM_PYARTS_WSC)
    list(APPEND PYARTS_WORKSPACE_CREATE "py_auto_workspace_split_create_${N}.cc")
    list(LENGTH PYARTS_WORKSPACE_CREATE N)
endwhile()

set(PYARTS_WORKSPACE_GEN "")
list(LENGTH PYARTS_WORKSPACE_GEN N)
while (N LESS NUM_PYARTS_WSG)
    list(APPEND PYARTS_WORKSPACE_GEN "py_auto_workspace_split_generic_${N}.cc")
    list(LENGTH PYARTS_WORKSPACE_GEN N)
endwhile()

########################################################################################
# Build the generator
add_executable(gen_auto_py gen_auto_py.cpp)
target_link_libraries(gen_auto_py PUBLIC methods)
target_include_directories(gen_auto_py PUBLIC ${ARTS_SOURCE_DIR}/src)

# Generate the actual files with this command
add_custom_command(OUTPUT py_auto_interface.h
                          py_auto_interface.cc
                          py_auto_workspace.cc
                          ${PYARTS_WORKSPACE_METHODS}
                          ${PYARTS_WORKSPACE_VARS}
                          ${PYARTS_WORKSPACE_GEN}
                          ${PYARTS_WORKSPACE_CREATE}
                   COMMAND gen_auto_py ${NUM_PYARTS_WSM} ${NUM_PYARTS_WSV} ${NUM_PYARTS_WSC} ${NUM_PYARTS_WSG}
                   DEPENDS gen_auto_py)
########################################################################################

########### pybind11 class interface
pybind11_add_module(pyarts_cpp THIN_LTO
    py_module.cpp
    py_matpack.cpp
    py_basic.cpp
    py_docserver.cpp
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
    py_physics.cpp
    py_spectroscopy.cpp
    py_jac.cpp
    py_agenda.cpp
    py_workspace.cpp
    py_global.cpp
    py_xsec_fit.cpp
    py_nlte.cpp
    py_constant.cpp
    py_convert.cpp
    py_physics.cpp
    py_predefined.cpp
    py_star.cpp
    py_math.cpp
    py_options.cpp
    py_hitran.cpp
    py_atm.cpp

    py_auto_interface.h
    py_auto_interface.cc
    py_auto_workspace.cc
    ${PYARTS_WORKSPACE_METHODS}
    ${PYARTS_WORKSPACE_VARS}
    ${PYARTS_WORKSPACE_GEN}
    ${PYARTS_WORKSPACE_CREATE}
)

if (ENABLE_PCH)
  target_precompile_headers(pyarts_cpp PRIVATE python_interface.h)
endif()
set_target_properties(pyarts_cpp PROPERTIES OUTPUT_NAME arts)
target_link_libraries(pyarts_cpp PUBLIC artscore)
target_include_directories(pyarts_cpp PUBLIC ${ARTS_SOURCE_DIR}/src)
target_include_directories(pyarts_cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pyarts_cpp PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(pyarts_cpp PUBLIC PYBIND11_DETAILED_ERROR_MESSAGES)
########################################################################################
