# #######################################################################################
if(NOT NUM_PYARTS_WSM)
    set(NUM_PYARTS_WSM 20)
endif()

message(STATUS "Using ${NUM_PYARTS_WSM} python interface automatic workspace method files (-DNUM_PYARTS_WSM=N to change)")

if(NOT NUM_PYARTS_WSV)
    set(NUM_PYARTS_WSV 3)
endif()

message(STATUS "Using ${NUM_PYARTS_WSV} python interface automatic workspace variable files (-DNUM_PYARTS_WSV=N to change)")

set(PYARTS_WORKSPACE_METHODS "")
list(LENGTH PYARTS_WORKSPACE_METHODS N)

while(N LESS NUM_PYARTS_WSM)
    list(APPEND PYARTS_WORKSPACE_METHODS "py_auto_wsm_${N}.cpp")
    list(LENGTH PYARTS_WORKSPACE_METHODS N)
endwhile()

set(PYARTS_WORKSPACE_VARS "")
list(LENGTH PYARTS_WORKSPACE_VARS N)

while(N LESS NUM_PYARTS_WSV)
    list(APPEND PYARTS_WORKSPACE_VARS "py_auto_wsv_${N}.cpp")
    list(LENGTH PYARTS_WORKSPACE_VARS N)
endwhile()

# #######################################################################################
# Build the generator
add_executable(gen_auto_py gen_auto_py.cpp)

set(PYARTS_AUTO_INTERFACE
    py_auto_interface.cpp
)

# Generate the actual files with this command
arts_generate(gen_auto_py
    ARGS ${NUM_PYARTS_WSV} ${NUM_PYARTS_WSM}
    OUTPUT ${PYARTS_AUTO_INTERFACE}
    COMMENT "Generating automatic python interface"
)

add_library(gen_helpers STATIC gen_auto_py_helpers.cpp pydocs.cpp)
target_link_libraries(gen_helpers PUBLIC artsworkspace)
target_include_directories(gen_helpers PUBLIC ${ARTS_SOURCE_DIR}/src)

add_executable(gen_auto_py_methods gen_auto_py_methods.cpp)
target_link_libraries(gen_auto_py_methods PUBLIC artsworkspace gen_helpers)
target_include_directories(gen_auto_py_methods PUBLIC ${ARTS_SOURCE_DIR}/src)

arts_generate(gen_auto_py_methods
    ARGS ${NUM_PYARTS_WSM}
    OUTPUT ${PYARTS_WORKSPACE_METHODS}
    COMMENT "Generating python workspace methods"
)

add_executable(gen_auto_py_options gen_auto_py_options.cpp)
target_link_libraries(gen_auto_py_options PUBLIC artsworkspace gen_helpers)
target_include_directories(gen_auto_py_options PUBLIC ${ARTS_SOURCE_DIR}/src)

set(PYARTS_WORKSPACE_OPTIONS
    py_auto_options.cpp
)

arts_generate(gen_auto_py_options
    OUTPUT ${PYARTS_WORKSPACE_OPTIONS}
    COMMENT "Generating python options"
)

add_executable(gen_auto_py_groups gen_auto_py_groups.cpp)
target_link_libraries(gen_auto_py_groups PUBLIC artsworkspace gen_helpers)
target_include_directories(gen_auto_py_groups PUBLIC ${ARTS_SOURCE_DIR}/src)

set(PYARTS_WORKSPACE_GROUPS
    py_auto_wsg_from_py_object.cpp
    py_auto_wsg_from_const_py_object.cpp
    py_auto_wsg_string_type.cpp
    py_auto_wsg_convert.cpp
    py_auto_wsg_to_py_wsv.cpp
    py_auto_wsg_wsv_implicit.cpp
    py_auto_wsg_wsv_implicit.h
    py_auto_wsg.h
    py_auto_wsgdocs.h
    py_auto_wsgdocs.cpp
    py_auto_agenda_operators.cpp
)

arts_generate(gen_auto_py_groups
    OUTPUT ${PYARTS_WORKSPACE_GROUPS}
    COMMENT "Generating python workspace groups"
)

add_executable(gen_auto_py_variables gen_auto_py_variables.cpp)
target_link_libraries(gen_auto_py_variables PUBLIC artsworkspace gen_helpers)
target_include_directories(gen_auto_py_variables PUBLIC ${ARTS_SOURCE_DIR}/src)

arts_generate(gen_auto_py_variables
    ARGS ${NUM_PYARTS_WSV}
    OUTPUT ${PYARTS_WORKSPACE_VARS}
    COMMENT "Generating python workspace variables"
)

# #######################################################################################

set(PYARTS_AUTO_FILES
    ${PYARTS_AUTO_INTERFACE}
    ${PYARTS_WORKSPACE_METHODS}
    ${PYARTS_WORKSPACE_VARS}
    ${PYARTS_WORKSPACE_GROUPS}
    ${PYARTS_WORKSPACE_OPTIONS}
)

# #######################################################################################

if (IPO_SUPPORTED)
  set(NANOBIND_LTO "LTO")
endif()

# NOTE: Sort these by how long they take to compile, to help multicore compilation
nanobind_add_module(pyarts_cpp NB_STATIC NB_DOMAIN ARTS ${NANOBIND_LTO} STABLE_ABI
    hpy_numpy_impl.cpp

    # py_griddedfield.cpp, py_matpack.cpp and py_rtepack.cpp were previously
    # single monolithic files (and the three slowest, most memory-hungry TUs
    # in the whole build, since nanobind's binding machinery is instantiated
    # independently per bound type). Each is now a small hub plus a few
    # sibling files grouping independent types, so the work compiles in
    # parallel with a much lower peak-memory ceiling per TU.
    py_griddedfield.cpp
    py_griddedfield_named.cpp
    py_griddedfield_stokvec.cpp
    py_griddedfield_extra.cpp
    py_matpack.cpp
    py_matpack_fixed.cpp
    py_matpack_complex.cpp
    py_rtepack.cpp
    py_rtepack_stokvec.cpp
    py_rtepack_propmat.cpp
    py_rtepack_muelmat.cpp
    py_rtepack_specmat.cpp

    ${PYARTS_AUTO_FILES}

    py_workspace.cpp
    py_agenda.cpp

    py_basic.cpp
    py_sensor.cpp
    py_lbl.cpp
    py_disort.cpp
    py_time.cpp
    py_species.cpp
    py_sparse.cpp
    py_scattering.cpp
    py_scattering_species.cpp
    py_psd.cpp
    py_quantum.cpp
    py_physics.cpp
    py_jac.cpp
    py_global.cpp
    py_lookup.cpp
    py_xsec_fit.cpp
    py_constant.cpp
    py_retrieval.cpp
    py_convert.cpp
    py_path.cpp
    py_predefined.cpp
    py_star.cpp
    py_math.cpp
    py_hitran.cpp
    py_atm.cpp
    py_surf.cpp
    py_fwd.cpp
    py_cia.cpp
    py_operators.cpp
    py_interp.cpp
    py_igrf.cpp
    py_zeeman.cpp
    py_module.cpp
    py_file.cpp
    py_nlte.cpp
    py_geodesy.cpp
    py_planets.cpp
    py_montecarlo.cpp
    py_rng.cpp
    py_tmatrix.cpp
)

if(ENABLE_PCH)
    target_precompile_headers(pyarts_cpp PRIVATE python_interface.h)
endif()

set_target_properties(pyarts_cpp PROPERTIES OUTPUT_NAME arts)
target_link_libraries(pyarts_cpp PUBLIC artsworkspace gen_helpers vdisort-cpp)
target_include_directories(pyarts_cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pyarts_cpp PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(
    pyarts_cpp PROPERTIES LIBRARY_OUTPUT_DIRECTORY
    "${ARTS_BINARY_DIR}/python/src/pyarts3/")

if (APPLE
    AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
    AND CMAKE_CXX_STANDARD GREATER_EQUAL 26
    AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 16)
  message(STATUS "WARNING: Linking pyarts_cpp on GCC16 on macOS requires -static-libstdc++ and -static-libgcc to avoid runtime errors. Adding these flags to the link command.  Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")
  target_link_options(pyarts_cpp PRIVATE
      -static-libstdc++ -static-libgcc)
  target_link_options(pyarts_cpp PRIVATE
      -Wl,-reexport_library,/usr/lib/libiconv.2.dylib)
endif()

include (ArtsNanobindStubs)
arts_add_cpp_stubs(
  constants
  convert
  disort
  vdisort
  file
  geodetic
  globals
  hitran
  interp
  lbl
  math
  physics
  predef
  rtepack
  zeeman
  tmatrix
)
