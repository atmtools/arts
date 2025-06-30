add_subdirectory(partfun)
add_subdirectory(python_interface)
add_subdirectory(core)

if(NETCDF_FOUND)
  set(NETCDF_CC_FILES
    nc_io.cc
    nc_io_array_types.cc
    nc_io_basic_types.cc
  )
endif(NETCDF_FOUND)

if(ENABLE_MPI)
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif(ENABLE_MPI)

file(GLOB HEADERFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
add_custom_target(UtilityHeadersArts SOURCES ${HEADERFILES})

# ####
add_library(workspace_core STATIC
  workspace_agendas.cpp
  workspace_groups.cpp
  workspace_variables.cpp
)
target_link_libraries(workspace_core PRIVATE arts_options)

add_executable(make_auto_wsg make_auto_wsg.cpp)
set (ARTS_WORKSPACE_GROUPS
  auto_workspace.cpp
  auto_agenda_operators.cpp
)
set (ARTS_WORKSPACE_INTERNALS_GROUPS
  auto_wsv_value_wrapper_copied.cpp
  auto_wsv_value_wrapper_vformat.cpp
  auto_wsv_value_wrapper_read_from_stream.cpp
  auto_wsv_value_wrapper_write_to_stream.cpp
  auto_wsv_value_wrapper_from_named_type.cpp
  auto_wsv_value_wrapper_type_name.cpp
  auto_wsg_docs.cpp
  auto_wsg_valid_wsg.cpp
)
add_custom_command(OUTPUT
    ${ARTS_WORKSPACE_GROUPS}
    ${ARTS_WORKSPACE_INTERNALS_GROUPS}
    auto_wsg.h
    auto_wsv_value_wrapper.h
    auto_agenda_operators.h
  COMMAND make_auto_wsg
  DEPENDS make_auto_wsg
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating automatic workspace groups"
)
target_link_libraries(make_auto_wsg PRIVATE util workspace_core)


# #######################################################################################

# copy auto generated files if they are different

set(ARTS_WORKSPACE_INTERNALS_FILES "")
foreach(file ${ARTS_WORKSPACE_INTERNALS_GROUPS})
add_custom_command(
    COMMENT "Copying internal workspace file (if required): ${file}"
    OUTPUT x${file}
    DEPENDS ${file}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${file}
        x${file}
    )
  list(APPEND ARTS_WORKSPACE_INTERNALS_FILES x${file})
endforeach()

# ####
add_library(workspace_internals STATIC
  workspace_methods.cpp
  workspace_agenda_class_default_constructor.cpp
  workspace_meta_methods.cpp
  workspace_variables_keywords.cpp
  ${ARTS_WORKSPACE_INTERNALS_FILES}
)
target_link_libraries(workspace_internals PUBLIC workspace_core artscore scattering)
target_include_directories(workspace_internals PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(workspace_internals PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# ####
add_executable(make_auto_wsv make_auto_wsv.cpp)
set (ARTS_WORKSPACE_VARIABLES auto_wsv.cpp)
add_custom_command(OUTPUT
    ${ARTS_WORKSPACE_VARIABLES}
    auto_wsv.h
  COMMAND make_auto_wsv
  DEPENDS make_auto_wsv
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating automatic workspace variables"
)
target_link_libraries(make_auto_wsv PRIVATE workspace_internals)

# ####
if(NOT NUM_ARTS_WSM)
  set(NUM_ARTS_WSM 10)
endif()

message(STATUS "Using ${NUM_ARTS_WSM} automatic workspace method files (--DNUM_ARTS_WSM=N to change)")

set(ARTS_WORKSPACE_METHODS "")
list(LENGTH ARTS_WORKSPACE_METHODS N)

while(N LESS NUM_ARTS_WSM)
  list(APPEND ARTS_WORKSPACE_METHODS "auto_wsm_${N}.cpp")
  list(LENGTH ARTS_WORKSPACE_METHODS N)
endwhile()

list(APPEND ARTS_WORKSPACE_METHODS auto_wsm.cpp auto_wsmmeta.cpp)

add_executable(make_auto_wsm make_auto_wsm.cpp)
add_custom_command(OUTPUT
    ${ARTS_WORKSPACE_METHODS}
    auto_wsm.h
  COMMAND make_auto_wsm ${NUM_ARTS_WSM}
  DEPENDS make_auto_wsm
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating automatic workspace methods"
)
target_link_libraries(make_auto_wsm PRIVATE workspace_internals)

# ####

set (ARTS_WORKSPACE_AGENDAS auto_wsa.cpp)

add_executable(make_auto_wsa make_auto_wsa.cpp)
add_custom_command(OUTPUT
    ${ARTS_WORKSPACE_AGENDAS}
    auto_wsa.h
    auto_wsa_options.h
  COMMAND make_auto_wsa
  DEPENDS make_auto_wsa
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating automatic workspace agendas"
)
target_link_libraries(make_auto_wsa PRIVATE workspace_internals)

# #######################################################################################

# copy auto generated files if they are different

set(ARTS_WORKSPACE_FILES "")
foreach(file ${ARTS_WORKSPACE_METHODS}
             ${ARTS_WORKSPACE_VARIABLES}
             ${ARTS_WORKSPACE_GROUPS}
             ${ARTS_WORKSPACE_AGENDAS})
add_custom_command(
    COMMENT "Copying workspace file (if required): ${file}"
    OUTPUT x${file}
    DEPENDS ${file}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${file}
        x${file}
    )
  list(APPEND ARTS_WORKSPACE_FILES x${file})
endforeach()

# ########## next target ###############
add_library(artsworkspace STATIC
  ${ARTS_WORKSPACE_FILES}
  callback.cpp
  workspace_agenda_class.cpp
  workspace_agenda_creator.cpp
  workspace_method_class.cpp
  workspace_method_extra_doc.cpp
  workspace_class.cpp
  m_abs.cc
  m_atm.cc
  m_atm_profile.cc
  m_background.cc
  m_cat.cc
  m_cia.cc
  m_covmat.cc
  m_disort.cc
  m_disort_settings.cc
  m_frequency_grid.cc
  m_fwd.cc
  m_predefined_absorption_models.cc
  m_xsec_fit.cc
  m_jactargets.cc
  m_lbl.cc
  m_linemixing.cc
  m_lookup.cc
  m_measurement_vector.cc
  m_model_state.cc
  m_nc.cc
  m_nlte.cc
  m_obsel.cc
  m_oem.cc
  m_operators.cc
  m_partfun.cc
  m_planets.cc
  m_ppvar.cc
  m_propagation_path.cc
  m_propagation_path_observer.cc
  m_propmat.cc
  m_retrieval.cc
  m_rad.cc
  m_scattering_species.cc
  m_spectral_flux.cc
  m_spectral_radiance.cc
  m_spectral_radiance_field.cc
  m_sun.cc
  m_surf.cc
  m_surface_spectral_radiance.cc
  m_wigner.cc
  m_xml.cc
  sun_methods.cc
  version.cc
  ${NETCDF_CC_FILES}
)

target_include_directories(artsworkspace PUBLIC
  ${ARTS_SOURCE_DIR}/3rdparty/invlib/src
  ${NETCDF_INCLUDE_DIR}
)

if(ENABLE_PCH)
  target_precompile_headers(artsworkspace PRIVATE workspace.h)
endif()

target_link_libraries(artsworkspace PUBLIC auto_partfun artscore workspace_internals)

if(NETCDF_FOUND)
  target_link_libraries(artsworkspace PUBLIC ${NETCDF_LIBRARIES})
endif(NETCDF_FOUND)

if(ENABLE_MPI)
  target_link_libraries(artsworkspace PUBLIC ${MPI_CXX_LIBRARIES})
endif(ENABLE_MPI)

if(ENABLE_RT4)
  target_link_libraries(artsworkspace PUBLIC rt4)
endif(ENABLE_RT4)

if(ENABLE_FASTEM)
  target_link_libraries(artsworkspace PUBLIC fastem)
endif(ENABLE_FASTEM)

if(ENABLE_TMATRIX)
  target_link_libraries(artsworkspace PUBLIC tmatrix)
endif(ENABLE_TMATRIX)

# ########## next target ###############


# #######################################################################################

# ########## IPO / LTO support
if(IPO_SUPPORTED)
  set_property(TARGET artsworkspace PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET binio PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET species PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_subdirectory(tests)
