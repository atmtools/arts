add_subdirectory(partfun)
add_subdirectory(python_interface)
add_subdirectory(core)

if(NETCDF_FOUND)
  set(NETCDF_CC_FILES
    nc_io.cc
    nc_io_array_types.cc
    nc_io_basic_types.cc
    nc_io_compound_types.cc
  )
endif(NETCDF_FOUND)

if(ENABLE_MPI)
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif(ENABLE_MPI)

file(GLOB HEADERFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
add_custom_target(UtilityHeadersArts SOURCES ${HEADERFILES})

# ########## target to extract version number from ChangeLog ###############
add_custom_target(auto_version_cc DEPENDS UpdateAutoVersion)

add_custom_command(
  OUTPUT auto_version.cc
  COMMAND ${CMAKE_COMMAND} -D "ARTS_BINARY_DIR:STRING=\"${ARTS_BINARY_DIR}\"" -D "ARTS_SOURCE_DIR:STRING=\"${ARTS_SOURCE_DIR}\"" -P ${ARTS_SOURCE_DIR}/cmake/scripts/update_auto_version_cc.cmake
  BYPRODUCTS auto_version_cc.txt
  DEPENDS ${ARTS_BINARY_DIR}/auto_version.txt
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Updating ARTS version number"
)

# ####
add_library(workspace_group STATIC
  workspace_groups.cpp
)

add_executable(make_auto_wsg make_auto_wsg.cpp)
add_custom_command(OUTPUT
  auto_wsg.h
  auto_wsg_init.cpp
  auto_wsg_share.cpp
  auto_workspace.cpp
  COMMAND make_auto_wsg
  DEPENDS make_auto_wsg
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating automatic workspace groups"
)
target_link_libraries(make_auto_wsg PRIVATE workspace_group)

# ####
add_library(workspace_internals STATIC
  workspace_agenda_class_default_constructor.cpp
  workspace_agendas.cpp
  workspace_meta_methods.cpp
  workspace_methods.cpp
  workspace_variables.cpp
  workspace_variables_keywords.cpp
  workspace_wsv.cpp
  auto_wsg_init.cpp
  auto_wsg_share.cpp
)
target_link_libraries(workspace_internals PUBLIC workspace_group artscore scattering)
target_include_directories(workspace_internals PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(workspace_internals PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# ####
add_executable(make_auto_wsv make_auto_wsv.cpp)
add_custom_command(OUTPUT auto_wsv.h auto_wsv.cpp
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

add_executable(make_auto_wsm make_auto_wsm.cpp auto_wsv.cpp)
add_custom_command(OUTPUT auto_wsm.h auto_wsm.cpp auto_wsmmeta.cpp ${ARTS_WORKSPACE_METHODS}
  COMMAND make_auto_wsm ${NUM_ARTS_WSM}
  DEPENDS make_auto_wsm
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating automatic workspace methods"
)
target_link_libraries(make_auto_wsm PRIVATE workspace_internals)

# ####
add_executable(make_auto_wsa make_auto_wsa.cpp)
add_custom_command(OUTPUT auto_wsa.h auto_wsa.cpp
  COMMAND make_auto_wsa
  DEPENDS make_auto_wsa
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating automatic workspace agendas"
)
target_link_libraries(make_auto_wsa PRIVATE workspace_internals)

# ########## next target ###############
add_library(artsworkspace STATIC
  auto_wsv.cpp
  auto_wsm.cpp
  auto_wsmmeta.cpp
  auto_wsa.cpp
  auto_workspace.cpp
  ${ARTS_WORKSPACE_METHODS}
  callback.cpp
  workspace_agenda_class.cpp
  workspace_agenda_creator.cpp
  workspace_method_class.cpp
  workspace_class.cpp
  auto_version.cc
  compile_time_tests.cc
  m_abs.cc
  m_absorptionlines.cc
  m_agenda_set.cc
  m_atm.cc
  m_background.cc
  m_cat.cc
  m_cia.cc
  m_covmat.cc
  m_disort.cc
  m_disort_settings.cc
  m_fwd.cc
  m_predefined_absorption_models.cc
  m_xsec_fit.cc
  m_jactargets.cc
  m_lbl.cc
  m_linemixing.cc
  m_model_state.cc
  m_nc.cc
  m_obsel.cc
  m_oem.cc
  m_operators.cc
  m_partfun.cc
  m_planets.cc
  m_poslos.cc
  m_ppvar.cc
  m_propagation_path.cc
  m_retrieval.cc
  m_rad.cc
  m_spectral_radiance.cc
  m_sun.cc
  m_surf.cc
  m_wigner.cc
  m_xml.cc
  sun_methods.cc
  xml_io.cc
  xml_io_array_types.cc
  xml_io_basic_types.cc
  xml_io_compound_types.cc
  xml_io_rtepack.cc
  xml_operators.cc
  ${NETCDF_CC_FILES}
)

target_include_directories(artsworkspace PUBLIC
  ${ARTS_SOURCE_DIR}/3rdparty/invlib/src
  ${NETCDF_INCLUDE_DIR}
)

add_dependencies(artsworkspace auto_version_cc)

if(ENABLE_PCH)
  target_precompile_headers(artsworkspace PRIVATE workspace.h)
endif()

target_link_libraries(artsworkspace PUBLIC auto_partfun xmliobase artscore workspace_internals)

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
add_library(xmliobase STATIC
  xml_io_base.cc
  xml_io_general_types.cc
)

target_link_libraries(xmliobase PUBLIC binio matpack)

# #######################################################################################

# ########## IPO / LTO support
if(IPO_SUPPORTED)
  set_property(TARGET artsworkspace PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET binio PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET species PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET xmliobase PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_subdirectory(tests)
