add_subdirectory(gui)
add_subdirectory(partfun)
add_subdirectory(python_interface)
add_subdirectory(workspace)
add_subdirectory(core)

if (NETCDF_FOUND)
  set (NETCDF_CC_FILES
    nc_io.cc
    nc_io_array_types.cc
    nc_io_basic_types.cc
    nc_io_compound_types.cc
    )
endif (NETCDF_FOUND)

if (ENABLE_MPI)
  include_directories (${MPI_CXX_INCLUDE_PATH})
endif (ENABLE_MPI)

file (GLOB HEADERFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
add_custom_target(UtilityHeadersArts SOURCES ${HEADERFILES})

########### target to extract version number from ChangeLog ###############

add_custom_target (auto_version_cc DEPENDS UpdateAutoVersion)

add_custom_command (
  OUTPUT auto_version.cc
  COMMAND ${CMAKE_COMMAND} -D "ARTS_BINARY_DIR:STRING=\"${ARTS_BINARY_DIR}\"" -D "ARTS_SOURCE_DIR:STRING=\"${ARTS_SOURCE_DIR}\"" -P ${ARTS_SOURCE_DIR}/cmake/scripts/update_auto_version_cc.cmake
  DEPENDS ${ARTS_BINARY_DIR}/auto_version.txt
)

########### next target ###############

add_library (artscore STATIC
  arts_omp.cc
  auto_version.cc
  compile_time_tests.cc
  disort.cc
  doit.cc
  gui_propmat.cc
  m_abs.cc
  m_absorptionlines.cc
  m_abs_lookup.cc
  m_agenda_set.cc
  m_append.cc
  m_array.cc
  m_artstime.cc
  m_atmosphere.cc
  m_atm.cc
  m_basic_types.cc
  m_batch.cc
  m_checked.cc
  m_cia.cc
  m_cloudbox.cc
  m_cloudradar.cc
  m_conversion.cc
  m_doit.cc
  m_extract.cc
  m_fluxes.cc
  m_fwd.cc
  m_gridded_fields.cc
  m_predefined_absorption_models.cc
  m_propmat.cc
  m_gas_scattering.cc
  m_general.cc
  m_geodetic.cc
  m_xsec_fit.cc
  m_jacobian.cc
  m_linemixing.cc
  m_microphysics.cc
  m_montecarlo.cc
  m_nc.cc
  m_nlte.cc
  m_oem.cc
  m_operators.cc
  m_optproperties.cc
  m_partfun.cc
  m_physics.cc
  m_planets.cc
  m_poslos.cc
  m_ppath.cc
  m_psd.cc
  m_raw.cc
  m_reduce.cc
  m_refraction.cc
  m_retrieval.cc
  m_rt4.cc
  m_rte.cc
  m_retrieval.cc
  m_select.cc
  m_sensor.cc
  m_sizes.cc
  m_sun.cc
  m_surf.cc
  m_surface.cc
  m_telsem.cc
  m_tessem.cc
  m_tmatrix.cc
  m_transmitter.cc
  m_wigner.cc
  m_xml.cc
  m_zeeman.cc
  montecarlo.cc
  ppath.cc
  propmat_field.cc
  psd.cc
  rt4.cc
  rte.cc
  sourcetext.cc
  sun_methods.cc
  surface.cc
  variousZZZ.cc
  xml_io.cc
  xml_io_array_types.cc
  xml_io_basic_types.cc
  xml_io_compound_types.cc
  xml_io_rtepack.cc
  ${NETCDF_CC_FILES}
        m_fluxes.cc)

target_include_directories (artscore PUBLIC 
  ${ARTS_SOURCE_DIR}/3rdparty/cdisort
  ${ARTS_SOURCE_DIR}/3rdparty/invlib/src
  ${NETCDF_INCLUDE_DIR}
  )

add_dependencies (artscore auto_version_cc)
if (ENABLE_PCH)
  target_precompile_headers(artscore PRIVATE pch_artscore.h)
endif()

target_link_libraries(artscore PUBLIC auto_partfun cdisort xmliobase artsfun workspace)

if (ENABLE_GUI)
  target_link_libraries (artscore PUBLIC artsgui)
  target_compile_definitions(artscore PUBLIC ARTS_GUI_ENABLED)
endif()

if (NETCDF_FOUND)
  target_link_libraries (artscore PUBLIC ${NETCDF_LIBRARIES})
endif (NETCDF_FOUND)

if (ENABLE_MPI)
  target_link_libraries (artscore PUBLIC ${MPI_CXX_LIBRARIES})
endif (ENABLE_MPI)

if (OpenMP_CXX_FOUND)
  target_link_libraries (artscore PUBLIC OpenMP::OpenMP_CXX)
endif ()

if (ENABLE_RT4)
  target_link_libraries (artscore PUBLIC rt4)
endif (ENABLE_RT4)

if (ENABLE_FASTEM)
  target_link_libraries (artscore PUBLIC fastem)
endif (ENABLE_FASTEM)

if (ENABLE_TMATRIX)
  target_link_libraries (artscore PUBLIC tmatrix)
endif (ENABLE_TMATRIX)

########### next target ###############

target_link_libraries(artscore PUBLIC rtepack)

########### next target ###############

add_library (xmliobase STATIC
        xml_io_base.cc
        xml_io_general_types.cc
        )

target_link_libraries (xmliobase PUBLIC binio matpack)

########################################################################################

########### IPO / LTO support
if (IPO_SUPPORTED)
  set_property(TARGET arts PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET artscore PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET binio PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET species PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET xmliobase PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_subdirectory(gpl)
target_link_libraries(artscore PUBLIC artscore_gpl_selector)

add_subdirectory(tests)

target_link_libraries(artscore PUBLIC fwd)
