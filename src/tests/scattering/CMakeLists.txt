# ########## testcases ###############

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/test_utils.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/test_utils.h
  @ONLY
)

add_executable(test_phase_matrix test_phase_matrix.cc)
target_link_libraries(test_phase_matrix scattering)
target_include_directories(test_phase_matrix PRIVATE ${ARTS_SOURCE_DIR}/src/core ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME "cpp.fast.scattering.test_phase_matrix" COMMAND test_phase_matrix)
add_dependencies(check-deps test_phase_matrix)

add_executable(test_extinction_matrix test_extinction_matrix.cc)
target_link_libraries(test_extinction_matrix scattering)
target_include_directories(test_extinction_matrix PRIVATE ${ARTS_SOURCE_DIR}/src/core ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME "cpp.fast.scattering.test_extinction_matrix" COMMAND test_extinction_matrix)
add_dependencies(check-deps test_extinction_matrix)

add_executable(test_absorption_vector test_absorption_vector.cc)
target_link_libraries(test_absorption_vector scattering)
target_include_directories(test_absorption_vector PRIVATE ${ARTS_SOURCE_DIR}/src/core ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME "cpp.fast.scattering.test_absorption_vector" COMMAND test_absorption_vector)
add_dependencies(check-deps test_absorption_vector)

add_executable(test_sht test_sht.cc)
target_link_libraries(test_sht scattering)
target_include_directories(test_sht PRIVATE ${ARTS_SOURCE_DIR}/src/core ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME "cpp.fast.scattering.test_sht" COMMAND test_sht)
add_dependencies(check-deps test_sht)

add_executable(test_angle_integration test_integration.cc)
target_link_libraries(test_angle_integration scattering)
target_include_directories(test_angle_integration PRIVATE ${ARTS_SOURCE_DIR}/src/core ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME "cpp.fast.scattering.angle_integration" COMMAND test_angle_integration)
add_dependencies(check-deps test_angle_integration)

add_executable(test_single_scattering_data test_single_scattering_data.cc)
target_link_libraries(test_single_scattering_data scattering artsworkspace)
target_include_directories(test_single_scattering_data PRIVATE ${ARTS_SOURCE_DIR}/src/core ${CMAKE_CURRENT_BINARY_DIR})

if(DOWNLOAD_XMLDATA_DEPS)
  add_dependencies(test_single_scattering_data arts_testdata_folder_download_xml)
endif()

ARTS_TEST_RUN_EXE("fast.scattering.test_single_scattering_data" test_single_scattering_data)
add_dependencies(check-deps test_single_scattering_data)

add_executable(test_particle_habit test_particle_habit.cc)
target_link_libraries(test_particle_habit scattering artsworkspace)
target_include_directories(test_particle_habit PRIVATE ${ARTS_SOURCE_DIR}/src/core ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_scattering_habit test_scattering_habit.cc)
target_link_libraries(test_scattering_habit scattering)
target_include_directories(test_scattering_habit PRIVATE ${ARTS_SOURCE_DIR}/src/core ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_psd test_psd.cc)
target_link_libraries(test_psd scattering)
target_include_directories(test_psd PRIVATE ${ARTS_SOURCE_DIR}/src/core ${CMAKE_CURRENT_BINARY_DIR})

if(DOWNLOAD_XMLDATA_DEPS)
  add_dependencies(test_particle_habit arts_testdata_folder_download_xml)
endif()

# add_test(NAME "cpp.fast.scattering.test_particle_habit" COMMAND test_particle_habit)
#add_dependencies(check-deps test_particle_habit)
